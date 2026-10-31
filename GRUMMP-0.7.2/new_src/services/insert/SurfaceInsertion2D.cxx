#include "GR_SurfaceInsertion.h"

#include <cmath>

#include "GR_LengthAniso2D.h"
#include "GR_Mesh2D.h"
#include "GR_Geometry.h"
#include "GR_GRCurve.h"
#include "GR_GRCurveGeom.h"

namespace GRUMMP {
/// SurfaceInserter2D begins here

/// Constructor, sets up all objects, calls base constructor
SurfaceInserter2D::SurfaceInserter2D(Mesh2D* mesh2D,
		std::shared_ptr<GRUMMP::Length> len, const double scale,
		const double grading, const char strSurfParams[]) :
		SurfaceInserter(len, grading, scale, strSurfParams), m_mesh(mesh2D), m_geom(
		NULL), m_pointsToInsert(), m_offsetCurveMap(), SwapDec2D(NULL), SmoothMan2D(
				m_mesh), RefMan2D(NULL) {
	m_anisoLength = new LengthAniso2D(LengthAniso2D::eMesh, 1, 1, LARGE_DBL);
	m_anisoLength->provideMesh(m_mesh);

	GR_index_t iLen = strlen(strSurfParams), iParam = 0;
	while (iParam < iLen) {
		char cType = strSurfParams[iParam++];
		switch (cType) {
		// swap decider type, right now only one type...
		case 'd': {
			char sdType = strSurfParams[iParam++];
			switch (sdType) {
			case 'd':
				SwapDec2D = new DelaunaySwapDecider2D();
				break;
			}
			break;
		}
			// aniso
		case 'a': {
			logMessage(1, "Using anisotropic insertion\n");
			m_isAnisotropic = true;
			// stems from the way AnisoRefinement currently done

			ARMesh2D = new AnisoRefinement(m_mesh,
					dynamic_cast<LengthAniso2D*>(m_anisoLength), 2, 0);
			m_mesh->attachAnisoRef(ARMesh2D);

			//setVertLengthScaleMetrics();
			break;
		}
			// o and p are offset handling
		case 'o': {
			m_useOffsets = true;
			char temp1, temp2, temp3;
			int curvenum;
			double x, y;
			sscanf(&strSurfParams[iParam], "%d%c%lf%c%lf%c%d", &curvenum,
					&temp1, &x, &temp2, &y, &temp3, &m_numberOfOffsets);
			m_offsetCoordinate = std::pair<int, CubitVector>(curvenum,
					CubitVector(x, y, 0));
			if (temp1 != 'x' || temp2 != 'y' || temp3 != 'n') {
				logMessage(1, "Warning, incorrect offset input string\n");
			}
			if (curvenum < 0) {
				m_offsetByNormals = true;
			}
		}
			break;
		case 'p': {

			m_useOffsets = true;
			char * temp = const_cast<char*>(&strSurfParams[iParam]);
			m_numberOfOffsets = strtoul(temp, &temp, 10);
			iParam += temp - &strSurfParams[iParam];
			char * pairs = strtok(temp, "c");
			for (GR_index_t iP = 0; iP < m_numberOfOffsets; iP++) {
				// collect pairs
				int c1 = strtol(pairs, &temp, 10);
				pairs = strtok(NULL, "c");
				int c2 = strtol(pairs, &temp, 10);
				pairs = strtok(NULL, "c");
				m_offsetPairs.push_back(std::pair<int, int>(c1, c2));
			}
			m_offsetByNormals = true;
		}
			break;
		case 'q': // quasi-structured
		{
			m_isStructured = true;
		}
			break;
		}
	}

	// safeguard incase we don't set it
	if (SwapDec2D == NULL) {
		SwapDec2D = new DelaunaySwapDecider2D();
	}

	SwapMan2D = new SwapManager2D(SwapDec2D, m_mesh);

	if (!m_isAnisotropic) {
		RefMan2D = new RefinementManager2D(m_mesh, m_IPC, m_isoLength, NULL);
	}
	// this class only watches when vertices are moved, for stat tracking
	m_mesh->addObserver(this, Observable::vertMoved);

	initializeStats(m_mesh);

}

SurfaceInserter2D::~SurfaceInserter2D() {

	finalizeStats(*m_mesh);

	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); cpm_it++) {
		cpm_it->second->clear();
		delete cpm_it->second;
	}
	m_mesh->removeObserver(this);
	if (SwapMan2D) {
		delete SwapMan2D;
	}
	if (SwapDec2D) {
		delete SwapDec2D;
	}
	if (m_isAnisotropic) {
		delete ARMesh2D;
	} else {
		delete RefMan2D;
	}

}

void SurfaceInserter2D::insertCurveInMesh(const char strCurveFilename[],
		const char strParams[]) {
	bool wasSimplicial = m_mesh->isSimplicial();
	if (!wasSimplicial)
		m_mesh->makeSimplicial();

	m_geom = new GRGeom2D(strCurveFilename);

#ifndef NDEBUG
	// writes the geometry to VTK file to check
	m_geom->writeGeometryVTK("surfins_geom");
#endif

	logMessage(1, "Reading Geometry File: %s \n", strCurveFilename);
	// Samples the surface
	(!m_isStructured) ? sampleSurface() : sampleQuasiStructured();

	// Removes vertices that encroach on the surface
	removeVertsNearSurface();

	writeVTKLegacyWithoutPurge(*m_mesh, "surfins_removemesh");

	// Insert sampled points into mesh
	insertPointsIntoMesh();

	writeVTKLegacyWithoutPurge(*m_mesh, "surfins_newinsert");

	assert(m_mesh->isValid());

	// re-numbers the new region(s)
	setRegions();

	// post-processing to ensure mesh is of adequate quality
	refineAndSmooth(strParams);

	// Re-combine tris into quads; there's currently no data available
	// to guide this, since the layer data that the advancing quad code
	// uses doesn't currently exist.

	assert(m_mesh->isValid());

	delete m_geom;
}
void SurfaceInserter2D::refineBoundaryForLength(GRGeom2D* geom) {
	m_geom = geom;
	m_isAnisotropic = false;
	for (GR_index_t i = 0; i < m_mesh->getNumTotalBdryFaces(); ++i) {
		BdryEdgeBase* bdryEdge =
				dynamic_cast<BdryEdgeBase*>(m_mesh->getBFace(i));
		double startParam = bdryEdge->getVert0Param();
		double endParam = bdryEdge->getVert1Param();
		double minParam = bdryEdge->getCurve()->get_curve_geom()->min_param();
		double maxParam = bdryEdge->getCurve()->get_curve_geom()->max_param();
		if (std::fabs(startParam - minParam) > 1e-13
				&& std::fabs(startParam - maxParam) > 1e-13)
			m_leftoverVertsToRemove.insert(bdryEdge->getVert(0));
		if (std::fabs(endParam - minParam) > 1e-13
				&& std::fabs(endParam - maxParam) > 1e-13)
			m_leftoverVertsToRemove.insert(bdryEdge->getVert(1));

	}
	sampleSurfaceForBoundaryPoints();
	writeSamplePointsToFile("testpoints");
	insertBoundaryPointsIntoMesh();
	removeVertsAfterInsertion();

}

static double boundaryLengthScale(GRCurveGeom * curve, double param,
		double scale, double grading) {
	if (curve->get_curve_type() == GRCurveGeom::ARC) {
		double length = curve->arc_length();
		double lengthscale = length / scale * 2.0;
		assert(lengthscale > 0);
		return lengthscale;
	} else {
		double length = curve->arc_length();

		double curvature = 1. / curve->curvature(param) / scale;
		CubitVector coord;
		curve->coord_at_param(param, coord);
		double lengthscale = length / scale;
		double lengthAtParam = curve->arc_length(curve->min_param(), param);
		double endscale = length * 0.25;
		double lengthFromEnd = std::min(length - lengthAtParam, lengthAtParam);
		if (lengthFromEnd < endscale) {
//			lengthscale = 0.5*(length/scale-length/scale/grading)*(sin(M_PI*lengthFromEnd/endscale-M_PI_2)+1.)+length/scale/grading;
//			lengthscale = (length/scale-length/scale/grading)*lengthFromEnd/endscale+length/scale/grading;
			lengthscale = (length / scale - length / scale / grading)
					* (sin(lengthFromEnd / endscale) * M_PI_2)
					+ length / scale / grading;
		}
		/////

		double minCurvature = curve->getMinRadiusOfCurvature() / scale;
		minCurvature = std::min(curvature, minCurvature);
		assert(iFuzzyComp(curvature, minCurvature) >= 0);
		if (curvature < lengthscale) {
//			lengthscale = std::max(curvature,length/scale/grading);
//			lengthscale = (lengthscale - length/scale/grading)*(curvature - minCurvature)/(lengthscale - minCurvature) + length/scale/grading;
			lengthscale = (lengthscale - length / scale / grading)
					* sin(
							(curvature - minCurvature)
									/ (lengthscale - minCurvature) * M_PI_2)
					+ length / scale / grading;
//			lengthscale = (lengthscale - length/scale/grading)*0.5*(sin((curvature - minCurvature)/(lengthscale - minCurvature)*M_PI-M_PI_2)+1.) + length/scale/grading;
		}
		assert(iFuzzyComp(lengthscale, length / scale / grading) >= 0);

		assert(lengthscale > 0);
		return lengthscale;
	}
}

void SurfaceInserter2D::equidistributeArcLength(
		std::vector<double>& paramToInsert,
		const std::vector<double>& newPointsArclength, GRCurveGeom* curveGeom,
		insPointSet& pointSet) {
	// otherwise, we need to spread points out and optimize point locations
	// number of sampling points for this curve.
	GR_index_t numPoints = paramToInsert.size();
	// lengthscale at the midpoint of the edges
	double* midpointLS = new double[numPoints - 1];
	// ntemp arclength for our jacobi iteration
	double* tempArclength = new double[numPoints];
	// initialize our arclengths
	for (GR_index_t index = 0; index < numPoints; index++) {
		tempArclength[index] = newPointsArclength[index];
	}
	// initialize our lengthscales at the midpoints of arcs
	for (GR_index_t index = 0; index < numPoints - 1; index++) {
		midpointLS[index] = boundaryLengthScale(curveGeom,
				curveGeom->param_at_arc_length(paramToInsert[index],
						0.5
								* (newPointsArclength[index + 1]
										- newPointsArclength[index])), m_scale,
				m_grading);
	}
	// lets smooth internal points using equidistribution
	// of length scale
	// uses jacobi on arclength
	double oldArclength = 0.;
	GR_index_t iter = 0;
	// initialize this, to keep track of convergence
	double change_in_arclen = LARGE_DBL;
	//
	while (change_in_arclen > 1.e-2 && ++iter < 1000) {
		change_in_arclen = 0.0;
		for (GR_index_t index = 1; index < numPoints - 1; index++) {
			oldArclength = tempArclength[index]; // used for determining change
			// smoothing
			tempArclength[index] = (midpointLS[index - 1]
					* tempArclength[index + 1]
					+ midpointLS[index] * tempArclength[index - 1])
					/ (midpointLS[index] + midpointLS[index - 1]);
			// compute parameter at that arclength
			paramToInsert[index] = curveGeom->param_at_arc_length(
					paramToInsert[index - 1],
					(tempArclength[index] - tempArclength[index - 1]));
			// compute the lengthscale at the new location
			// TODO: Should this be tempAL instead of newPointsAL?  The latter seems like
			// it wouldn't account for changes in the position of the points.  Not likely
			// to make a very big difference, but still, it would nice if it were right.
			midpointLS[index - 1] = boundaryLengthScale(curveGeom,
					curveGeom->param_at_arc_length(paramToInsert[index - 1],
							0.5
									* (newPointsArclength[index]
											- newPointsArclength[index - 1])),
					m_scale, m_grading);
			// update the total change
			change_in_arclen += fabs(oldArclength - tempArclength[index])
					/ oldArclength;
		}
		logMessage(MSG_SERVICE, "iter: %d \t change_in_arclen: %6.4f out of %6.4f\n", iter,
				change_in_arclen, oldArclength);
	}
	// inserts all the points based on the computed values.
	for (GR_index_t index = 1; index < numPoints - 1; index++) {
		double currentParam = curveGeom->param_at_arc_length(0.0,
				tempArclength[index]);
		pointSet.insert(insPoint2D(currentParam, NULL));
	}
	// clean up the temporary arrays
	delete[] midpointLS;
	delete[] tempArclength;
}

int SurfaceInserter2D::sampleSurfaceForBoundaryPoints() {

	/// get all the curves to sample
	std::list<GRCurve*> curve_list;
	std::list<GRCurve*>::iterator itc;
	m_geom->get_curves(curve_list);

	logMessage(1, "Initializing Sample Locations\n");

	// use of CubitVectors due to original geometry implementation
	CubitVector tempCubitVector;
	int iCurveNum = 0;
	/// sample each curve individually
	for (itc = curve_list.begin(); itc != curve_list.end();
			itc++, iCurveNum++) {

		GRCurve * curveToInsert = (*itc);

		// create the pointset
		insPointSet * pointSet = new insPointSet();

		GRCurveGeom * curveGeom = curveToInsert->get_curve_geom();

		assert(curveGeom);

		// store the arclength at the new points, and their parameters
		std::vector<double> newPointsArclength;
		std::vector<double> paramToInsert;

		double arcLength = curveGeom->arc_length();
		// parameters
		double minParam = curveGeom->min_param();
		double maxParam = curveGeom->max_param();

		// corresponding length scales
		double length = 0.;

		// set the tolerance of distance functions.
		curveGeom->set_tolerance(1.e-6);

		// get start and end lengthscales
		double startLS = boundaryLengthScale(curveGeom, minParam, m_scale, m_grading);
		double endLS = boundaryLengthScale(curveGeom, maxParam, m_scale, m_grading);

		// discrete length between coordinates, used as a comparison
		CubitVector cv1, cv2;
		for (GR_index_t iParam = minParam; iParam < maxParam; iParam++) {
			curveGeom->get_coord(iParam, cv1);
			curveGeom->get_coord(iParam + 1, cv2);
			length += (cv1 - cv2).length();
		}

		if (startLS > 0.95 * length) {
			m_pointsToInsert.insert(
					std::pair<GRCurve*, insPointSet*>(curveToInsert, pointSet));

			return 0;
		}

		paramToInsert.push_back(0.0);
		newPointsArclength.push_back(0.0);

		double currentParam = 0;
		// if the start length is less than the arc length at the parameter, insert a point
		// or if its a straight line, due to curvature issues
		if (iFuzzyComp(startLS, length) <= 0
				|| curveGeom->get_curve_type() == GRCurveGeom::LINE) {
			currentParam = curveGeom->coord_at_dist(startLS, true,
					tempCubitVector);
		} else {
			currentParam = maxParam;
		}

		// while we keep adding initial sampling points

		// if the currentLength at our initial point + the next length to be
		// added is less than the total length - half the endpoint length scale
		// keep going, otherwise stop

		// I added an outer loop here to try to more accurately produce a
		// discretization that is smooth and covers the curve exactly.  In
		// the end, that wasn't wildly successful; what sped the code up
		// dramatically was fixing an N^2 implementation of getting length
		// scale for curved surfaces.  I'm leaving the outer loop in here
		// just in case we may someday hit other geometries for which it has
		// a benefit.
		double scale = 1;
		double tol = 1.e-6;
		double error = 1;
		int passes = 0;
		while (fabs(error) > tol && passes < 6) {
			passes++;
			startLS = scale * boundaryLengthScale(curveGeom, minParam, m_scale, m_grading);
			endLS = scale * boundaryLengthScale(curveGeom, maxParam, m_scale, m_grading);
			double nextLS = startLS;
			currentParam = minParam;
			newPointsArclength.clear();
			paramToInsert.clear();
			paramToInsert.push_back(0.0);
			newPointsArclength.push_back(0.0);

			double currentLength = 0;
			while (currentLength + nextLS < arcLength - 0.5 * endLS) {
				// updates the currParam
				// This is a predictor-corrector RK2-like scheme, which should give us
				// a better initial guess at point distribution than the old
				// code (which was essentially explicit Euler-like).
				double predictedParam = curveGeom->param_at_arc_length(currentParam,
						nextLS);
				double LSAtPredParam = scale
						* boundaryLengthScale(curveGeom, predictedParam, m_scale,
								m_grading);

				currentParam = curveGeom->param_at_arc_length(currentParam,
						(nextLS + LSAtPredParam) / 2);
				// updates cvToInsert
				// recomputes the arc_length
				currentLength = curveGeom->arc_length(0, currentParam);

				paramToInsert.push_back(currentParam);
				newPointsArclength.push_back(currentLength);
//				paramToInsert2.push_back(maxParam - currentParam);
//				newPointsArclength2.push_back(arcLength - currentLength);
				// update the next lengthscale
				nextLS = scale
						* boundaryLengthScale(curveGeom, currentParam, m_scale,
								m_grading);
				logMessage(MSG_DEBUG, "currentPar: %6.8f nextLS: %6.4g\n", currentParam, nextLS);
			}
			double totalLengthUsed = currentLength + endLS;
			error = 1 - totalLengthUsed / arcLength;
			// Predict the next param value indirectly: find the parameter change
			// with half the next step, double, and add to the current parameter.
			// This avoids problems with slight overruns of the max param for
			// curves.
			double predictedParam = 2 * (curveGeom->param_at_arc_length(currentParam,
					nextLS/2) - currentParam) + currentParam;
			logMessage(MSG_DETAIL, "Last param:  %.8f  Pred param: %.8f  Max param: %.8f\n",
					currentParam, predictedParam, maxParam);
			scale = scale * maxParam / predictedParam;

			logMessage(MSG_DETAIL,
					"Total AL: %.8f.  AL used: %.8f. Scale: %.8f.  Error: %8.4G\n",
					arcLength, totalLengthUsed, scale, error);
		}
//		printf("INDEX1 = %d\n",index);
//		nextLS = startLS;
//		currentParam = minParam;
//		currentLength = 0.;

//////		index = 0;
////		while(currentLength+nextLS < 0.5*arcLength){
//////			// updates the currParam
////			currentParam = curveGeom->param_at_arc_length(currentParam, nextLS);
//////			// updates cvToInsert
//////			// recomputes the arc_length
////			currentLength = curveGeom->arc_length(0,currentParam);
//////
//////
//////			// update the next lengthscale
////			nextLS = boundaryLengthScale(curveGeom,currentParam, m_scale, m_grading);
////////			index++;
////		}
////		printf("INDEX2 = %d\n",index);
//		std::vector<int>::size_type sz = paramToInsert2.size();
//		std::vector<int>::size_type sz2 = newPointsArclength2.size();
////		paramToInsert.push_back(maxParam/2.);
////		newPointsArclength.push_back(arcLength/2.);
//		assert(sz == sz2);
//		for (int i = sz - 1; i >= 0; i--) {
//			paramToInsert.push_back(paramToInsert2[i]);
//			newPointsArclength.push_back(newPointsArclength2[i]);
//		}

		paramToInsert.push_back(maxParam);
		newPointsArclength.push_back(arcLength);
		// if there are no points added, we are done
		if (paramToInsert.size() == 2) {
			pointSet->insert(insPoint2D(minParam, NULL));
			m_pointsToInsert.insert(
					std::pair<GRCurve*, insPointSet*>(curveToInsert, pointSet));
		} else {
			equidistributeArcLength(paramToInsert, newPointsArclength,
					curveGeom, *pointSet);
		}
		m_pointsToInsert.insert(
				std::pair<GRCurve*, insPointSet*>(curveToInsert, pointSet));
	}
	// brief output, return the total number of points to insert
	// well aware if a point shares multiple curves, it is counted twice
	GR_index_t np = 0;

	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); cpm_it++) {
		np += cpm_it->second->size();
	}
	logMessage(1, "%u Point Locations Finalized\n", np);
	return np;
}

int SurfaceInserter2D::sampleSurface() {

	/// gets lengthscale inside any cell on the domain
	BaryLengthScale2D BLS(m_mesh, m_isAnisotropic);

	/// get all the curves to sample
	std::list<GRCurve*> curve_list;
	std::list<GRCurve*>::iterator itc;
	m_geom->get_curves(curve_list);

	logMessage(1, "Initializing Sample Locations\n");

	std::set<int> curveNumsToIgnore;
	// should be empty anyways, so don't need to check if I'm using offsetsth
	for (GR_index_t iP = 0; iP < m_offsetPairs.size(); iP++) {
		curveNumsToIgnore.insert(m_offsetPairs[iP].second);
	}

	// use of CubitVectors due to original geometry implementation
	CubitVector tempCubitVector;
	int iCurveNum = 0;
	/// sample each curve individually
	for (itc = curve_list.begin(); itc != curve_list.end();
			itc++, iCurveNum++) {

		// skip the curve if we don't need to sample;
		if (m_useOffsets && curveNumsToIgnore.count(iCurveNum))
			continue;

		GRCurve * curveToInsert = (*itc);

		// create the pointset
		insPointSet * pointSet = new insPointSet();

		GRCurveGeom * curveGeom = curveToInsert->get_curve_geom();

		assert(curveGeom);

		BLS.setCurve(curveToInsert);
		// store the arclength at the new points, and their parameters
		std::vector<double> newPointsArclength;
		std::vector<double> paramToInsert;

		// store the cells, so we don't need to look up
		std::vector<Cell*> containingCells;

		double arcLength = curveGeom->arc_length();
		// parameters
		double minParam = curveGeom->min_param();
		double maxParam = curveGeom->max_param();

		// corresponding length scales
		double currentParam, startLS, nextLS, endLS; // lengthscales
		Cell * startCell, *nextCell, *endCell;
		double length = 0., currentLength = 0.;

		// set the tolerance of distance functions.
		curveGeom->set_tolerance(1.e-6);

		// get start and end lengthscales
		startLS = BLS.getLengthScale(minParam);
		startCell = BLS.getCell();

		endLS = BLS.getLengthScale(maxParam);
		endCell = BLS.getCell();

		// discrete length between coordinates, used as a comparison
		CubitVector cv1, cv2;
		for (GR_index_t iParam = minParam; iParam < maxParam; iParam++) {
			curveGeom->get_coord(iParam, cv1);
			curveGeom->get_coord(iParam + 1, cv2);
			length += (cv1 - cv2).length();
		}

		// This is all initialization for point locations
		// if we want to sample by curvature, we can set lengthscale as the
		// minimum of the lengthscale, and the radius of curvature
		// unless the curvature is zero, in which case, run away.
		//		if (m_sampleByCurvature){
		//			curvature = curveGeom->curvature(minParam);
		//			if (iFuzzyComp(curvature, 0) != 0) startLS = min(startLS, 1./curvature);
		//		}
		paramToInsert.push_back(0.0);
		newPointsArclength.push_back(0.0);

		// if the start length is less than the arc length at the parameter, insert a point
		// or if its a straight line, due to curvature issues
		if (iFuzzyComp(startLS, length) <= 0
				|| curveGeom->get_curve_type() == GRCurveGeom::LINE) {
			currentParam = curveGeom->coord_at_dist(startLS, true,
					tempCubitVector);
		} else {
			currentParam = maxParam;
		}

		// update the params and arclengths, carrying extra data to avoid recomputing
		containingCells.push_back(startCell);

		nextLS = startLS;
		nextCell = startCell;
		currentParam = minParam;
		// while we keep adding initial sampling points

		// if the currentLength at our initial point + the next length to be
		// added is less than the total length - half the endpoint length scale
		// keep going, otherwise stop

		while (currentLength + nextLS < arcLength - 0.5 * endLS) {
			// updates the currParam
			currentParam = curveGeom->param_at_arc_length(currentParam, nextLS);
			// updates cvToInsert
			// recomputes the arc_length
			currentLength = curveGeom->arc_length(0, currentParam);

			paramToInsert.push_back(currentParam);
			newPointsArclength.push_back(currentLength);
			containingCells.push_back(nextCell);
			// update the next lengthscale
			nextLS = BLS.getLengthScale(currentParam);
			nextCell = BLS.getCell();
			//			if (m_sampleByCurvature){
			//				curvature = curveGeom->curvature(currentParam);
			//				if (iFuzzyComp(curvature, 0) != 0)
			//					nextLS = min(startLS, 1./curvature);
			//			}
		}
		paramToInsert.push_back(maxParam);
		newPointsArclength.push_back(arcLength);
		// if there are no points added, we are done
		if (paramToInsert.size() == 2) {
			pointSet->insert(insPoint2D(minParam, startCell));
			//if its not a line, split it once, to ensure we have more than just the end points (helps for splines, circles, etc)
			if (curveGeom->get_curve_type() != GRCurveGeom::LINE) {
				double param = curveGeom->mid_TVT(minParam, maxParam);
				BLS.getLengthScale(param);
				pointSet->insert(insPoint2D(param, BLS.getCell()));
			}
			pointSet->insert(insPoint2D(maxParam, endCell));
			m_pointsToInsert.insert(
					std::pair<GRCurve*, insPointSet*>(curveToInsert, pointSet));

		} else {
			// otherwise, we need to spread points out and optimize point locations

			{
				// number of sampling points for this curve.
				GR_index_t numPoints = paramToInsert.size();
				// lengthscale at the midpoint of the edges
				double * midpointLS = new double[numPoints - 1];
				Cell ** midpointCell = new Cell*[numPoints - 1];
				// ntemp arclength for our jacobi iteration
				double * tempArclength = new double[numPoints];
				// initialize our arclengths
				for (GR_index_t index = 0; index < numPoints; index++) {
					tempArclength[index] = newPointsArclength[index];
					//				printf("%d %f\n",index,newPointsArclength[index]);
				}
				for (GR_index_t index = 0; index < numPoints - 1; index++) {
					midpointCell[index] = containingCells[index];
				}
				// initialize our lengthscales at the midpoints of arcs
				for (GR_index_t index = 0; index < numPoints - 1; index++) {

					midpointLS[index] =
							BLS.getLengthScale(
									curveGeom->param_at_arc_length(
											paramToInsert[index],
											0.5
													* (newPointsArclength[index
															+ 1]
															- newPointsArclength[index])),
									midpointCell[index]);
					midpointCell[index] = BLS.getCell();
					// again, if we sample by curvature, do this

					//				if (m_sampleByCurvature){
					//					curvature = curveGeom->curvature(curveGeom->param_at_coord(tempCubitVector));
					//					if (iFuzzyComp(curvature, 0) != 0)
					//						midpointLS[index] = min(midpointLS[index], 1./curvature);
					//				}
				}

				// lets smooth internal points using equidistribution
				// of length scale
				// uses jacobi on arclength

				double oldArclength = 0.;
				GR_index_t iter = 0;

				// initialize this, to keep track of convergence
				double change_in_arclen = LARGE_DBL;

				while (change_in_arclen > 1.e-2 && ++iter < 1000) {
					change_in_arclen = 0.0;
					for (GR_index_t index = 1; index < numPoints - 1; index++) {
						oldArclength = tempArclength[index]; // used for determining change

						// smoothing
						tempArclength[index] = (midpointLS[index - 1]
								* tempArclength[index + 1]
								+ midpointLS[index] * tempArclength[index - 1])
								/ (midpointLS[index] + midpointLS[index - 1]);
						// compute parameter at that arclength
						paramToInsert[index] = curveGeom->param_at_arc_length(
								paramToInsert[index - 1],
								tempArclength[index]
										- tempArclength[index - 1]);

						// compute the lengthscale at the new location
						midpointLS[index - 1] = BLS.getLengthScale(
								curveGeom->param_at_arc_length(
										paramToInsert[index - 1],
										0.5
												* (tempArclength[index]
														- tempArclength[index
																- 1])),
								midpointCell[index - 1]);
						midpointCell[index - 1] = BLS.getCell();
						//					if (m_sampleByCurvature){
						//						curvature = curveGeom->curvature(curveGeom->param_at_coord(tempCubitVector));
						//						if (iFuzzyComp(curvature, 0) != 0)
						//							midpointLS[index-1] = min(midpointLS[index-1], 1./curvature);;
						//					}
						// update the total change
						change_in_arclen += fabs(
								oldArclength - tempArclength[index])
								/ oldArclength;
					}

				}

				// inserts all the points based on the computed values.
				pointSet->insert(insPoint2D(minParam, startCell));
				for (GR_index_t index = 1; index < numPoints - 1; index++) {
					currentParam = curveGeom->param_at_arc_length(0.0,
							tempArclength[index]);
					BLS.getLengthScale(currentParam, midpointCell[index]);
					pointSet->insert(insPoint2D(currentParam, BLS.getCell()));
				}
				pointSet->insert(insPoint2D(maxParam, endCell));
				m_pointsToInsert.insert(
						std::pair<GRCurve*, insPointSet*>(curveToInsert,
								pointSet));
				// clean up the temporary arrays
				delete[] midpointCell;
				delete[] midpointLS;
				delete[] tempArclength;
			}
		}
	}

	// in debug mode, spit out the set of points
#ifndef NDEBUG
	writeSamplePointsToFile("surfins_samplepoints_initial");
#endif

	if (m_useOffsets)
		sampleOffsetCurves();

	// brief output, return the total number of points to insert
	// well aware if a point shares multiple curves, it is counted twice
	GR_index_t np = 0;

	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); cpm_it++) {
		np += cpm_it->second->size();
	}
	logMessage(1, "%u Point Locations Finalized\n", np);
	return np;
}
// lets get offsets
void SurfaceInserter2D::sampleOffsetCurves() {

	// set up curves, then worry about direction
	std::vector<GRCurve*> firstCurves;
	std::list<GRCurve*> curve_list;
	std::list<GRCurve*>::iterator itc;

	int iCurve = m_offsetCoordinate.first;
	CubitVector oldPoint, newPoint, direction = m_offsetCoordinate.second;
	double offsetDistance, projection, projectionDistance;
	m_geom->get_curves(curve_list);

	if (m_offsetPairs.empty()) {
		int numCurves = m_geom->num_curves();
		for (GR_index_t iOffset = 0; iOffset < m_numberOfOffsets; iOffset++) {
			if (iCurve < 0) {
				std::map<GRCurve *, GRCurve *> offsetcurves =
						m_geom->offset_curves(curve_list, direction.length());
				m_offsetCurveMap.insert(offsetcurves.begin(),
						offsetcurves.end());
				curve_list.clear();
				for (std::map<GRCurve *, GRCurve*>::iterator it =
						offsetcurves.begin(); it != offsetcurves.end(); ++it) {
					curve_list.push_back(it->second);
					firstCurves.push_back(it->first);
				}

			} else {
				std::pair<GRCurve *, GRCurve *> curvePair;
				if (iOffset == 0)
					curvePair = m_geom->offset_curve_linear(iCurve, direction);
				else
					curvePair = m_geom->offset_curve_linear(
							numCurves + iOffset - 1, direction);
				m_offsetCurveMap.insert(curvePair);
				firstCurves.push_back(curvePair.first);
			}
		}
	} else {
		// sets up vector of curves, have to do a little trickery to get things in order
		std::vector<GRCurve*> orderedCurves;

		for (itc = curve_list.begin(); itc != curve_list.end(); itc++) {
			orderedCurves.push_back(*itc);
		}
		for (GR_index_t iP = 0; iP < m_offsetPairs.size(); iP++) {
			GRCurve * oldCurve = orderedCurves[m_offsetPairs[iP].first];
			GRCurve * newCurve = orderedCurves[m_offsetPairs[iP].second];
			m_offsetCurveMap.insert(
					std::pair<GRCurve*, GRCurve*>(oldCurve, newCurve));
			firstCurves.push_back(oldCurve);
		}

	}
//	printf("%d %d\n",firstCurves.size(),m_offsetCurveMap.size());
//	printf("%d %d")

	// the curve vector contains the first curve. This makes sure we sample in order.
	for (GR_index_t iOldCurve = 0; iOldCurve < firstCurves.size();
			iOldCurve++) {
		GRCurve * oldCurve = firstCurves[iOldCurve];
		GRCurve * newCurve = m_offsetCurveMap.find(oldCurve)->second;
		printf("%p %p\n", oldCurve, newCurve);
		// now sample it using the previous curves points, endpoints remain the same
		insertionCurveMap_it cpm_it = m_pointsToInsert.find(oldCurve);
		assert(cpm_it != m_pointsToInsert.end());

		insPointSet * newPointSet = new insPointSet();

		// this becomes trickier since can't project past first one.

		for (insPointSet_it pointSet_it = cpm_it->second->begin();
				pointSet_it != cpm_it->second->end(); ++pointSet_it) {
			// for each point in the set, match up the end points, otherwise, corect based on the normal
			GRCurveGeom * newCurveGeom = newCurve->get_curve_geom();
			double param = pointSet_it->getParam();
			if (pointSet_it == cpm_it->second->begin()) {
				param = newCurveGeom->min_param();
			} else if (pointSet_it == --(cpm_it->second->end())) {
				param = newCurveGeom->max_param();
			} else {
				// lets correct based on the normal.
				if (m_offsetByNormals)
					oldCurve->get_curve_geom()->unit_normal(param, direction);

				oldCurve->get_curve_geom()->coord_at_param(param, oldPoint);

				newCurveGeom->closest_coord_on_curve(oldPoint, newPoint,
						newCurveGeom->min_param(), newCurveGeom->max_param());

				// find the intersect between the direction from the old point to the new curve,
				// using the closest point on the new curve as an estimate of how far away we are
				// and then computing the intersection between a line segment normal distance away
				// and the new curve
				offsetDistance = (newPoint - oldPoint).length();

				projection = direction % (newPoint - oldPoint);
				projectionDistance = projection / fabs(projection) * 2.0
						* offsetDistance;
				// this point is 4 times the distance away from the closest point on the curve
				newPoint = oldPoint + projectionDistance * direction;
				double coords1[3], coords2[3];
				oldPoint.get_xyz(coords1);
				newPoint.get_xyz(coords2);
				std::vector<double> intersectParams;
				bool qIntersect = newCurveGeom->line_intersect(coords1, coords2,
						newCurveGeom->min_param(), newCurveGeom->max_param(),
						intersectParams);
				// if this doesn't produce an intersect, use the closest point
				// otherwise use the first parameter, which is the smallest
				if (qIntersect) {
					// pick the closest point
					newCurveGeom->coord_at_param(intersectParams[0], newPoint);
					param = intersectParams[0];
					double minLength = (newPoint - oldPoint).length();

					for (GR_index_t iParam = 1; iParam > intersectParams.size();
							iParam++) {
						newCurveGeom->coord_at_param(intersectParams[iParam],
								newPoint);
						if ((newPoint - oldPoint).length() < minLength) {
							minLength = (newPoint - oldPoint).length();
							param = intersectParams[iParam];
						}
					}
				} else {
					param = newCurve->get_curve_geom()->closest_param_on_curve(
							oldPoint);
				}

			}
			newPointSet->insert(insPoint2D(param, pointSet_it->getCell()));
		}
		m_pointsToInsert.insert(
				std::pair<GRCurve*, insPointSet*>(newCurve, newPointSet));
	}

#ifndef NDEBUG
	writeSamplePointsToFile("surfins_samplepoints_offset");
#endif

}
int SurfaceInserter2D::removeVertsNearSurface() {

	// make sure we can remove boundary verts
	m_mesh->allowBdryChanges();

	int numVertsToRemove;
	std::vector<Vert*> vertsToRemove;

	// create BLS object -> its small so don't carry one in the class
	BaryLengthScale2D BLS(m_mesh, m_isAnisotropic);

	logMessage(1, "Removing points near the curve\n");

	double radius, LSatPoint1, LSatPoint2;
	double distToEdgeCenter, distToPoint1, distToPoint2;
	insPoint2D insertionPoint1, insertionPoint2;
	CubitVector point1, point2, pointToCheck, center;
	Cell * cell1 = NULL;
	Cell * cell2 = NULL;

	std::set<Cell*> spCInc;
	std::set<Vert*> spVTemp;

	// loop over pairs of potential points on the same curve,
	// these will be internal boundary edges in the final mesh
	// remove all points that would encroach when insertion is performed
	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); cpm_it++) {
		GRCurve * curveToInsert = cpm_it->first;
		BLS.setCurve(curveToInsert);

		// iterators
		insPointSet * pointSet = cpm_it->second;
		insPointSet_it pointSet_end = --pointSet->end();
		insPointSet_it pointSet_it = pointSet->begin();

		// set up pairs of points to iterate over
		insertionPoint1 = (*pointSet_it);
		curveToInsert->position_from_u(insertionPoint1.getParam(), point1);
		LSatPoint1 = BLS.getLengthScale(insertionPoint1.getParam(),
				insertionPoint1.getCell());
		cell1 = BLS.getCell();

		for (; pointSet_it != pointSet_end;) {
			insertionPoint2 = *(++pointSet_it);
			curveToInsert->position_from_u(insertionPoint2.getParam(), point2);
			LSatPoint2 = BLS.getLengthScale(insertionPoint2.getParam(),
					insertionPoint2.getCell());
			cell2 = BLS.getCell();

			center = 0.5 * (point1 + point2);

			if (m_isAnisotropic) {
				// scaled arbitrarily, TODO fix this
				radius = 0.5 * BLS.getAnisoLength(point1, point2, cell1, cell2);
			} else {
				radius = 0.5 * (point1 - point2).length();
			}
			// using the cell the point would be inside
			// look at all three vertices, and collect their
			// neighboring vertices to form a local set of vertices
			// with no vertex more than two vertices away
			std::set<Vert*> neighVerts;

			for (int iV = 0; iV < cell1->getNumVerts(); iV++) {
				findNeighborhoodInfo(cell1->getVert(iV), spCInc, spVTemp);
				neighVerts.insert(spVTemp.begin(), spVTemp.end());
				neighVerts.insert(cell1->getVert(iV));

			}
			for (int iV = 0; cell1 != cell2 && iV < cell2->getNumVerts();
					iV++) {
				findNeighborhoodInfo(cell2->getVert(iV), spCInc, spVTemp);
				neighVerts.insert(spVTemp.begin(), spVTemp.end());
				neighVerts.insert(cell2->getVert(iV));
			}
			// this is checking every vertex

			//			for(GR_index_t iVert = 0; iVert < m_mesh->getNumVerts(); iVert++){
			//							Vert * vertToCheck = m_mesh->getVert(iVert);
			for (std::set<Vert*>::iterator itv = neighVerts.begin();
					itv != neighVerts.end(); ++itv) {
				Vert * vertToCheck = *itv;
				if (vertToCheck->getVertType() != Vert::eBdryApex
						&& !vertToCheck->isDeletionRequested()) {
					// This is where the magic happens
					pointToCheck.set(vertToCheck->getCoords());

					distToPoint1 = (pointToCheck - point1).length();
					distToPoint2 = (pointToCheck - point2).length();
					// treat anisotropic differently
					if (!m_isAnisotropic) {
						distToEdgeCenter = dDIST2D(pointToCheck, center);
					} else {
						distToEdgeCenter = BLS.getAnisoLength(center,
								pointToCheck, cell1, cell2);
						CubitVector vertToPoint1 = ~(pointToCheck - point1);
						CubitVector vertToPoint2 = ~(pointToCheck - point2);
						double scale = vertToCheck->isBdryVert() ? 0.25 : 1.0;
						LSatPoint1 = scale
								* BLS.getAnisoLengthScale(
										insertionPoint1.getParam(),
										insertionPoint1.getCell(),
										&vertToPoint1);
						LSatPoint2 = scale
								* BLS.getAnisoLengthScale(
										insertionPoint2.getParam(),
										insertionPoint2.getCell(),
										&vertToPoint2);
					}
					// check if the point encroaches, if it encroaches then flag it
					// else check the lengthscales
					if (iFuzzyComp(distToEdgeCenter, radius) < 1
							|| LSatPoint1 * 0.5 > distToPoint1
							|| LSatPoint2 * 0.5 > distToPoint2) {
						vertsToRemove.push_back(vertToCheck);
						vertToCheck->markToDelete(); // mark it so we don't need to check again.
					}
				}
			}
			// set the 1st point to the 2nd and continue
			insertionPoint1 = insertionPoint2;
			LSatPoint1 = LSatPoint2;
			point1 = point2;
		}
	}

	// now that we are done flagging verts for removal, set up the metrics for swapping and refinement
	if (m_isAnisotropic) {
		double metric[5];
		setVertOriginalMetrics(); //TODO this function is global, Switch this to local at some point

		// we also need to set the metrics at the sample points, based on the original mesh
		for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
				cpm_it != m_pointsToInsert.end(); cpm_it++) {
			GRCurve * curveToInsert = cpm_it->first;
			for (insPointSet_it pointSet_it = cpm_it->second->begin();
					pointSet_it != cpm_it->second->end(); ++pointSet_it) {
				curveToInsert->position_from_u(pointSet_it->getParam(), point1);
				BLS.getMetric(point1, metric, pointSet_it->getCell());
				pointSet_it->setMetric(metric);
			}
		}
	}

	numVertsToRemove = vertsToRemove.size();

#ifndef NDEBUG
	// outputs the points to be read into paraview
	std::vector<CubitVector> tempPoints;
	for (GR_index_t iV = 0; iV < vertsToRemove.size(); iV++) {
		tempPoints.push_back(CubitVector(vertsToRemove[iV]->getCoords()));
	}
	writeVTKPointSet(tempPoints, NULL, "surfins_removepts", "");
	tempPoints.clear();
#endif

	int numSwaps;
	for (GR_index_t i = 0; i < vertsToRemove.size(); i++) {
		Vert * vert = vertsToRemove[i];
		if (!m_mesh->removeVert(vert, numSwaps)) {
			m_leftoverVertsToRemove.insert(vert);
		}
	}
	m_mesh->sendEvents();

	// Reset queue, clean up mesh
	std::map<Vert*, Vert *> *vertMap = new std::map<Vert*, Vert *>();
	m_mesh->purgeAllEntities(vertMap);
	SmoothMan2D.emptyQueue();
	updateVertsToRemoveAfterInsertion(*vertMap);
	delete vertMap;

	// if its anisotropic, we refine differently and need to
	// clear the queue to prevent buildup
	if (!m_isAnisotropic)
		RefMan2D->emptyQueue();

	m_mesh->writeTempMesh();

	//update the cells stored with the sample points to contain the current cells
	updateSamplePoints();

	// gets internal boundary intersections if there are any
	if (m_mesh->getNumIntBdryEdgesInUse() > 0) {
		sampleForIntBdryIntersections();
	}

	// in debug mode, output lots the final point locations
#ifndef NDEBUG
	writeSamplePointsToFile("surfins_samplepoints");
#endif

	// Update the statistics
	numRemovedVerts = numInitVerts - m_mesh->getNumVerts();
	return numVertsToRemove;
}
void SurfaceInserter2D::updateSamplePoints() {

	logMessage(1, "Updating Sample Points\n");
	CubitVector point;
	Cell * containingCell = m_mesh->getCell(0);
	BaryLengthScale2D BLS(m_mesh);

	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); cpm_it++) {
		GRCurve * curveToInsert = cpm_it->first;
		insPointSet * pointSet = cpm_it->second;
		for (insPointSet_it pointSet_it = pointSet->begin();
				pointSet_it != pointSet->end(); ++pointSet_it) {

			curveToInsert->position_from_u(pointSet_it->getParam(), point);
//			containingCell = pointSet_it->getCell();
//			if(!containingCell->isValid() || !containingCell->doFullCheck())
			containingCell = BLS.findCell(point, NULL);
//			else
//				containingCell = BLS.findCell(point,containingCell);

			// this should be really fast, cause we are moving along the curve.
			// still might benefit from a search tree here,
			assert(!containingCell->isDeleted());

			pointSet_it->setCell(containingCell);

		}
	}
}

void SurfaceInserter2D::sampleForIntBdryIntersections() {

	CubitVector point, point1, point2;
	int numViolated = 0, numTies = 0;
	Face *facesViolated[3], *facesTied[3];

	logMessage(1, "Looking for intersection points with internal boundaries\n");
	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); ++cpm_it) {

		std::vector<insPoint2D> pointsToDelete;
		std::vector<insPoint2D> pointsToInsert;

		GRCurve * curveToInsert = cpm_it->first;
		GRCurveGeom * curveGeom = curveToInsert->get_curve_geom();
		insPointSet * pointSet = cpm_it->second;

		insPointSet_it pointSet_last = --pointSet->end();

		insPoint2D insPoint1 = *(pointSet->begin());
		insPoint2D insPoint2;

		for (insPointSet_it pointSet_it = pointSet->begin();
				pointSet_it != pointSet_last;) {
			insPoint2 = *(++pointSet_it);
			// if consecutive points have different regions, then we have an intersection
			// lets find this easily
			if (insPoint1.getCell()->getRegion()
					!= insPoint2.getCell()->getRegion()) {
				// need to find the intersections.
				// use the mesh's version, which should hit closest boundary face;
				double coords1[2], coords2[2];
				curveGeom->coord_at_param(insPoint1.getParam(), point1);
				point1.get_xyz(coords1);
				curveGeom->coord_at_param(insPoint2.getParam(), point2);
				point2.get_xyz(coords2);
				Cell * cell = m_mesh->findCell(coords2, insPoint1.getCell(),
						numViolated, facesViolated, numTies, facesTied);
				if (numViolated == 0)
					assert(0);
				cell = facesViolated[0]->getOppositeCell(cell);
				assert(facesViolated[0]->isBdryFace());
				IntBdryEdge * IBE = dynamic_cast<IntBdryEdge*>(cell);
				std::vector<double> intersectParams;
				bool qIntersect = curveGeom->line_intersect(
						IBE->getFace(0)->getVert(0)->getCoords(),
						IBE->getFace(0)->getVert(1)->getCoords(),
						curveGeom->min_param(), curveGeom->max_param(),
						intersectParams);

				if (qIntersect == false) {
					// we might have found the wrong boundary face
					// look at its neighbors, walk along them to find intersections
					// used for findNeighborhoodInfo
					std::set<Cell*> spCIncident;
					std::set<Vert*> spVCand;
					std::set<BFace*> spBFIncident;

					Vert *pV0 = IBE->getVert(0), *pV1 = IBE->getVert(1);
					IntBdryEdge * IBE0 = IBE, *IBE1 = IBE;
					bool qBdryVert = true;
					// while we are in our chain, should be quick to find the nearest interesecting
					// boundary face, only looking at sides with valid boundary faces.

					while (pV0->isValid() || pV1->isValid()) {
						// each vert corresponds to side of the chain
						if (pV0) {
							findNeighborhoodInfo(pV0, spCIncident, spVCand,
									&spBFIncident, &qBdryVert);
							if (spBFIncident.size() != 2) {
								pV0 = NULL;
								continue;
							}
							// continue our chain along the BFace we aren't already on
							IntBdryEdge * incBFace =
									dynamic_cast<IntBdryEdge*>(*spBFIncident.begin());
							if (incBFace == IBE0 && spBFIncident.size() == 2)
								incBFace =
										dynamic_cast<IntBdryEdge*>(*(++spBFIncident.begin()));
							IBE0 = incBFace;
							pV0 = IBE0->getFace(0)->getOppositeVert(pV0);
							qIntersect = curveGeom->line_intersect(
									IBE0->getFace(0)->getVert(0)->getCoords(),
									IBE0->getFace(0)->getVert(1)->getCoords(),
									curveGeom->min_param(),
									curveGeom->max_param(), intersectParams);
							if (qIntersect) {
								IBE = IBE0;
								break;
							}
						}
						if (pV1) {
							findNeighborhoodInfo(pV1, spCIncident, spVCand,
									&spBFIncident, &qBdryVert);
							if (spBFIncident.size() != 2) {
								pV1 = NULL;
								continue;
							}
							IntBdryEdge * incBFace =
									dynamic_cast<IntBdryEdge*>(*spBFIncident.begin());
							if (incBFace == IBE1 && spBFIncident.size() == 2)
								incBFace =
										dynamic_cast<IntBdryEdge*>(*(++spBFIncident.begin()));
							IBE1 = incBFace;
							pV1 = IBE1->getFace(0)->getOppositeVert(pV1);
							qIntersect = curveGeom->line_intersect(
									IBE1->getFace(0)->getVert(0)->getCoords(),
									IBE1->getFace(0)->getVert(1)->getCoords(),
									curveGeom->min_param(),
									curveGeom->max_param(), intersectParams);
							if (qIntersect) {
								IBE = IBE1;
								break;
							}
						}
					}
				}

				for (GR_index_t i = 0; i < intersectParams.size(); i++) {
					curveGeom->coord_at_param(intersectParams[i], point);
					if ((point2 - point).length() > (point1 - point).length()) {
						pointsToDelete.push_back(insPoint1);
					} else {
						pointsToDelete.push_back(insPoint2);
					}
					pointsToInsert.push_back(
							insPoint2D(intersectParams[i], IBE));

				}
			}

			insPoint1 = insPoint2;
		} // pointSet
		  // wait till we are done.
		for (GR_index_t iP = 0; iP < pointsToDelete.size(); iP++) {
			pointSet->erase(pointsToDelete[iP]);
		}
		for (GR_index_t iP = 0; iP < pointsToInsert.size(); iP++) {
			pointSet->insert(pointsToInsert[iP]);
		}
	}
}

int SurfaceInserter2D::insertBoundaryPointsIntoMesh() {

	logMessage(1, "Inserting sample point locations\n");

	// use watson inserter for boundaries
	WatsonInserter WI(m_mesh);

	// need to know length
	BaryLengthScale2D BLS(m_mesh, m_isAnisotropic);
	// The previous point is checked to see if its been inserted in the previous location
	// as for a set of single connected curves
	double newPtLoc[2];
	Vert * newVert = pVInvalidVert;

	// used in a lot of places
	CubitVector point;

	// iterator
	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); cpm_it++) {
		GRCurve * curveToInsert = cpm_it->first;
		insPointSet * pointSet = cpm_it->second;
		for (insPointSet_it pointSet_it = pointSet->begin();
				pointSet_it != pointSet->end(); ++pointSet_it) {

			curveToInsert->position_from_u(pointSet_it->getParam(), point);
			point.get_xyz(newPtLoc);

			double splitParam = pointSet_it->getParam();
			// this should be really fast, cause we are moving along the curve.
			// still might benefit from a search tree here,

			BdryEdgeBase * bdryCell = NULL;

			for (GR_index_t i = 0; i < m_mesh->getNumTotalBdryFaces(); ++i) {
				BdryEdgeBase* bdryEdge =
						dynamic_cast<BdryEdgeBase*>(m_mesh->getBFace(i));
				if (!bdryEdge->doFullCheck())
					continue;
				double startParam = bdryEdge->getVert0Param();
				double endParam = bdryEdge->getVert1Param();
				if (bdryEdge->getCurve() == curveToInsert
						&& ((startParam < splitParam && endParam > splitParam)
								|| (startParam > splitParam
										&& endParam < splitParam))) {
					bdryCell = bdryEdge;
					break;
				}
			}
			if (bdryCell == NULL)
				continue;
			assert(bdryCell->doFullCheck());
			if (iFuzzyComp(bdryCell->getVert0Param(), splitParam) == 0) {
				newVert = bdryCell->getVert(0);
				m_leftoverVertsToRemove.erase(newVert);
				continue;
			}
			if (iFuzzyComp(bdryCell->getVert1Param(), splitParam) == 0) {
				newVert = bdryCell->getVert(1);
				m_leftoverVertsToRemove.erase(newVert);
				continue;
			}
			// we need to split at a boundary cell

			WI.computeHull(newPtLoc,
					bdryCell->getFace(0)->getOppositeCell(bdryCell),
					dynamic_cast<BFace*>(bdryCell), splitParam);
			newVert = WI.insertPointInHull();
			newVert->setType(Vert::eBdry);
//			newVert->setType(Vert::eBdryApex);
//			m_mesh->writeTempMesh();

			assert(newVert->isValid());
			newVert->setParentEntity(curveToInsert);
			pointSet_it->setVert(newVert); // change this to point to the vertex
		}
	}
	m_mesh->sendEvents();

#ifndef NDEBUG
//	writeVTKLegacy(*m_mesh,"surfins_recover");
//	m_mesh->writeTempMesh();
#endif

	numInsertedVerts = m_mesh->getNumVerts() - numInitVerts + numRemovedVerts;
	return numInsertedVerts;

}
int SurfaceInserter2D::insertPointsIntoMesh() {

	logMessage(1, "Inserting sample point locations\n");

	// use watson inserter for boundaries
	WatsonInserter WI(m_mesh);

	// use swapping inserter elsewhere, cheaper for now
//	LengthAniso* pLength = new LengthAniso2D(LengthAniso2D::eMesh,1,1,LARGE_DBL);
//	pLength->provideMesh(m_mesh);
	LengthAniso2D* pLA2D = dynamic_cast<LengthAniso2D*>(m_anisoLength);
	assert(pLA2D);
	SwappingInserter2D SwapIns2D(m_mesh, SwapMan2D, pLA2D);

	// need to know length
	BaryLengthScale2D BLS(m_mesh, m_isAnisotropic);
	// The previous point is checked to see if its been inserted in the previous location
	// as for a set of single connected curves
	double newPtLoc[] = { 0, 0, 0 };
	Vert * newVert = pVInvalidVert;

	std::vector<std::pair<Vert*, Vert*> > edgesToRecover;

	// used in a lot of places
	Face * face = pFInvalidFace;
	Cell * containingCell = m_mesh->getCell(0);
	CubitVector point;

	// Set Up Insertion Objects
	SwapIns2D.setForcedInsertion(true);
	SwapMan2D->clearQueue();
	// iterator
	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); cpm_it++) {
		GRCurve * curveToInsert = cpm_it->first;
		insPointSet * pointSet = cpm_it->second;
		insPointSet_it pointSet_last = --pointSet->end();
		for (insPointSet_it pointSet_it = pointSet->begin();
				pointSet_it != pointSet->end(); ++pointSet_it) {

			bool qEndPoint = (pointSet_it == pointSet->begin()
					|| pointSet_it == pointSet_last);
			curveToInsert->position_from_u(pointSet_it->getParam(), point);
			point.get_xyz(newPtLoc);

			// if the stored containingCell is valid, use it as the first guess
			// might not be cause of point insertion
			if (pointSet_it->getCell()->isValid()
					&& !pointSet_it->getCell()->isDeleted()
					&& pointSet_it->getCell()->getType() == Cell::eTriCell)
				containingCell = pointSet_it->getCell();
			// perform the insert
			assert(containingCell->doFullCheck());

			// this should be really fast, cause we are moving along the curve.
			// still might benefit from a search tree here,

			containingCell = BLS.findCell(point, containingCell);
			Cell * bdryCell = NULL;

			if (pointSet_it->getCell()->getType() == Cell::eIntBdryEdge)
				bdryCell = pointSet_it->getCell();
			// gets the nearest boundary containingCell, if there is one

			// for endpoints, we want to check if they intersect with a boundary face
			// as we can't check boundary faces prior to removal (they may change)
			// we can quickly check from the insertion cell
			// this is quick, and only done at endpoints of curves
			if (qEndPoint) {
				for (int iF = 0; iF < 3; iF++) {
					if (containingCell->getFace(iF)->getOppositeCell(
							containingCell)->isBdryCell()) {
						// project the point onto the boundary face
						double projection[] = { newPtLoc[0], newPtLoc[1] };
						containingCell->getFace(iF)->projectOntoFace(
								projection);
						// if its at a vertex continue, this will be quickly dealt with upon insertion;
						if (dDIST2D(projection,
								containingCell->getFace(iF)->getVert(0)->getCoords())
								< 1.e-13
								|| dDIST2D(projection,
										containingCell->getFace(iF)->getVert(1)->getCoords())
										< 1.e-13)
							break;

						// if we aren't at a vertex, compare the distance between the projection and the face
						if (dDIST2D(projection,newPtLoc) < 1.e-6) {
							bdryCell =
									containingCell->getFace(iF)->getOppositeCell(
											containingCell);
							break;
						}
					}
				}
			}
			if (bdryCell) {
				// we need to split at a boundary cell
				double splitParam =
						dynamic_cast<BdryEdgeBase*>(bdryCell)->getCurve()->u_from_position(
								CubitVector(newPtLoc));
				WI.computeHull(newPtLoc,
						bdryCell->getFace(0)->getOppositeCell(bdryCell),
						dynamic_cast<BFace*>(bdryCell), splitParam);
				newVert = WI.insertPointInHull();
				newVert->setType(Vert::eBdryApex);

			} else {
				// its inside a cell
				assert(
						isInTri(newPtLoc,
								containingCell->getVert(0)->getCoords(),
								containingCell->getVert(1)->getCoords(),
								containingCell->getVert(2)->getCoords()) > 0);

				newVert = SwapIns2D.insertPoint(newPtLoc, containingCell);

				if (qEndPoint)
					newVert->setType(Vert::eBdryApex);
				else
					newVert->setType(Vert::eBdryTwoSide);

			}

			if (!newVert) {
				// if for whatever reason, swappping inserter didn't work, then use watson inserter
				// not sure if this is ever called, but its a safety precaution
				WI.computeHull(newPtLoc, containingCell);
				newVert = WI.insertPointInHull();
				newVert->setType(Vert::eBdryApex);
			}
			m_mesh->writeTempMesh();

			assert(newVert->isValid());

			newVert->setParentEntity(curveToInsert);
			pointSet_it->setVert(newVert); // change this to point to the vertex

			if (m_isAnisotropic) {
				// Setting metric value only matters if you're going to apply them to the vertex.
				// That's done in the setLengthScale call, apparently.
//				double metric[5];
//				for(int iM = 0; iM < 3; iM++)
//					metric[iM] = pointSet_it->getMetric(iM);
				// never use the last two entries, so ignore
				//newVert->setMetric(metric);
				m_anisoLength->setLengthScale(newVert);
			}
		}

		// iterate across pairs
		insPoint2D pt1 = (*pointSet->begin()), pt2;
		for (insPointSet_it pointSet_it = pointSet->begin();
				pointSet_it != pointSet_last;) {
			pt2 = *(++pointSet_it);

			if (pt1.getVert() == pt2.getVert())
				continue;

			// if we can't recover for some reason, then try again later. Almost always works.
			if (m_mesh->recoverEdge(pt1.getVert(), pt2.getVert(), face))
				face->lock();
			else
				edgesToRecover.push_back(
						std::pair<Vert*, Vert*>(pt1.getVert(), pt2.getVert()));

#ifndef NDEBUG
			m_mesh->writeTempMesh();
#endif

			pt1 = pt2;
		}
		m_mesh->sendEvents();
	}
	// now remove all verts we couldn't, and attempt to recover edges we could not.
	// has never not worked.
	removeVertsAfterInsertion();
	for (GR_index_t iEdge = 0; iEdge < edgesToRecover.size(); iEdge++) {
		if (m_mesh->recoverEdge(edgesToRecover[iEdge].first,
				edgesToRecover[iEdge].second, face))
			face->lock();
		else {
			assert(0);
		}
	}

	writeVTKLegacyWithoutPurge(*m_mesh, "surfins_insert");

	m_mesh->sendEvents();
	// this handles offsets and quad cells.
	if (m_useOffsets) {
		createOffsetQuadCells();
	}
	// create internal boundary edges, formally.
	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); ++cpm_it) {
		GRCurve * curveToInsert = cpm_it->first;
		insPointSet * pointSet = cpm_it->second;
		insPointSet_it pointSet_last = --pointSet->end();
		insPoint2D insPoint1 = (*pointSet->begin());
		insPoint2D insPoint2;
		for (insPointSet_it pointSet_it = pointSet->begin();
				pointSet_it != pointSet_last;) {
			insPoint2 = *(++pointSet_it);
			if (insPoint1.getVert() == insPoint2.getVert())
				continue;
			assert(insPoint1.getVert()->isValid());
			assert(insPoint2.getVert()->isValid());

			IntBdryEdge * IBE = m_mesh->createIntBdryEdge(insPoint1.getVert(),
					insPoint2.getVert(), curveToInsert, insPoint1.getParam(),
					insPoint2.getParam());
			assert(IBE->doFullCheck());
			insPoint1 = insPoint2;
		}
	}
	m_mesh->sendEvents();

#ifndef NDEBUG
	writeVTK(*m_mesh, "surfins_recover");
	m_mesh->writeTempMesh();
#endif

	numInsertedVerts = m_mesh->getNumVerts() - numInitVerts + numRemovedVerts;
	return numInsertedVerts;

}
void SurfaceInserter2D::refineAndSmooth(const char strParams[]) {

	// lets build queues if we are anisotropic, since aniso stuff doesn't have
	// the ability to use a local queue.

	std::set<Face*> swapQueue;

	if (m_isAnisotropic) {

		for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
				cpm_it != m_pointsToInsert.end(); cpm_it++) {
			//		GRCurve * curveToInsert = cpm_it->first;
			insPointSet * pointSet = cpm_it->second;

			for (insPointSet_it pointSet_it = pointSet->begin();
					pointSet_it != pointSet->end(); ++pointSet_it) {
				Vert * vert = pointSet_it->getVert();
				if (!vert)
					continue;
				std::set<Cell*> spCInc;
				std::set<Vert*> spVTemp;
				findNeighborhoodInfo(vert, spCInc, spVTemp);
				for (std::set<Cell*>::iterator itc = spCInc.begin();
						itc != spCInc.end(); ++itc) {
					for (int iV = 0; iV < (*itc)->getNumFaces(); iV++) {
						Face * face = (*itc)->getFace(iV);
						if (face->isBdryFace())
							continue;

						swapQueue.insert(face);
					}
				}
			}
		}
	}

	// rest is same as normal GRUMMP
	GR_index_t iLen = strlen(strParams), iParam = 0;
	//	setVertMetrics(*m_mesh);
	numRefinedVerts = m_mesh->getNumVerts();

	// clear events queue
	m_mesh->sendEvents();
	while (iParam < iLen) {
		char cType = strParams[iParam++];
		assert(cType);
		switch (cType) {
		case 'W':
		case 'w':  // Do some swapping
		{
			char cMeasure = strParams[iParam++];
			GRUMMP::SwapDecider2D *pSD2D = NULL;
			switch (cMeasure) {
			case 'I':
			case 'i': // Insphere criterion
				logMessage(1, "Doing Delaunay face swapping...\n");
				pSD2D = new GRUMMP::DelaunaySwapDecider2D();
				break;
			case 'M':
			case 'm':
			case 'S':
			case 's': // MinmaxAngle
				if (!m_isAnisotropic) {
					logMessage(1,
							"Doing face swapping to decrease max angle...\n");
					pSD2D = new GRUMMP::MinMaxAngleSwapDecider2D();
				} else {
					logMessage(1, "Doing aniso face swapping...\n");
					pSD2D = SwapDec2D;
				}
				break;
			default:
				vFatalError("Bad option for swapping criterion", "main()");
				break;
			}
			if (m_isAnisotropic) {
				int nSwaps = 0;
				for (std::set<Face*>::iterator itf = swapQueue.begin();
						itf != swapQueue.end(); ++itf) {
					Face * face = *itf;
					nSwaps += ARMesh2D->iAniSwap(face);
				}
				logMessage(1, "Swapped %d times anisotropically\n", nSwaps);
			} else {
				SwapMan2D->changeSwapDecider(pSD2D);
				int nSwaps = SwapMan2D->swapAllFaces(false);
				logMessage(1, "Swapped %d faces\n", nSwaps);
				m_mesh->allowBdryChanges();
				SwapMan2D->clearQueue();
			}
			if (pSD2D != SwapDec2D)
				delete pSD2D;

		}
			break;
		case 'M':
		case 'm': // Do some smoothing
		{
			iParam = iParam + 2;
			if (m_isAnisotropic) {
				//				GRUMMP::SmoothingManagerMesquite2D *pSMM2D = new GRUMMP::SmoothingManagerMesquite2D(m_mesh);
				//				m_mesh->disallowBdryChanges();
				//				//				SmoothMan2D.setSmoothingGoal(5);
				//				//				int nsmooth = 0;
				//				//				for(std::set<Vert*>::iterator itv = smoothQueue.begin(); itv != smoothQueue.end(); ++itv){
				//				//					nsmooth += pSMM2D->smoothVert(*itv);
				//				//				}
				////				pSMM2D->vSetTargetfromMetric();
				////				int nsmooth = pSMM2D->smoothAllVerts(1);
				////				printf("aniso smoothed %d verts\n",nsmooth);
				//				delete pSMM2D;
				//				m_mesh->sendEvents();
				//				m_mesh->allowBdryChanges();

			} else {
				//			m_mesh->disallowBdryChanges();
				//			SmoothMan2D.setSmoothingGoal(5);
				SmoothMan2D.setSmoothingThreshold(30.);
				SmoothMan2D.smoothAllQueuedVerts();
				m_mesh->sendEvents();
				//			m_mesh->allowBdryChanges();
			}
		}
			break;
		case 'r':
		case 'R': {
			if (m_isAnisotropic) {
				// currently don't do this
				/*
				 double dRho=1.5;
				 double dThetaMin=28.0*M_PI/180.0;
				 double dAdes=ARMesh2D->dAvgArea()/dRho;
				 double dIdeal=2*sqrt(dAdes/sqrt(3));
				 double dTau=dIdeal*(7.0-4.0*sqrt(3.0)*sin(dThetaMin)-4.0*cos(dThetaMin)*cos(dThetaMin))/(4.0*cos(dThetaMin)*cos(dThetaMin)-1.0);
				 double dMaxCircRad=1.0/sqrt(3)*(dIdeal+(dTau));
				 double MinLength=dIdeal-dTau;

				 assert(dMaxCircRad>0);
				 assert(MinLength>0);
				 CellQueue* pCellQueue = new CellQueue(dynamic_cast<Mesh*>(m_mesh),dMaxCircRad,MinLength,true);
				 ARMesh2D->vImproveMesh(pCellQueue,200000000,0.4,false); //0.2 works well
				 delete pCellQueue;
				 */
			} else {

				// protect first
				if (!m_useOffsets) {
					logMessage(1, "Protecting small angle vertices\n");
					for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
							cpm_it != m_pointsToInsert.end(); cpm_it++) {
						insPointSet * pointSet = cpm_it->second;
						for (insPointSet_it pointSet_it = pointSet->begin();
								pointSet_it != pointSet->end(); ++pointSet_it) {
							if (!pointSet_it->getVert()->isSmallAngleVert()
									&& pointSet_it->getVert()->getVertType()
											== Vert::eBdryApex) {
								identifyAndFixSmallAngle(
										pointSet_it->getVert());
							}
						}
					}
				}
				m_mesh->purgeAllEntities();

				m_mesh->evaluateQuality();

				writeVTKLegacyWithoutPurge(*m_mesh, "surfins_prerefine");

				RefMan2D->fixBoundaryEncroachment();
				RefMan2D->refineAllQueuedCells();

				writeVTKLegacyWithoutPurge(*m_mesh, "surfins_postrefine");

				m_mesh->purgeAllEntities();
				m_mesh->evaluateQuality(); //post refine
			}
		}
			break;
		}

	}

	if (SwapMan2D)
		SwapMan2D->changeSwapDecider(SwapDec2D);
	numRefinedVerts = m_mesh->getNumVerts() - numRefinedVerts;

}

void SurfaceInserter2D::setRegions() {
	IntBdryEdge * IBE = NULL;
	Cell * cell0 = NULL;
	Cell * cell1 = NULL;
	// first find max region, then use one higher
	int iRegion, iLeftReg, iRightReg;
	iRegion = iLeftReg = iRightReg = 0;
	for (GR_index_t iC = 0; iC < m_mesh->getNumCells(); iC++) {
		iRegion = max(iRegion, m_mesh->getCell(iC)->getRegion());
	}
	iRegion++;

	// not interested in optimal coloring, so potentially could give large numbers of regions
	// and take a long time
	for (GR_index_t iB = 0; iB < m_mesh->getNumIntBdryFaces(); iB++) {

		IBE = m_mesh->getIntBFace(iB);

		cell0 = IBE->getFace(0)->getOppositeCell(IBE);
		cell1 = IBE->getFace(1)->getOppositeCell(IBE);
		iLeftReg = IBE->getCurve()->region_left();
		iRightReg = IBE->getCurve()->region_right();
		// this case we have prescribed regions, which we will attempt to obey
		if (iLeftReg != iRightReg) {
			// make sure the regions are correct, with a quick check
			if (cell0 == IBE->getFace(0)->getLeftCell()) {
				if (cell0->getRegion() != iLeftReg)
					setRegionFromCell(cell0, iLeftReg, m_mesh);
				if (cell1->getRegion() != iRightReg)
					setRegionFromCell(cell1, iRightReg, m_mesh);
			} else {
				if (cell1->getRegion() != iLeftReg)
					setRegionFromCell(cell1, iLeftReg, m_mesh);
				if (cell0->getRegion() != iRightReg)
					setRegionFromCell(cell0, iRightReg, m_mesh);
			}
		} else if (cell0->getRegion() == cell1->getRegion()) {

			// if they have the same regions
			// traverse the mesh and reset
			// otherwise, assume things are good

			if (cell0->getRegion() != iRegion)
				setRegionFromCell(cell0, iRegion, m_mesh);
			else if (cell1->getRegion() != iRegion)
				setRegionFromCell(cell1, iRegion, m_mesh);

			logMessage(1, "Relabelling new region as region %d\n", iRegion);
			if (iRegion++ == iMaxRegion - 1) {
				logMessage(1, "Curves are not closed\n");
				break;
			}
		}
	}

}

// this is basically identical to processing input angles.
void SurfaceInserter2D::identifyAndFixSmallAngle(Vert * apex_vert) {

	//Every apex vertex potentially harbors one or many small angle complexes.
	//Find them and protect them by a shell.

	std::set<BFace*> small_angle_complex;

	//Find the boundary edges attached to apex_vert. Construct a map which
	//orders these edges as function of the plane angle formed by the tangent.

	std::set<std::pair<double, BFace*> > angle_edge_pairs;

	CubitVector apex_coords(apex_vert->x(), apex_vert->y(), 0.);

	//Since edge to vert connectivity is not stored, we have to go through all the edges
	//to find those bounded by apex_vert. Fortunately, at this point, we have a small amount
	//of edges so this should not be too expensive.

	double min_length = LARGE_DBL, check_length;
	std::set<Cell*> spCInc;
	std::set<Vert*> spVTemp;
	std::set<BFace*> pspBFInc;
	bool qBdryVert = true;
	findNeighborhoodInfo(apex_vert, spCInc, spVTemp, &pspBFInc, &qBdryVert);
	std::set<BFace*>::iterator ite, ite_end = pspBFInc.end();

	if (pspBFInc.size() < 2)
		return;
	for (ite = pspBFInc.begin(); ite != ite_end; ++ite) {
		BFace * bface = *ite;
		if (bface->isDeleted())
			continue;

		check_length = bface->calcSize();
		min_length = check_length < min_length ? check_length : min_length;
		small_angle_complex.insert(bface);

	}
	//	double check_dist = 0.01 * min_length;
	ite_end = small_angle_complex.end();

	for (ite = small_angle_complex.begin(); ite != ite_end; ++ite) {

		BFace * bface = *ite;
		if (bface->isDeleted())
			continue;
		Vert * oppVert =
				bface->getVert(0) == apex_vert ?
						bface->getVert(1) : bface->getVert(0);

		CubitVector check_coords(oppVert->getCoords());

		double angle = atan2(check_coords.y() - apex_coords.y(),
				check_coords.x() - apex_coords.x());

		angle_edge_pairs.insert(std::make_pair(angle, bface));

	}

	if (angle_edge_pairs.size() < 2)
		return;

	//Add the first entry to the end of the map (with an offset of -2*M_PI)
	//to be able to complete a cycle around the apex_vertex when looking for
	//small angles.

	angle_edge_pairs.insert(
			std::make_pair(angle_edge_pairs.begin()->first + 2. * M_PI,
					angle_edge_pairs.begin()->second));

	//Go through the map and identify pairs of boundary edges forming small angles.
	std::vector<std::pair<BFace*, BFace*>> edges_around_small_angle;

	std::set<std::pair<double, BFace*> >::iterator ita, ita_end =
			--angle_edge_pairs.end();

	for (ita = angle_edge_pairs.begin(); ita != ita_end;) {
		std::pair<double, BFace*> pair1 = *ita, pair2 = *(++ita);
		double angle1 = pair1.first, angle2 = pair2.first;

		//If the boundary edges are separated by an angle of more that 60 degrees,
		//then they cannot form a small angle. Otherwise, we need to protect the small angle.

		if ((angle2 - angle1) > M_PI / 3.)
			continue;

		//If the boundary edges do not share a common region where the small angle
		//is located, then they do not form a small angle.

		BFace *bface1 = pair1.second, *bface2 = pair2.second;

		//Found a pair of boundary edges forming a small angle. Insert them in
		//a small angle map (make sure no duplicate entries exist in the map).
		//Second entry is always on CW side of first entry.
		edges_around_small_angle.push_back(std::make_pair(bface1, bface2));

	}
	if (edges_around_small_angle.empty())
		return;

	//edges_of_small_angle contains pairs of boundary edges sharing a small angle.
	//cross product of tangent(first) and tangent(second) is always positive.
	//Now time to form small angle complexes and split edges on concentric shells.
	while (!edges_around_small_angle.empty()) {

		small_angle_complex.clear();

		BFace* bface0 = pBFInvalidBFace, *bface1 = pBFInvalidBFace;

		do {

			bface0 = edges_around_small_angle.back().first;
			bface1 = edges_around_small_angle.back().second;

			edges_around_small_angle.pop_back();

			small_angle_complex.insert(bface0);
			small_angle_complex.insert(bface1);

		} while (bface1 == edges_around_small_angle.front().first
				&& !edges_around_small_angle.empty());

		if (small_angle_complex.empty())
			return;

		apex_vert->markAsSmallAngleVert(true);
		logMessage(2, "small angle vert %f %f\n", apex_vert->x(),
				apex_vert->y());

		// We'll use a watson inserter to split the boundary
		WatsonInserter WI(m_mesh);

		std::set<BFace*> all_edges;

		double shell_radius = LARGE_DBL;

		BFace* bface = pBFInvalidBFace;

		std::set<BFace*>::const_iterator it = small_angle_complex.begin(),
				it_end = small_angle_complex.end();

		for (; it != it_end; ++it) {

			bface = *it;
			if (bface->isDeleted())
				continue;
			if (bface->getVert(0)->getVertType() == Vert::eBdryApex
					&& bface->getVert(1)->getVertType() == Vert::eBdryApex) {
				BdryEdgeBase * bdryFace = dynamic_cast<BdryEdgeBase*>(bface);
				double newLoc[3];
				bool concentricSplit;
				CubitVector split_location;

				double splitParam = bdryFace->getSplitData(
						BdryEdgeBase::MID_PARAM, split_location,
						concentricSplit);
				newLoc[0] = split_location.x();
				newLoc[1] = split_location.y();
				newLoc[2] = 0.0;
				bool qIntBdry = bface->getType() == Cell::eIntBdryEdge;
				WI.computeHull(newLoc,
						bface->getFace(0)->getOppositeCell(bface), bface,
						splitParam);

				Vert * newVert = WI.insertPointInHull();

				if (qIntBdry)
					newVert->setType(Vert::eBdryTwoSide);
				else
					newVert->setType(Vert::eBdry);

				newVert->markAsShellVert(true);
				m_isoLength->setLengthScale(newVert);
				m_mesh->writeTempMesh();

				// figure out which new edge to add
				findNeighborhoodInfo(newVert, spCInc, spVTemp, &pspBFInc,
						&qBdryVert);
				if ((*pspBFInc.begin())->getVert(0) == apex_vert
						|| (*pspBFInc.begin())->getVert(1) == apex_vert) {
					bface = (*pspBFInc.begin());
				} else
					bface = *(++pspBFInc.begin());

				all_edges.insert(bface);

			} else {
				// bface already been split, so just pass it as is
				all_edges.insert(bface);

			}

			double bfacesize = bface->calcSize();
			shell_radius = std::min(bfacesize, shell_radius);

		}

		std::set<BFace*>::iterator itv, itv_end = all_edges.end();

		double distance = 0.5 * shell_radius;
		for (itv = all_edges.begin(); itv != itv_end; ++itv) {

			bface = *itv;
			// don't split if its too long

			assert(
					bface->getVert(0)->isSmallAngleVert()
							|| bface->getVert(1)->isSmallAngleVert());

			BdryEdgeBase * bdryFace = dynamic_cast<BdryEdgeBase*>(bface);

			GRCurveGeom * cg = bdryFace->getCurve()->get_curve_geom();
			Vert * oppVert =
					bface->getVert(0) == apex_vert ?
							bface->getVert(1) : bface->getVert(0);

			double apexParam = cg->param_at_coord(
					CubitVector(apex_vert->getCoords()));
			double oppVertParam = cg->param_at_coord(
					CubitVector(oppVert->getCoords()));

			double newLoc[3];
			CubitVector split_location;
			double splitParam;
			assert(bface->calcSize() > distance);
			// check which direction we need to figure things out in.
			if (iFuzzyComp(apexParam, oppVertParam) < 1) {
				splitParam = cg->param_at_arc_length(apexParam, distance);
				cg->coord_at_param(splitParam, split_location);
			} else {
				splitParam = cg->param_at_arc_length(oppVertParam,
						bface->calcSize() - distance);
				cg->coord_at_param(splitParam, split_location);
			}

			newLoc[0] = split_location.x();
			newLoc[1] = split_location.y();
			newLoc[2] = split_location.z();

			bool qIntBdry = bface->getType() == Cell::eIntBdryEdge;
			WI.computeHull(newLoc, bface->getFace(0)->getOppositeCell(bface),
					bface, splitParam);

			Vert * newVert = WI.insertPointInHull();

			if (qIntBdry)
				newVert->setType(Vert::eBdryTwoSide);
			else
				newVert->setType(Vert::eBdry);

			newVert->markAsShellVert(true);
			m_isoLength->setLengthScale(newVert);

		}
	}

}

int SurfaceInserter2D::removeVertsAfterInsertion() {
	int numRemoved = 0;
	for (std::set<Vert*>::iterator it = m_leftoverVertsToRemove.begin();
			it != m_leftoverVertsToRemove.end(); ++it) {
		int nSwaps;
		if (m_mesh->removeVert(*it, nSwaps))
			numRemoved++;
		else
			logMessage(1, "Failed to remove vert %f %f\n", (*it)->x(),
					(*it)->y());
	}
	return numRemoved;
}
// End of class

void SurfaceInserter2D::setVertLengthScaleMetrics() {

//	for(GR_index_t iV = 0; iV < m_mesh->getNumVerts(); iV++){
//		std::set<Cell*> spCInc;
//		std::set<Vert*> spVTemp;
//		findNeighborhoodInfo(m_mesh->getVert(iV),spCInc,spVTemp);
//		double metric[] = {0.0,0.0,0.0,0.0,0.0};
//
//		for(std::set<Cell*>::iterator itc = spCInc.begin(); itc != spCInc.end(); ++itc){
//			Cell * cell = *itc;
//			double tempMetric[3] = {0.0,0.0,0.0};
//			// use the normal metric, by swapping diagonal elements of the covariance metric,
//			// then unswap when need to measure length
//			calcCovarianceMatrix(cell->getVert(0),cell->getVert(1),cell->getVert(2),tempMetric);
//			// find longest edge in physical space, align metric with it
//			double scaling = 18./cell->calcSize();
//			for(int i = 0; i < 3; i++)
//				metric[i] += scaling*tempMetric[i];
//		}
//
//		for(int i = 0; i < 3; i++)
//			metric[i] *= 1./spCInc.size();
//		logMessage(4,"Aniso Metric for Vert %d | %f %f %f\n",iV,metric[0],metric[1],metric[2]);
//		m_mesh->getVert(iV)->setMetric(metric);
//	}

	m_anisoLength->setLengthScale(m_mesh);
}
void SurfaceInserter2D::setVertOriginalMetrics() {
	for (GR_index_t iV = 0; iV < m_mesh->getNumVerts(); iV++) {
		Vert * vert = m_mesh->getVert(iV);
		std::set<Cell*> spCInc;
		std::set<Vert*> spVTemp;
		findNeighborhoodInfo(vert, spCInc, spVTemp);

		double metric[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

		for (std::set<Cell*>::iterator itc = spCInc.begin();
				itc != spCInc.end(); ++itc) {
			Cell * cell = *itc;
			double tempMetric[3] = { 0.0, 0.0, 0.0 };
			// use the normal metric, by swapping diagonal elements of the covariance metric,
			// then unswap when need to measure length
			calcCovarianceMatrix(cell->getVert(0), cell->getVert(1),
					cell->getVert(2), tempMetric);
			double scaling = 18. / cell->calcSize() / cell->calcSize()
					/ cell->calcSize();

			for (int i = 0; i < 3; i++)
				metric[i] += scaling * tempMetric[2 - i];
		}
		for (int i = 0; i < 3; i++)
			metric[i] *= 1. / spCInc.size();

		metric[1] *= -1.0;

		vert->setMetric(metric);
	}

	//m_LengthAniso->setLengthScale(m_mesh);
}

int SurfaceInserter2D::sampleQuasiStructured() {
	//	// determine all intersection points with mesh
	std::list<GRCurve*> curve_list;
	std::list<GRCurve*>::iterator itc;
	m_geom->get_curves(curve_list);
	logMessage(1, "Quasi Structured Intersections\n");

	// use of CubitVectors due to original geometry implementation
	CubitVector tempCubitVector, cvToInsert;
	int iCurveNum = 0;
	/// sample each curve individually
	Cell * seedCell = m_mesh->getCell(0);
	CubitVector cv;
	double coord[3];
	double param = 0.;
	int numViolated = 0, numTies = 0;
	Face *facesViolated[3], *facesTied[3];
	for (itc = curve_list.begin(); itc != curve_list.end();
			itc++, iCurveNum++) {
		GRCurve * curveToInsert = (*itc);
		insPointSet * pointSet = new insPointSet();
		// find faces that intersect by traversal.
		GRCurveGeom * cg = curveToInsert->get_curve_geom();
		param = cg->min_param();
		cg->coord_at_param(param, cv);
		cv.get_xyz(coord);
		seedCell = m_mesh->findCell(coord, seedCell, numViolated, facesViolated,
				numTies, facesTied);

		// this is the first cell.
		// check its three faces for intersection
		std::vector<double> tempParams;

		pointSet->insert(insPoint2D(cg->min_param(), seedCell));
		for (GR_index_t iC = 0; iC < m_mesh->getNumCells(); iC++) {
			Cell * cell = m_mesh->getCell(iC);
			tempParams.clear();
			GR_index_t startIndex = 0;
			std::map<double, double> paramSet;

			for (GR_index_t iF = 0; iF < 3; iF++) {
				Face * face = cell->getFace(iF);
				CubitVector facedir = ~(CubitVector(
						face->getVert(0)->getCoords())
						- CubitVector(face->getVert(1)->getCoords()));
				CubitVector normal;
				cg->line_intersect(face->getVert(0)->getCoords(),
						face->getVert(1)->getCoords(), cg->min_param(),
						cg->max_param(), tempParams);

				for (GR_index_t i = startIndex; i < tempParams.size(); i++) {
					cg->unit_normal(tempParams[i], normal);
					double dot = fabs(facedir % normal);
					//				printf("temp %d\n",tempParams.size());
					if (dot > 1 / 2.) {
						cg->coord_at_param(tempParams[i], cv);
						cv.get_xyz(coord);
						//						printf("coord set\n");
						//						printf("dot %f param %f | %f %f \n",dot,tempParams[i],coord[0],coord[1]);
						//						printf("normal %f %f | %f %f\n",normal.x(),normal.y(),facedir.x(),facedir.y());
						//						face->getVert(0)->printVertInfo();
						//						face->getVert(1)->printVertInfo();

						//						assert(cg->coord_on_curve(cv));
						//						seedCell = m_mesh->findCell(coord,seedCell,numViolated,facesViolated,numTies,facesTied);
						paramSet.insert(
								std::pair<double, double>(dot, tempParams[i]));

					}

				}
				startIndex = tempParams.size();
				// grab the maximal one first.

			}
			if (!paramSet.empty()) {
				//					for(std::map<double,double>::iterator itd = paramSet.begin(); itd != paramSet.end(); ++itd){
				//						printf("pointset %.16f param %.16f\n",itd->first,itd->second);
				//					}
				//					printf("%f %f winner\n",paramSet.rbegin()->first,paramSet.rbegin()->second);
				pointSet->insert(insPoint2D(paramSet.rbegin()->second, cell));
			}

		}
		//		int ii = 0;
		//		// lets do this a different way. Find a seed cell.
		//		while((!qDone || seedCell->getType() != Cell::eTriCell) && ii++ < 20){
		//			std::set<Face*> neighFaces;
		//			for(GR_index_t iV = 0; iV < 3; iV++){
		//				std::set<Cell*> spCInc;
		//				std::set<Vert*> spVTemp;
		//				std::set<BFace*> pspBFInc;
		//				std::set<Face*> faces;
		//
		//				bool qBdryVert = true;
		//				findNeighborhoodInfo(seedCell->getVert(iV),spCInc,spVTemp,&pspBFInc,&qBdryVert,&faces);
		//				neighFaces.insert(faces.begin(),faces.end());
		//
		//			}
		//			// find a cell with two faces going through it
		//
		//			for(std::set<Face*>::iterator it = neighFaces.begin(); it != neighFaces.end(); ++it){
		//				Face * face = *it;
		//				tempParams.clear();
		//				int intersect = cg->line_intersect(face->getVert(0)->getCoords(),face->getVert(1)->getCoords(),tempParams);
		//				for(GR_index_t i = 0; i < tempParams.size(); i++){
		//					if(tempParams[i] > param){
		//						seedFace = face;
		//						seedCell = face->getCell(0) == seedCell ? face->getCell(0) : face->getCell(1);
		//						printf("%.16f tempParam\n",tempParams[i]);
		//						param = tempParams[i];
		//
		//					}
		//
		//				}
		//			}
		//		}

		//		for(std::set<double>::iterator itd = params.begin(); itd != params.end(); ++itd){
		//			printf("inserting %.16f param\n",*itd);
		//			pointSet->insert(insPoint(*itd,seedCell));
		//
		//		}
		pointSet->insert(insPoint2D(cg->max_param(), seedCell));
		//		for(insPointSet_it it = pointSet->begin(); it != pointSet->end(); ++it ){
		//			printf("%.16f param\n",it->getParam());
		//		}
		m_pointsToInsert.insert(
				std::pair<GRCurve*, insPointSet*>(curveToInsert, pointSet));

	}
	std::vector<CubitVector> tempCVs;
	//	CubitVector cv;
	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); cpm_it++) {
		GRCurve * curveToInsert = cpm_it->first;
		insPointSet * pointSet = cpm_it->second;
		insPointSet_it pointSet_end = pointSet->end();
		insPointSet_it pointSet_it = pointSet->begin();
		for (; pointSet_it != pointSet_end; ++pointSet_it) {
			curveToInsert->position_from_u(pointSet_it->getParam(), cv);
			tempCVs.push_back(cv);
		}
	}
	writeVTKPointSet(tempCVs, NULL, "surfins_samplepts_init", "");
	tempCVs.clear();

	return m_pointsToInsert.size();
}
//	// determine all intersection points with mesh
//	std::list<GRCurve*> curve_list;
//	std::list<GRCurve*>::iterator itc;
//	m_geom->get_curves(curve_list);
//
//	logMessage(1,"Quasi Structured Intersections\n");
//
//	// use of CubitVectors due to original geometry implementation
//	CubitVector tempCubitVector, cvToInsert;
//	int iCurveNum = 0;
//	/// sample each curve individually
//	Cell * seedCell = m_mesh->getCell(0);
//	CubitVector cv;
//	double coord[3];
//	double param = 0.;
//	int numViolated = 0, numTies = 0;
//	Face *facesViolated[3], *facesTied[3];
//	Face * seedFace = NULL;
//	for(itc = curve_list.begin(); itc != curve_list.end(); itc++, iCurveNum++) {
//		GRCurve * curveToInsert = (*itc);
//
//		insPointSet * pointSet = new insPointSet();
//		// find faces that intersect by traversal.
//		GRCurveGeom * cg = curveToInsert->get_curve_geom();
//		param = cg->min_param();
//		cg->coord_at_param(param,cv);
//		cv.get_xyz(coord);
//		seedCell = m_mesh->findCell(coord,seedCell,numViolated,facesViolated,numTies,facesTied);
//
//
//		// this is the first cell.
//		// check its three faces for intersection
//		std::vector<double> tempParams;
//		bool qDone = false, qDiagonal = true;
//		Vert * seedVert = NULL;
//		double testParam = 0.0;
//		pointSet->insert(insPoint(cg->min_param(),seedCell));
//		int ii = 0;
//		while((!qDone || seedCell->getType() != Cell::eTriCell) && ii < 20){
//			ii++;
//			if(!qDiagonal){
//				pointSet->insert(insPoint(testParam,seedCell));
//			}
//			tempParams.clear();
//
//			double longest = 0.0;
//			for(int iF = 0; iF < 3; iF++){
//				Face * face = seedCell->getFace(iF);
//				longest = max(longest,face->calcSize());
//			}
//			testParam = 10.0*cg->max_param();
//			// go through every face, find the lowest param thats greater than the current param
//			for(int iF = 0; iF < 3; iF++){
//				Face * face = seedCell->getFace(iF);
////				if(face->isBdryFace()) continue;
//				tempParams.clear();
//				int intersect = cg->line_intersect(face->getVert(0)->getCoords(),face->getVert(1)->getCoords(),tempParams);
//				for(GR_index_t i = 0; i < tempParams.size(); i++){
////					printf("%d %d\n",testParam > tempParams[i], iFuzzyComp(param,tempParams[i]));
////					printf("%f %f\n",tempParams[i],param);
//
//					if(testParam > tempParams[i] && iFuzzyComp(param,tempParams[i]) == -1 ){
//						if(intersect == 0 && intersect == 1 && seedVert != face->getVert(intersect)){
//							seedVert = face->getVert(intersect);
//						} else {
//						testParam = tempParams[i];
//						seedFace = face;
//						}
//					}
//				}
//			}
//			param = testParam;
////			printf("testParam %f\n",testParam);
//			if(!seedFace->isBdryFace())
//			seedCell = seedFace->getOppositeCell(seedCell);
//
////				printf("intersect %d %d\n",intersect,tempParams.size());
////				if( intersect >= 2 ){
////
////					seedCell = face->getOppositeCell(seedCell);
////					break;
////				} else
////					if( intersect > -1 && seedVert != face->getVert(intersect)){
////
////						seedVert = face->getVert(intersect);
////						//					for(GR_index_t i = 0; i < tempParams.size(); i++){
////						//printf("tempParams[ %.16f\n",tempParams[i]);
////						//					}
////						double maxParam = 0.0;
////						for(GR_index_t i = 0; i < tempParams.size(); i++){
////							maxParam = max(maxParam,tempParams[i]);
////						}
////						maxParam += 0.00001;
////						if(maxParam < cg->max_param()){
////							cg->coord_at_param(maxParam,cv);
////							cv.get_xyz(coord);
////							seedCell = m_mesh->findCell(coord,seedCell,numViolated,facesViolated,numTies,facesTied);
////						}
////						else
////							tempParams.clear();
////
////					}
////			}
//			if(testParam > 5.0*cg->max_param()){
//					qDone = true;
//			}
//			if(iFuzzyComp(seedFace->calcSize(),longest) == 0){
//				qDiagonal = true;
//			} else
//				qDiagonal = false;
//
//		}
//
//		pointSet->insert(insPoint(cg->max_param(),seedCell));
//		for(insPointSet_it it = pointSet->begin(); it != pointSet->end(); ++it ){
//			printf("%.16f param\n",it->getParam());
//		}
//		m_pointsToInsert.insert(std::pair<GRCurve*,insPointSet*>(curveToInsert,pointSet));
//
//	}
//	return m_pointsToInsert.size();
//}
void SurfaceInserter2D::createOffsetQuadCells() {
	logMessage(1, "Creating quads between offset curves\n");
	m_mesh->allowNonSimplicial();
	//need to find pairs of surfaces. I know my pairs
	for (std::multimap<GRCurve*, GRCurve*>::iterator curveMap_it =
			m_offsetCurveMap.begin(); curveMap_it != m_offsetCurveMap.end();
			++curveMap_it) {
		insertionCurveMap_it cpm_it1 = m_pointsToInsert.find(
				curveMap_it->first);
		insertionCurveMap_it cpm_it2 = m_pointsToInsert.find(
				curveMap_it->second);

		assert(cpm_it1 != m_pointsToInsert.end());
		assert(cpm_it2 != m_pointsToInsert.end());

		// iterate across pairs to figure out what quads to create

		insPointSet * pointSet1 = cpm_it1->second;
		insPointSet * pointSet2 = cpm_it2->second;

		// they must match up
		assert(pointSet1->size() == pointSet2->size());
		insPointSet_it pointSet_it1 = pointSet1->begin();
		insPointSet_it pointSet_it2 = pointSet2->begin();

		// go until curve 1 hits the end until the end
		insPointSet_it pointSet_end1 = --pointSet1->end();
		// the first set of points
		insPoint2D pt1_1 = (*pointSet_it1), pt2_1 = (*pointSet_it2);
		// second set
		insPoint2D pt1_2, pt2_2;

		// iterate over two pairs, recover edge
		for (; pointSet_it1 != pointSet_end1;) {
			// the second set of points
			pt1_2 = *(++pointSet_it1);
			pt2_2 = *(++pointSet_it2);
			m_mesh->createQuadFromVerts(pt1_1.getVert(), pt1_2.getVert(),
					pt2_2.getVert(), pt2_1.getVert());
			// update points
			pt1_1 = pt1_2;
			pt2_1 = pt2_2;
		}

	}

}
void SurfaceInserter2D::writeSamplePointsToFile(const char * filename) {
	CubitVector point;
	std::vector<CubitVector> pointVector;

	for (insertionCurveMap_it cpm_it = m_pointsToInsert.begin();
			cpm_it != m_pointsToInsert.end(); cpm_it++) {
		GRCurve * curveToInsert = cpm_it->first;
		insPointSet * pointSet = cpm_it->second;
		insPointSet_it pointSet_end = pointSet->end();
		insPointSet_it pointSet_it = pointSet->begin();
		for (; pointSet_it != pointSet_end; ++pointSet_it) {
			curveToInsert->position_from_u(pointSet_it->getParam(), point);
			pointVector.push_back(point);
		}
	}
	writeVTKPointSet(pointVector, NULL, filename, "");
	pointVector.clear();
}

}		// end of namespace

void smoothBdryRefine(Mesh2D* pM2D, GRGeom2D* pG2D, double dBGrading,
		double dBScale, char strSurfInsParams[]) {
	std::shared_ptr<GRUMMP::Length> len2d = GRUMMP::Length::Create(pM2D,
			dBScale, dBGrading);
	GRUMMP::SurfaceInserter2D SurfIns2D(pM2D, len2d, dBGrading, dBScale,
			strSurfInsParams);
	//			GRUMMP::SurfaceInserter2D SurfIns2D(pM2D, dGrading, dScale, strSurfInsParams);
	SurfIns2D.refineBoundaryForLength(pG2D);
}

