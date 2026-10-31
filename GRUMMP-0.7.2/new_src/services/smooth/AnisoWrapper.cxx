/* ***************************************************************** 
 MESQUITE -- The Mesh Quality Improvement Toolkit

 Copyright 2010 Sandia National Laboratories.  Developed at the
 University of Wisconsin--Madison under SNL contract number
 624796.  The U.S. Government and the University of Wisconsin
 retain certain rights to this software.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License 
 (lgpl.txt) along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 (2010) kraftche@cae.wisc.edu    

 ***************************************************************** */

#include "GR_config.h"
#ifdef HAVE_MESQUITE

/** \file AnisoWrapper.cpp
 *  \brief Impement AnisoWrapper class
 *  \author Jason Kraftcheck 
 */

#include "AnisoWrapper.hpp"
#include "MsqDebug.hpp"
#include "LInfTemplate.hpp"
#include "LPtoPTemplate.hpp"
#include "PMeanPTemplate.hpp"
#include "TrustRegion.hpp"
#include "ConjugateGradient.hpp"
#include "TerminationCriterion.hpp"
#include "InstructionQueue.hpp"
#include "QualityAssessor.hpp"

#include "TPower2.hpp"
#include "TShapeB1.hpp"
#include "TShapeNB1.hpp"
#include "TShape2DNB2.hpp"
#include "TShapeSize2DNB1.hpp"
#include "TShapeSizeB1.hpp"
#include "TShapeOrientNB1.hpp"
#include "TShapeOrientB1.hpp"
#include "TShapeSizeOrientB1.hpp"
#include "TShapeSizeOrientNB1.hpp"

#include "TQualityMetric.hpp"
#include "IdealShapeTarget.hpp"
#include "TagVertexMesh.hpp"
#include "ReferenceMesh.hpp"

#include "TargetWriter.hpp"
#include "Settings.hpp"

#include <iostream>

namespace MESQUITE_NS
{

  const double DEFAULT_BETA = 0.005;
  const double DEFAULT_BETA_3D = 0.005;
  const int DEFUALT_PARALLEL_ITERATIONS = 10;

  /**\brief Targets combining ideal skew with reference aspect ratio */
  class RefAspect : public TargetCalculator {
  public:
    explicit
    RefAspect(ReferenceMeshInterface* ref_mesh) :
	refMesh(ref_mesh)
    {
    }

    virtual
    ~RefAspect()
    {
    }

    virtual bool
    get_3D_target(PatchData& pd, size_t element, Sample sample,
		  MsqMatrix<3, 3>& W_out, MsqError& err);

    virtual bool
    get_surface_target(PatchData& pd, size_t element, Sample sample,
		       MsqMatrix<3, 2>& W_out, MsqError& err);

    virtual bool
    get_2D_target(PatchData& pd, size_t element, Sample sample,
		  MsqMatrix<2, 2>& W_out, MsqError& err);

    virtual bool
    have_surface_orient() const
    {
      return false;
    }

  private:
    ReferenceMeshInterface* refMesh;
  };

  bool
  RefAspect::get_3D_target(PatchData& pd, size_t element, Sample /*sample*/,
			   MsqMatrix<3, 3>& W_out, MsqError& err)
  {
//  MsqMatrix<3,3> ref_J;
//  get_refmesh_Jacobian_3D( refMesh, pd, element, sample, ref_J, err ); MSQ_ERRZERO(err);
//  ideal_skew_3D( pd.element_by_index(element).get_element_type(),
//                 sample, pd, W_out, err ); MSQ_ERRZERO(err);
//  W_out = W_out * aspect( ref_J );

    std::vector<Mesquite::Mesh::ElementHandle> element_handle;

    Mesh *msqMesh = pd.get_mesh();
    msqMesh->get_all_elements(element_handle, err);

    TagHandle lid_tag1 = msqMesh->tag_get("LOCAL_ID1", err);
    TagHandle lid_tag2 = msqMesh->tag_get("LOCAL_ID2", err);
    TagHandle lid_tag3 = msqMesh->tag_get("LOCAL_ID3", err);
    TagHandle lid_tag4 = msqMesh->tag_get("LOCAL_ID4", err);
    TagHandle lid_tag5 = msqMesh->tag_get("LOCAL_ID5", err);
    TagHandle lid_tag6 = msqMesh->tag_get("LOCAL_ID6", err);
    TagHandle lid_tag7 = msqMesh->tag_get("LOCAL_ID7", err);
    TagHandle lid_tag8 = msqMesh->tag_get("LOCAL_ID8", err);
    TagHandle lid_tag9 = msqMesh->tag_get("LOCAL_ID9", err);
    std::vector<double> lidtest1(1);
    std::vector<double> lidtest2(1);
    std::vector<double> lidtest3(1);
    std::vector<double> lidtest4(1);
    std::vector<double> lidtest5(1);
    std::vector<double> lidtest6(1);
    std::vector<double> lidtest7(1);
    std::vector<double> lidtest8(1);
    std::vector<double> lidtest9(1);

    msqMesh->tag_get_element_data(lid_tag1, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest1), err);
    msqMesh->tag_get_element_data(lid_tag2, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest2), err);
    msqMesh->tag_get_element_data(lid_tag3, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest3), err);
    msqMesh->tag_get_element_data(lid_tag4, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest4), err);
    msqMesh->tag_get_element_data(lid_tag5, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest5), err);
    msqMesh->tag_get_element_data(lid_tag6, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest6), err);
    msqMesh->tag_get_element_data(lid_tag7, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest7), err);
    msqMesh->tag_get_element_data(lid_tag8, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest8), err);
    msqMesh->tag_get_element_data(lid_tag9, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest9), err);
    W_out(0, 0) = lidtest1[0];
    W_out(0, 1) = lidtest2[0];
    W_out(0, 2) = lidtest3[0];
    W_out(1, 0) = lidtest4[0];
    W_out(1, 1) = lidtest5[0];
    W_out(1, 2) = lidtest6[0];
    W_out(2, 0) = lidtest7[0];
    W_out(2, 1) = lidtest8[0];
    W_out(2, 2) = lidtest9[0];

    return true;
  }

  bool
  RefAspect::get_surface_target(PatchData& /*pd*/, size_t /*element*/,
				Sample /*sample*/, MsqMatrix<3, 2>& /*W_out*/,
				MsqError& err)
  {
    MSQ_SETERR(err)("Internal error requesting 3x2 surface target for non-oriented targets",
    MsqError::INVALID_STATE);
    return false;
  }
  bool
  RefAspect::get_2D_target(PatchData& pd, size_t element, Sample /*sample*/,
			   MsqMatrix<2, 2>& W_out, MsqError& err)
  {
    std::vector<Vector3D> Vcoords;
    // pd.get_element_vertex_coordinates(element,Vcoords,err);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // MsqMeshEntity element_Entity=pd.element_by_index(element);
    std::vector<Mesquite::Mesh::ElementHandle> element_handle;

    // const std::size_t *Vindex=element_Entity.get_vertex_index_array();
    Mesh *msqMesh = pd.get_mesh();
    msqMesh->get_all_elements(element_handle, err);
    // std::vector<size_t> vtx_offsets;
    //  msqMesh->elements_get_attached_vertices(&element_handle[element],1,vertices,vtx_offsets,err);
    TagHandle lid_tag1 = msqMesh->tag_get("LOCAL_ID1", err);
    TagHandle lid_tag2 = msqMesh->tag_get("LOCAL_ID2", err);
    TagHandle lid_tag3 = msqMesh->tag_get("LOCAL_ID3", err);
    TagHandle lid_tag4 = msqMesh->tag_get("LOCAL_ID4", err);
    std::vector<double> lidtest1(1);
    std::vector<double> lidtest2(1);
    std::vector<double> lidtest3(1);
    std::vector<double> lidtest4(1);

    msqMesh->tag_get_element_data(lid_tag1, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest1), err);
    msqMesh->tag_get_element_data(lid_tag2, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest2), err);
    msqMesh->tag_get_element_data(lid_tag3, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest3), err);
    msqMesh->tag_get_element_data(lid_tag4, 1, &element_handle[element],
				  Mesquite::arrptr(lidtest4), err);
    W_out(0, 0) = lidtest1[0];
    W_out(0, 1) = lidtest2[0];
    W_out(1, 0) = lidtest3[0];
    W_out(1, 1) = lidtest4[0];

//    msqMesh->tag_get_vertex_data(lid_tag1,3,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest1),err);
//    msqMesh->tag_get_vertex_data(lid_tag2,3,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest2),err);
//    msqMesh->tag_get_vertex_data(lid_tag3,3,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest3),err);
//    msqMesh->tag_get_vertex_data(lid_tag4,3,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest4),err);
//        W_out(0,0)=(lidtest1[0]+lidtest1[1]+lidtest1[2])/3;
//	W_out(0,1)=(lidtest2[0]+lidtest2[1]+lidtest2[2])/3;
//	W_out(1,0)=(lidtest3[0]+lidtest3[1]+lidtest3[2])/3;
//	W_out(1,1)=(lidtest4[0]+lidtest4[1]+lidtest4[2])/3;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //msqMesh->tag_get_vertex_data( lid_tag1,num_verts,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest1),MesqErr);
    //msqMesh->tag_get_vertex_data( lid_tag2,num_verts,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest2),MesqErr);
    //msqMesh->tag_get_vertex_data( lid_tag3,num_verts,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest3),MesqErr);

//   if(Vcoords[0].y()<25 && Vcoords[0].y()>-25 )
// {
//	   W_out(0,0) = 1; W_out(0,1) =0;
//	   W_out(1,0) =0; W_out(1,1) = 10;
//  }
//   else
//   {
//	   W_out(0,0) = 1; W_out(0,1) =0.5;
//	   W_out(1,0) =0; W_out(1,1) = sqrt(3)/2;
//   }
    return true;
  }

//bool RefAspect::get_2D_target( PatchData& pd,
//                          size_t element,
//                          Sample sample,
//                          MsqMatrix<2,2>& W_out,
//                          MsqError& err )
//{
//	std::vector<Vector3D> Vcoords;
//   // pd.get_element_vertex_coordinates(element,Vcoords,err);
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////
//    MsqMeshEntity element_Entity=pd.element_by_index(element);
//    std::vector<Mesquite::Mesh::ElementHandle> element_handle;
//    const std::size_t *Vindex=element_Entity.get_vertex_index_array();
//    Mesh *msqMesh=pd.get_mesh();
//    msqMesh->get_all_elements(element_handle,err);
//    std::vector<Mesquite::Mesh::VertexHandle>vertices;
//    std::vector<size_t> vtx_offsets;
//    msqMesh->elements_get_attached_vertices(&element_handle[element],1,vertices,vtx_offsets,err);
//    TagHandle lid_tag1=msqMesh->tag_get("LOCAL_ID1",err);
//    TagHandle lid_tag2=msqMesh->tag_get("LOCAL_ID2",err);
//    TagHandle lid_tag3=msqMesh->tag_get("LOCAL_ID3",err);
//    TagHandle lid_tag4=msqMesh->tag_get("LOCAL_ID4",err);
//    std::vector<double> lidtest1(3);
//    std::vector<double> lidtest2(3);
//    std::vector<double> lidtest3(3);
//    std::vector<double> lidtest4(3);
//    msqMesh->tag_get_vertex_data(lid_tag1,3,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest1),err);
//    msqMesh->tag_get_vertex_data(lid_tag2,3,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest2),err);
//    msqMesh->tag_get_vertex_data(lid_tag3,3,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest3),err);
//    msqMesh->tag_get_vertex_data(lid_tag4,3,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest4),err);
//        W_out(0,0)=(lidtest1[0]+lidtest1[1]+lidtest1[2])/3;
//	W_out(0,1)=(lidtest2[0]+lidtest2[1]+lidtest2[2])/3;
//	W_out(1,0)=(lidtest3[0]+lidtest3[1]+lidtest3[2])/3;
//	W_out(1,1)=(lidtest4[0]+lidtest4[1]+lidtest4[2])/3;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    //msqMesh->tag_get_vertex_data( lid_tag1,num_verts,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest1),MesqErr);
//    //msqMesh->tag_get_vertex_data( lid_tag2,num_verts,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest2),MesqErr);
//    //msqMesh->tag_get_vertex_data( lid_tag3,num_verts,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest3),MesqErr);
//
//
////   if(Vcoords[0].y()<25 && Vcoords[0].y()>-25 )
//// {
////	   W_out(0,0) = 10; W_out(0,1) =0;
////	   W_out(1,0) =0; W_out(1,1) = 1;
////  }
////   else
////   {
////	   W_out(0,0) = 1; W_out(0,1) =0.5;
////	   W_out(1,0) =0; W_out(1,1) = sqrt(3)/2;
////   }
//  return true;
//}

  /*
   * For the local smoothing function of SmoothingManagerMesquite2D, the maxTime must be set really low (1 or 2)
   */
  AnisoWrapper::AnisoWrapper() :
      maxTime(100), mBeta(DEFAULT_BETA), parallelIterations(
	  DEFUALT_PARALLEL_ITERATIONS)
  {
  }

  AnisoWrapper::~AnisoWrapper()
  {
  }

  void
  AnisoWrapper::set_vertex_movement_limit_factor(double beta)
  {
    assert(beta > 0.0);
    assert(beta < 1.0);
    mBeta = beta;
  }

  void
  AnisoWrapper::set_cpu_time_limit(double limit)
  {
    assert(limit >= 0.0);
    maxTime = limit;
  }

  void
  AnisoWrapper::set_parallel_iterations(int count)
  {
    if (count < 1) {
      assert(false);
      count = 1;
    }
    parallelIterations = count;
  }

//void AnisoWrapper::run_wrapper( Mesh* mesh,
//                                ParallelMesh* pmesh,
//                                MeshDomain* domain,
//                                Settings* settings,
//                                QualityAssessor* qa,
//                                MsqError& err )
//{
//    // Use initial vertex coordinates as reference mesh
//    // NOTE: need to put init_mesh in instruction queue
//    //       to get it initialized with initial coords.
//  TagVertexMesh init_mesh( err, mesh ); MSQ_ERRRTN(err);
//  ReferenceMesh ref_mesh( &init_mesh );
//    // Quality Metrics
//  RefAspect target( &ref_mesh );
//    // No-barrier phase
//  TShape2DNB2 mu_no;
//  TPower2 mu_tp2(&mu_no);
//  TQualityMetric metric_no( &target, &mu_no );
//    // Barrier phase
//  TShapeB1 mu_b;
//  TQualityMetric metric_b( &target, &mu_b );
//    // QualityAssessor
//  qa->add_quality_assessment( &metric_no );
//  qa->add_quality_assessment( &metric_b );
//
//    // Phase 1: use non-barrier metric until mesh does not
//    //          have any inverted elements
//  PMeanPTemplate obj_func_no( 1.0, &metric_no );
//ConjugateGradient improver_no( &obj_func_no );
////TrustRegion improver_no( &obj_func_no );
//  improver_no.use_global_patch();
//  TerminationCriterion inner_no;
//  inner_no.add_absolute_vertex_movement_edge_length( mBeta );
//  inner_no.add_untangled_mesh();
//  if (maxTime > 0.0)
//    inner_no.add_cpu_time( maxTime );
//  improver_no.set_inner_termination_criterion( &inner_no );
//  InstructionQueue q_no;
//  Timer totalTimer;
//  q_no.add_tag_vertex_mesh( &init_mesh, err ); MSQ_ERRRTN(err);
//  q_no.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
//  q_no.set_master_quality_improver( &improver_no, err ); MSQ_ERRRTN(err);
//  q_no.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
//  q_no.run_common( mesh, pmesh, domain, settings, err ); MSQ_ERRRTN(err);
//    // Phase 2: use-barrer metric on untangled mesh
//  LPtoPTemplate obj_func_b( 1, &metric_b );
////   PMeanPTemplate obj_func_b( 1.0, &metric_b );
//ConjugateGradient improver_b( &obj_func_b );
// //TrustRegion improver_b( &obj_func_b );
//  improver_b.use_global_patch();
////improver_no.use_element_on_vertex_patch();
//  TerminationCriterion inner_b;
// inner_b.add_absolute_gradient_L2_norm( 5e-7);
// //  inner_b.add_relative_vertex_movement (0.001);
///*
//add_absolute_gradient_L2_norm( double eps )
//add_absolute_gradient_inf_norm( double eps )
//add_relative_gradient_L2_norm( double eps )
//add_relative_gradient_inf_norm( double eps )
//add_absolute_quality_improvement( double eps )
//add_relative_quality_improvement( double eps )
//add_absolute_vertex_movement( double eps )
//add_absolute_vertex_movement_edge_length( double beta )
//add_relative_vertex_movement( double eps )
//add_absolute_successive_improvement( double eps )
//add_relative_successive_improvement( double eps )*/
// ////////////////////////////optimization criteria//////////////////////////////
// //inner_b.add_absolute_gradient_L2_norm( 5e-3);
// if (maxTime > 0.0) {
//    double remaining = (maxTime - totalTimer.since_birth());
//    if (remaining <= 0.0 ){
//      //MSQ_DBGOUT(2) << "Optimization is terminating without perfoming shape improvement." << std::endl;
//      remaining = 0.0;
//    }
//    inner_b.add_cpu_time( remaining );
//  }
//  improver_b.set_inner_termination_criterion( &inner_b );
//  InstructionQueue q_b;
//  // NOTE: Don't add init_mesh to this queue because we to use
//  //       the initial initial mesh (not the output of the non-barrier
//  //       phase).
//  q_b.set_master_quality_improver( &improver_b, err ); MSQ_ERRRTN(err);
//  q_b.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
//      //MsqMeshEntity element_Entity=pd.element_by_index(element);
//      //const std::size_t *Vindex=element_Entity.get_vertex_index_array();
//      //Mesh *msqMesh=pd.get_mesh();
//      std::vector<Mesquite::Mesh::VertexHandle>vertices;
//      mesh->get_all_vertices(vertices, err);
//      int num_verts = vertices.size();
//      TagHandle lid_tag1=mesh->tag_get("LOCAL_ID1",err);
//      TagHandle lid_tag2=mesh->tag_get("LOCAL_ID2",err);
//      TagHandle lid_tag3=mesh->tag_get("LOCAL_ID3",err);
//      TagHandle lid_tag4=mesh->tag_get("LOCAL_ID4",err);
//      std::vector<double> lidtest1(num_verts);
//      std::vector<double> lidtest2(num_verts);
//      std::vector<double> lidtest3(num_verts);
//      std::vector<double> lidtest4(num_verts);
//      mesh->tag_get_vertex_data(lid_tag1,num_verts,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest1),err);
//      mesh->tag_get_vertex_data(lid_tag2,num_verts,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest2),err);
//      mesh->tag_get_vertex_data(lid_tag3,num_verts,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest3),err);
//      mesh->tag_get_vertex_data(lid_tag4,num_verts,Mesquite::arrptr(vertices),Mesquite::arrptr(lidtest4),err);
//      q_b.run_common( mesh, pmesh, domain, settings, err ); MSQ_ERRRTN(err);
//
//  TargetWriter writer( &target );
//  Settings linearMaps;
//  writer.loop_over_mesh( mesh, domain, &linearMaps, err );
//}
  void
  AnisoWrapper::run_wrapper(Mesh* mesh, ParallelMesh* pmesh, MeshDomain* domain,
			    Settings* settings, QualityAssessor* qa,
			    MsqError& err)
  {
    // Use initial vertex coordinates as reference mesh
    // NOTE: need to put init_mesh in instruction queue
    //       to get it initialized with initial coords.
    TagVertexMesh init_mesh(err, mesh);
    MSQ_ERRRTN(err);
    ReferenceMesh ref_mesh(&init_mesh);

    // Quality Metrics
    RefAspect target(&ref_mesh);

    // No-barrier phase
    TShape2DNB2 mu_no;
    TPower2 mu_tp2(&mu_no);
    TQualityMetric metric_no(&target, &mu_no);

    // Barrier phase
    TShapeB1 mu_b;
    TShapeOrientB1 mu_b2;
    TShapeSizeB1 mu_b3;
    TShapeSizeOrientB1 mu_b4;
    TQualityMetric metric_b(&target, &mu_b);

// if decide to use an ideal isotorpic target
    IdealShapeTarget Itarget;
    TQualityMetric metric_bIdeal(&Itarget, &mu_b);

    // QualityAssessor
    qa->add_quality_assessment(&metric_no);
    qa->add_quality_assessment(&metric_b);

    // Phase 1: use non-barrier metric until mesh does not
    //          have any inverted elements
    PMeanPTemplate obj_func_no(1.0, &metric_no);
    ConjugateGradient improver_no(&obj_func_no);
//TrustRegion improver_no( &obj_func_no );
    //improver_no.use_global_patch();
    TerminationCriterion inner_no;
    inner_no.add_absolute_vertex_movement_edge_length(mBeta);
    inner_no.add_untangled_mesh();
    if (maxTime > 0.0)
      inner_no.add_cpu_time(maxTime);
    improver_no.set_inner_termination_criterion(&inner_no);
    InstructionQueue q_no;
    Timer totalTimer;
    q_no.add_tag_vertex_mesh(&init_mesh, err);
    MSQ_ERRRTN(err);
    q_no.add_quality_assessor(qa, err);
    MSQ_ERRRTN(err);
    q_no.set_master_quality_improver(&improver_no, err);
    MSQ_ERRRTN(err);
    q_no.add_quality_assessor(qa, err);
    MSQ_ERRRTN(err);
    q_no.run_common(mesh, pmesh, domain, settings, err);
    MSQ_ERRRTN(err);

    // Phase 2: use-barrer metric on untangled mesh
    LPtoPTemplate obj_func_b(1, &metric_b);
//   PMeanPTemplate obj_func_b( 1.0, &metric_b );
//ConjugateGradient improver_b( &obj_func_b );
    TrustRegion improver_b(&obj_func_b);
    improver_b.use_global_patch();
//  improver_b.use_element_on_vertex_patch();
    TerminationCriterion inner_b;
//   inner_b.add_absolute_gradient_L2_norm(5e-5);
    inner_b.add_relative_gradient_L2_norm(0.01);
    /*
     //   inner_b.add_iteration_limit(5);
     //   inner_b.add_absolute_vertex_movement_edge_length( 0.0000001 );
     //   inner_b.add_relative_vertex_movement( 0.000001 );
     //add_absolute_gradient_L2_norm( double eps )
     //add_absolute_gradient_inf_norm( double eps )
     //add_relative_gradient_L2_norm( 1 );
     //   inner_b.add_relative_successive_improvement( 0.1 );
     //   inner_b.add_absolute_quality_improvement(100);
     //   inner_b.add_relative_quality_improvement(100);
     //   inner_b.add_absolute_successive_improvement(0.1);
     //add_absolute_vertex_movement( double eps )
     //add_absolute_vertex_movement_edge_length( double beta )
     //add_relative_vertex_movement( double eps )
     //add_absolute_successive_improvement( double eps )
     //add_relative_successive_improvement( double eps )
     ////////////////////////////optimization criteria//////////////////////////////
     //inner_b.add_absolute_gradient_L2_norm( 5e-3);
     */
    if (maxTime > 0.0) {
      double remaining = (maxTime - totalTimer.since_birth());
      if (remaining <= 0.0) {
	//MSQ_DBGOUT(2) << "Optimization is terminating without perfoming shape improvement." << std::endl;
	remaining = 0.0;
      }
      inner_b.add_cpu_time(remaining);
    }
    improver_b.set_inner_termination_criterion(&inner_b);
    InstructionQueue q_b;
    // NOTE: Don't add init_mesh to this queue because we to use
    //       the initial initial mesh (not the output of the non-barrier
    //       phase).
    q_b.set_master_quality_improver(&improver_b, err);
    MSQ_ERRRTN(err);
    q_b.add_quality_assessor(qa, err);
    MSQ_ERRRTN(err);
    //MsqMeshEntity element_Entity=pd.element_by_index(element);
    //const std::size_t *Vindex=element_Entity.get_vertex_index_array();
    //Mesh *msqMesh=pd.get_mesh();
    std::vector<Mesquite::Mesh::ElementHandle> elements;
    mesh->get_all_elements(elements, err);
    // mesh->get_all_vertices(elements, err);
    int num_elements = elements.size();
    TagHandle lid_tag1 = mesh->tag_get("LOCAL_ID1", err);
    TagHandle lid_tag2 = mesh->tag_get("LOCAL_ID2", err);
    TagHandle lid_tag3 = mesh->tag_get("LOCAL_ID3", err);
    TagHandle lid_tag4 = mesh->tag_get("LOCAL_ID4", err);
    std::vector<double> lidtest1(num_elements);
    std::vector<double> lidtest2(num_elements);
    std::vector<double> lidtest3(num_elements);
    std::vector<double> lidtest4(num_elements);
    mesh->tag_get_vertex_data(lid_tag1, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest1), err);
    mesh->tag_get_vertex_data(lid_tag2, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest2), err);
    mesh->tag_get_vertex_data(lid_tag3, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest3), err);
    mesh->tag_get_vertex_data(lid_tag4, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest4), err);
    q_b.run_common(mesh, pmesh, domain, settings, err);
    MSQ_ERRRTN(err);

    TargetWriter writer(&target);
    Settings linearMaps;
    writer.loop_over_mesh(mesh, domain, &linearMaps, err);
  }

  /*
   * 3D Mesquite Smoothing Version
   */

  AnisoWrapper3D::AnisoWrapper3D() :
      maxTime(100), mBeta(DEFAULT_BETA_3D), parallelIterations(
	  DEFUALT_PARALLEL_ITERATIONS)
  {
  }

  AnisoWrapper3D::~AnisoWrapper3D()
  {
  }

  void
  AnisoWrapper3D::set_vertex_movement_limit_factor(double beta)
  {
    assert(beta > 0.0);
    assert(beta < 1.0);
    mBeta = beta;
  }

  void
  AnisoWrapper3D::set_cpu_time_limit(double limit)
  {
    assert(limit >= 0.0);
    maxTime = limit;
  }

  void
  AnisoWrapper3D::set_parallel_iterations(int count)
  {
    if (count < 1) {
      assert(false);
      count = 1;
    }
    parallelIterations = count;
  }

  void
  AnisoWrapper3D::run_wrapper(Mesh* mesh, ParallelMesh* pmesh,
			      MeshDomain* domain, Settings* settings,
			      QualityAssessor* qa, MsqError& err)
  {
    // Use initial vertex coordinates as reference mesh
    // NOTE: need to put init_mesh in instruction queue
    //       to get it initialized with initial coords.
//  TagVertexMesh init_mesh( err, mesh ); MSQ_ERRRTN(err);
    ReferenceMesh ref_mesh(mesh);

    // Quality Metrics
    RefAspect target(&ref_mesh);

    // No-barrier phase
    TShapeNB1 mu_no;
    TQualityMetric metric_no(&target, &mu_no);

    // Barrier phase
    TShapeB1 mu_b;
    TShapeOrientB1 mu_b2;
    TQualityMetric metric_b(&target, &mu_b2);

// if decide to use an ideal isotorpic target
//    IdealShapeTarget Itarget;
//    TQualityMetric metric_bIdeal( &Itarget, &mu_b );

    // QualityAssessor
    qa->add_quality_assessment(&metric_no);
    qa->add_quality_assessment(&metric_b);

    // Phase 1: use non-barrier metric until mesh does not
    //          have any inverted elements
    //PMeanPTemplate obj_func_no( 1.0, &metric_no );
    LPtoPTemplate obj_func_no(1.0, &metric_no);
    //ConjugateGradient improver_no( &obj_func_no );
    TrustRegion improver_no(&obj_func_no);
    improver_no.use_global_patch();
    //improver_no.use_element_on_vertex_patch();
    TerminationCriterion inner_no;
    //inner_no.add_relative_gradient_L2_norm( 0.01 );
    inner_no.add_absolute_vertex_movement_edge_length(mBeta);
    inner_no.add_untangled_mesh();
    if (maxTime > 0.0)
      inner_no.add_cpu_time(maxTime);
    improver_no.set_inner_termination_criterion(&inner_no);
    InstructionQueue q_no;
    Timer totalTimer;
    //q_no.add_tag_vertex_mesh( &init_mesh, err ); MSQ_ERRRTN(err);
    q_no.add_quality_assessor(qa, err);
    MSQ_ERRRTN(err);
    q_no.set_master_quality_improver(&improver_no, err);
    MSQ_ERRRTN(err);
    q_no.add_quality_assessor(qa, err);
    MSQ_ERRRTN(err);
    q_no.run_common(mesh, pmesh, domain, settings, err);
    MSQ_ERRRTN(err);

    // Phase 2: use-barrer metric on untangled mesh
    //PMeanPTemplate obj_func_b( 1.0, &metric_b );
    LPtoPTemplate obj_func_b(1.0, &metric_b);
    //ConjugateGradient improver_b( &obj_func_b );
    TrustRegion improver_b(&obj_func_b);

    improver_b.use_element_on_vertex_patch();
//  improver_b.use_global_patch();
    TerminationCriterion inner_b;
    inner_b.add_relative_gradient_L2_norm(0.01);

    inner_b.add_absolute_vertex_movement_edge_length(mBeta);
    if (maxTime > 0.0) {
      double remaining = maxTime - totalTimer.since_birth();
      if (remaining <= 0.0) {
	//MSQ_DBGOUT(2) << "Optimization is terminating without perfoming shape improvement." << std::endl;
	remaining = 0.0;
      }
      inner_b.add_cpu_time(remaining);
    }

    improver_b.set_inner_termination_criterion(&inner_b);

    InstructionQueue q_b;
    //q_b.add_tag_vertex_mesh( &init_mesh, err ); MSQ_ERRRTN(err);
    q_b.add_quality_assessor(qa, err);
    MSQ_ERRRTN(err);
    q_b.set_master_quality_improver(&improver_b, err);
    MSQ_ERRRTN(err);
    q_b.add_quality_assessor(qa, err);
    MSQ_ERRRTN(err);
    //MsqMeshEntity element_Entity=pd.element_by_index(element);
    //const std::size_t *Vindex=element_Entity.get_vertex_index_array();
    //Mesh *msqMesh=pd.get_mesh();
    std::vector<Mesquite::Mesh::ElementHandle> elements;
    mesh->get_all_elements(elements, err);
    // mesh->get_all_vertices(elements, err);
    int num_elements = elements.size();
    TagHandle lid_tag1 = mesh->tag_get("LOCAL_ID1", err);
    TagHandle lid_tag2 = mesh->tag_get("LOCAL_ID2", err);
    TagHandle lid_tag3 = mesh->tag_get("LOCAL_ID3", err);
    TagHandle lid_tag4 = mesh->tag_get("LOCAL_ID4", err);
    TagHandle lid_tag5 = mesh->tag_get("LOCAL_ID5", err);
    TagHandle lid_tag6 = mesh->tag_get("LOCAL_ID6", err);
    TagHandle lid_tag7 = mesh->tag_get("LOCAL_ID7", err);
    TagHandle lid_tag8 = mesh->tag_get("LOCAL_ID8", err);
    TagHandle lid_tag9 = mesh->tag_get("LOCAL_ID9", err);
    std::vector<double> lidtest1(num_elements);
    std::vector<double> lidtest2(num_elements);
    std::vector<double> lidtest3(num_elements);
    std::vector<double> lidtest4(num_elements);
    std::vector<double> lidtest5(num_elements);
    std::vector<double> lidtest6(num_elements);
    std::vector<double> lidtest7(num_elements);
    std::vector<double> lidtest8(num_elements);
    std::vector<double> lidtest9(num_elements);
    mesh->tag_get_vertex_data(lid_tag1, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest1), err);
    mesh->tag_get_vertex_data(lid_tag2, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest2), err);
    mesh->tag_get_vertex_data(lid_tag3, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest3), err);
    mesh->tag_get_vertex_data(lid_tag4, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest4), err);
    mesh->tag_get_vertex_data(lid_tag5, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest5), err);
    mesh->tag_get_vertex_data(lid_tag6, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest6), err);
    mesh->tag_get_vertex_data(lid_tag7, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest7), err);
    mesh->tag_get_vertex_data(lid_tag8, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest8), err);
    mesh->tag_get_vertex_data(lid_tag9, num_elements,
			      Mesquite::arrptr(elements),
			      Mesquite::arrptr(lidtest9), err);
    q_b.run_common(mesh, pmesh, domain, settings, err);
    MSQ_ERRRTN(err);
  }

} // namespace MESQUITE_NS
#endif
