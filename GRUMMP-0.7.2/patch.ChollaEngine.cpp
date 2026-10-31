--- /usr/anslab/src/cgma-13.1/geom/Cholla/ChollaEngine.cpp	2014-09-08 22:47:34.000000000 -0700
+++ cgma-13.1/geom/Cholla/ChollaEngine.cpp	2013-12-10 22:28:12.000000000 -0800
@@ -1168,61 +1168,6 @@
   for ( kk = cholla_surface_list.size(); kk > 0; kk-- )
   {
     ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
-    DLIList<FacetEntity*> facet_entity_list;
-    DLIList<CubitPoint*> point_list;
-    chsurf_ptr->get_points(point_list);
-    chsurf_ptr->get_facets(facet_entity_list);
-    DLIList<CubitFacet*> facet_list;
-    CAST_LIST( facet_entity_list, facet_list, CubitFacet );
-
-    FacetEvalTool *eval_tool_ptr = new FacetEvalTool( facet_list,point_list,
-                                                      interp_order,min_dot );
-    chsurf_ptr->assign_eval_tool(eval_tool_ptr);
-
-    // go through each of this surface's curves and create CurveFacetEvalTools
-
-    DLIList<ChollaCurve *> chcurv_list;
-    chsurf_ptr->get_curves( chcurv_list );
-    for (ii=0; ii<chcurv_list.size(); ii++)
-    {
-      ChollaCurve *chcurv_ptr = chcurv_list.get_and_step();
-      if (chcurv_ptr->get_eval_tool() == NULL)
-      {
-        CubitSense orientation_wrt_surface;
-        determine_curve_orientation( chsurf_ptr, chcurv_ptr, orientation_wrt_surface );
-        CubitPoint *start_ptr, *end_ptr;
-        chcurv_ptr->get_ends( start_ptr, end_ptr );
-
-        CurveFacetEvalTool *curv_eval_tool_ptr = new CurveFacetEvalTool;
-        stat = curv_eval_tool_ptr->initialize( eval_tool_ptr, start_ptr, end_ptr, orientation_wrt_surface );
-        if ( stat != CUBIT_SUCCESS )
-        {
-            return stat;
-        }
-        chcurv_ptr->assign_eval_tool( curv_eval_tool_ptr );
-      }
-    }
-  }
-  return stat;
-}
-/*{
-  CubitStatus stat = CUBIT_SUCCESS;
-  int ii, kk;
-
-  // make sure the facet flags have been reset
-
-  for ( kk = cholla_surface_list.size(); kk > 0; kk-- )
-  {
-    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
-    chsurf_ptr->reset_facet_flags();
-  }
-
-  // now loop through surfaces and create them
-
-  int mydebug = 0;
-  for ( kk = cholla_surface_list.size(); kk > 0; kk-- )
-  {
-    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
     DLIList<FacetEntity*> facet_entity_list;    
     chsurf_ptr->get_facets(facet_entity_list);
     DLIList<CubitFacet*> facet_list;
@@ -1282,7 +1227,7 @@
     }
   }
   return stat;
-}*/
+}
 
 CubitStatus ChollaEngine::rebuild_surface_and_curve_eval_tools(
   DLIList<ChollaSurface*> &cholla_surface_list,
