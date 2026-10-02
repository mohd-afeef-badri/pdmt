
clean-local:
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning files in :"; pwd
	@echo "*============================================================*"
	@echo ""
	rm  -f *~ *.log
	rm  -f *.tmp *.vtk *.mesh *.vtu *.med


maintainer-clean-local: clean-local
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning automake generated files"
	@echo "*============================================================*"
	@echo ""
	rm  -rf Makefile Makefile.in

EXTRA_DIST = check_2d_dual_vtu.py check_2d_viscous_vtu.py \
	check_3d_regularization.py check_3s_regularization.py \
	CheckViscousMedGroup.edp MakeViscousMed2D.edp \
	open_surface.msh3d viscous_square.msh

# These tests load build-tree plugins directly and bypass run/PDMT, which
# normally exports the MED runtime library paths. Use configure substitutions
# here so local builds and CI workspaces do not need a hard-coded path.
MED_TEST_ENV =
if HAVEDEPENDENCIES
MED_TEST_ENV += LD_LIBRARY_PATH="$(PDMT_HDF5_DIR)/lib:$(PDMT_MEDFILE_DIR)/lib:$(PDMT_MEDCOUPLING_DIR)/lib:$${LD_LIBRARY_PATH}"
endif
if HAVEMED
MED_TEST_ENV += LD_LIBRARY_PATH="$(hdf5)/lib:$(medfile)/lib:$(medcoupling)/lib:$${LD_LIBRARY_PATH}"
endif

check-local: test-vtk test-msh test-mesh test-med test-2d-modes test-2d-viscous test-med-2d-viscous test-3d test-conserve-edge test-3s test-med-3d
	@echo ""
	@echo "*============================================================*"
	@echo " performing check for pdmt"
	@echo "*============================================================*"
	@echo ""

UnitTestCovexHull: UnitTestCovexHull.edp
	$(MED_TEST_ENV) $(FFMPI) -np 1 UnitTestCovexHull.edp -v 0 -ns -nw > UnitTestCovexHull.log

test-2d-modes: check_2d_dual_vtu.py
	./../run/PDMT --dimension 2 --mesh square --square_mesh_size 4 --mode smooth_dual --out_mesh test_2d_smooth.vtu > test_2d_smooth.log
	python3 ./check_2d_dual_vtu.py test_2d_smooth_0000.vtu --expect-polygons 15 --expect-boundary 24
	./../run/PDMT --dimension 2 --mesh square --square_mesh_size 4 --mode smooth_dual --smooth_iterations 3 --smooth_relaxation 0.3 --out_mesh test_2d_smooth_regularized.vtu > test_2d_smooth_regularized.log
	python3 ./check_2d_dual_vtu.py test_2d_smooth_regularized_0000.vtu --expect-polygons 15 --expect-boundary 24 --regularized-from test_2d_smooth_0000.vtu
	./../run/PDMT --dimension 2 --mesh square --square_mesh_size 4 --mode subdivided_dual --out_mesh test_2d_subdivided_square.vtu > test_2d_subdivided_square.log
	python3 ./check_2d_dual_vtu.py test_2d_subdivided_square_0000.vtu --expect-polygons 15 --expect-boundary 24 --more-connectivity-than test_2d_smooth_0000.vtu
	./../run/PDMT --dimension 2 --mesh square --square_mesh_size 4 --mode subdivided_dual --smooth_iterations 3 --smooth_relaxation 0.3 --out_mesh test_2d_subdivided_regularized.vtu > test_2d_subdivided_regularized.log
	python3 ./check_2d_dual_vtu.py test_2d_subdivided_regularized_0000.vtu --expect-polygons 15 --expect-boundary 24 --regularized-from test_2d_subdivided_square_0000.vtu
	./../run/PDMT --dimension 2 --mesh circle --circle_mesh_size 24 --mode subdivided_dual --out_mesh test_2d_subdivided_circle.vtu > test_2d_subdivided_circle.log
	python3 ./check_2d_dual_vtu.py test_2d_subdivided_circle_0000.vtu

test-2d-viscous: check_2d_dual_vtu.py check_2d_viscous_vtu.py viscous_square.msh
	./../run/PDMT --dimension 2 --mesh ./viscous_square.msh --mode smooth_dual --viscous-layer yes --viscous-layer_groups wall_top --viscous-layer_count 10 --viscous-layer_thickness 0.6 --out_mesh test_2d_viscous.vtu > test_2d_viscous.log
	python3 ./check_2d_viscous_vtu.py test_2d_viscous_0000.vtu
	./../run/PDMT --dimension 2 --mesh ./viscous_square.msh --mode subdivided_dual --viscous-layer yes --viscous-layer_groups wall_top --viscous-layer_count 10 --viscous-layer_thickness 0.6 --out_mesh test_2d_viscous_subdivided.vtu > test_2d_viscous_subdivided.log
	python3 ./check_2d_dual_vtu.py test_2d_viscous_subdivided_0000.vtu --expect-polygons 25 --expect-boundary 28 --more-connectivity-than test_2d_viscous_0000.vtu
	@if ./../run/PDMT --dimension 2 --mesh ./viscous_square.msh --viscous-layer yes --viscous-layer_groups missing --viscous-layer_count 2 --viscous-layer_thickness 0.1 --out_mesh invalid.vtu > test_2d_viscous_invalid.log 2>&1; then echo "unknown viscous-layer group unexpectedly succeeded"; exit 1; fi
	grep -q "unknown or empty physical boundary group(s): missing" test_2d_viscous_invalid.log

test-3d: UnitTest3D.edp single_tetra.mesh3d offcenter_cube.mesh3d
	$(MED_TEST_ENV) $(FFMPI) -np 1 UnitTest3D.edp -v 0 -nw > UnitTest3D.log
	python3 ./check_3d_regularization.py offcenter_unsmoothed_0000.vtu offcenter_regularized_0000.vtu

test-conserve-edge: UnitTestConserveEdge.edp physical_edge.msh3d
	$(MED_TEST_ENV) $(FFMPI) -np 1 UnitTestConserveEdge.edp -v 0 -nw > UnitTestConserveEdge.log

test-3s: UnitTest3S.edp folded_surface_edge.msh3d closed_surface.msh3d open_surface.msh3d
	$(MED_TEST_ENV) $(FFMPI) -np 1 UnitTest3S.edp -v 0 -nw > UnitTest3S.log
	python3 ./check_3s_regularization.py open_surface_smooth_0000.vtu open_surface_smooth_regularized_0000.vtu
	python3 ./check_3s_regularization.py open_surface_subdivided_0000.vtu open_surface_subdivided_regularized_0000.vtu
	python3 ./check_3s_regularization.py open_surface_feature_0000.vtu open_surface_feature_regularized_0000.vtu

if WITHMEDIO
test-med-2d-viscous: CheckViscousMedGroup.edp MakeViscousMed2D.edp viscous_square.msh check_2d_dual_vtu.py check_2d_viscous_vtu.py
	$(MED_TEST_ENV) $(FFMPI) -np 1 MakeViscousMed2D.edp -v 0 -nw > MakeViscousMed2D.log
	$(MED_TEST_ENV) ./../run/PDMT --dimension 2 --mesh ./viscous_square_input.med --med_mesh_name TriangularMesh --mode smooth_dual --viscous-layer yes --viscous-layer_groups edge_group_12 --viscous-layer_count 10 --viscous-layer_thickness 0.6 --out_mesh test_med_2d_viscous.vtu > test_med_2d_viscous.log
	python3 ./check_2d_viscous_vtu.py test_med_2d_viscous_0000.vtu
	$(MED_TEST_ENV) ./../run/PDMT --dimension 2 --mesh ./viscous_square_input.med --med_mesh_name TriangularMesh --mode smooth_dual --viscous-layer yes --viscous-layer_groups edge_group_12 --viscous-layer_count 10 --viscous-layer_thickness 0.6 --out_mesh test_med_2d_viscous.med > test_med_2d_viscous_med.log
	$(MED_TEST_ENV) $(FFMPI) -np 1 CheckViscousMedGroup.edp -v 0 -nw > CheckViscousMedGroup.log
	$(MED_TEST_ENV) ./../run/PDMT --dimension 2 --mesh ./viscous_square_input.med --med_mesh_name TriangularMesh --mode subdivided_dual --viscous-layer yes --viscous-layer_groups edge_group_12 --viscous-layer_count 10 --viscous-layer_thickness 0.6 --out_mesh test_med_2d_viscous_subdivided.vtu > test_med_2d_viscous_subdivided.log
	python3 ./check_2d_dual_vtu.py test_med_2d_viscous_subdivided_0000.vtu --expect-polygons 25 --expect-boundary 28 --more-connectivity-than test_med_2d_viscous_0000.vtu
	@if $(MED_TEST_ENV) ./../run/PDMT --dimension 2 --mesh ./viscous_square_input.med --med_mesh_name TriangularMesh --viscous-layer yes --viscous-layer_groups missing --viscous-layer_count 2 --viscous-layer_thickness 0.1 --out_mesh invalid_med.vtu > test_med_2d_viscous_invalid.log 2>&1; then echo "unknown MED viscous-layer group unexpectedly succeeded"; exit 1; fi
	grep -q "unknown MED boundary group(s): missing" test_med_2d_viscous_invalid.log

test-med-3d: UnitTestMED3D.edp single_tetra.mesh3d closed_surface.msh3d folded_surface_edge.msh3d
	$(MED_TEST_ENV) $(FFMPI) -np 1 UnitTestMED3D.edp -v 0 -nw > UnitTestMED3D.log
else
test-med-2d-viscous:

test-med-3d:
endif

test-typ:
	./../run/PDMT --debug --mesh ./../msh/pentagon.mesh  --out_mesh out.typ2

test-vtk:
	./../run/PDMT --debug --mesh ./../msh/disk5holes.vtk

test-msh:
	./../run/PDMT --debug --mesh ./../msh/disk1hole.msh

test-mesh:
	./../run/PDMT --debug --mesh ./../msh/pentagon.mesh

test-med:
if WITHMEDIO
	./../run/PDMT --debug --mesh ./../msh/triangle1hole.med --med_mesh_name Mesh_1
endif
