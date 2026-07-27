
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

EXTRA_DIST = check_2d_dual_vtu.py

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

check-local: test-vtk test-msh test-mesh test-med test-2d-modes test-3d test-conserve-edge test-3s test-med-3d
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
	./../run/PDMT --dimension 2 --mesh square --square_mesh_size 4 --mode subdivided_dual --out_mesh test_2d_subdivided_square.vtu > test_2d_subdivided_square.log
	python3 ./check_2d_dual_vtu.py test_2d_subdivided_square_0000.vtu --expect-polygons 15 --expect-boundary 24 --more-connectivity-than test_2d_smooth_0000.vtu
	./../run/PDMT --dimension 2 --mesh circle --circle_mesh_size 24 --mode subdivided_dual --out_mesh test_2d_subdivided_circle.vtu > test_2d_subdivided_circle.log
	python3 ./check_2d_dual_vtu.py test_2d_subdivided_circle_0000.vtu

test-3d: UnitTest3D.edp single_tetra.mesh3d
	$(MED_TEST_ENV) $(FFMPI) -np 1 UnitTest3D.edp -v 0 -nw > UnitTest3D.log

test-conserve-edge: UnitTestConserveEdge.edp physical_edge.msh3d
	$(MED_TEST_ENV) $(FFMPI) -np 1 UnitTestConserveEdge.edp -v 0 -nw > UnitTestConserveEdge.log

test-3s: UnitTest3S.edp folded_surface_edge.msh3d closed_surface.msh3d
	$(MED_TEST_ENV) $(FFMPI) -np 1 UnitTest3S.edp -v 0 -nw > UnitTest3S.log

test-med-3d: UnitTestMED3D.edp single_tetra.mesh3d closed_surface.msh3d ../msh/cylinder_edge.med ../msh/cylinder_surface_edge.med
if WITHMEDIO
	$(MED_TEST_ENV) $(FFMPI) -np 1 UnitTestMED3D.edp -v 0 -nw > UnitTestMED3D.log
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
