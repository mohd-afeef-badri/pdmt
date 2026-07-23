
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

check-local: test-vtk test-msh test-mesh test-med test-3d test-conserve-edge test-3s test-med-3d
	@echo ""
	@echo "*============================================================*"
	@echo " performing check for pdmt"
	@echo "*============================================================*"
	@echo ""

UnitTestCovexHull: UnitTestCovexHull.edp
	$(FFMPI) -np 1 UnitTestCovexHull.edp -v 0 -ns -nw > UnitTestCovexHull.log

test-3d: UnitTest3D.edp single_tetra.mesh3d
	$(FFMPI) -np 1 UnitTest3D.edp -v 0 -nw > UnitTest3D.log

test-conserve-edge: UnitTestConserveEdge.edp physical_edge.msh3d
	$(FFMPI) -np 1 UnitTestConserveEdge.edp -v 0 -nw > UnitTestConserveEdge.log

test-3s: UnitTest3S.edp folded_surface_edge.msh3d closed_surface.msh3d
	$(FFMPI) -np 1 UnitTest3S.edp -v 0 -nw > UnitTest3S.log

test-med-3d: UnitTestMED3D.edp single_tetra.mesh3d closed_surface.msh3d
if WITHMEDIO
	$(FFMPI) -np 1 UnitTestMED3D.edp -v 0 -nw > UnitTestMED3D.log
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
