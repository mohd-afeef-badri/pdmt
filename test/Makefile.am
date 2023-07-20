
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

check-local: check-convex-hull test-vtk test-msh test-mesh test-med
	@echo ""
	@echo "*============================================================*"
	@echo " performing check for pdmt"
	@echo "*============================================================*"
	@echo ""

check-convex-hull: UnitTestCovexHull.edp
	$(FFMPI) -np 1 UnitTestCovexHull.edp -v 0 -ns -nw > UnitTestCovexHull.log

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
