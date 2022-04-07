
clean-local:
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning files in :"; pwd
	@echo "*============================================================*"
	@echo ""
	rm  -f *~


maintainer-clean-local: clean-local
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning automake generated files"
	@echo "*============================================================*"
	@echo ""
	rm  -rf Makefile Makefile.in

check-local: check-convex-hull
	@echo ""
	@echo "*============================================================*"
	@echo " performingg check for pdmt"
	@echo "*============================================================*"
	@echo ""

check-convex-hull: UnitTestCovexHull.edp
	./../../../psd-solve/PSD_Solve -np 1 UnitTestCovexHull.edp -v 0 -ns -nw > UnitTestCovexHull.log
