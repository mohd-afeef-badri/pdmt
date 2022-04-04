
MEDCOUPLING_ROOT_DIR=/home/mb258512/Salome/Salome-Test-Bench/SALOME-master-UB18.04/INSTALL/MEDCOUPLING
MEDFILE_ROOT_DIR=/home/mb258512/Salome/Salome-Test-Bench/SALOME-master-UB18.04/INSTALL/medfile
HDF5_ROOT_DIR=/home/mb258512/Salome/Salome-Test-Bench/SALOME-master-UB18.04/INSTALL/hdf5

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

check-medCoupling: medCoupling.cpp
	export LD_LIBRARY_PATH=$(MEDFILE_ROOT_DIR)/lib:$(HDF5_ROOT_DIR)/lib:$(MEDCOUPLING_ROOT_DIR)/lib
	$(CXX) medCoupling.cpp -I$(MEDCOUPLING_ROOT_DIR)/include -I$(MEDFILE_ROOT_DIR)/include -I$(HDF5_ROOT_DIR)/include -L$(MEDCOUPLING_ROOT_DIR)/lib -lmedloader -lmedcoupling -linterpkernel -lmedcouplingremapper $(PDMT_DEBUG) -DMEDCOUPLING
	@echo ""
