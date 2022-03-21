LIST_COMPILE= pdmt.cpp

MEDCOUPLING_ROOT_DIR=/home/mb258512/Salome/Salome-Test-Bench/SALOME-master-UB18.04/INSTALL/MEDCOUPLING
MEDFILE_ROOT_DIR=/home/mb258512/Salome/Salome-Test-Bench/SALOME-master-UB18.04/INSTALL/medfile
HDF5_ROOT_DIR=/home/mb258512/Salome/Salome-Test-Bench/SALOME-master-UB18.04/INSTALL/hdf5

DEP1 = pdmt.$(SO_SUFFIX)

#PDMT_DEBUG = -DDEBUG
PDMT_DEBUG = -DNDEBUG
	
all-local: $(DEP1)
	@echo ""
	@echo "*============================================================*"
	@echo " Compilation done in :"; pwd
	@echo "*============================================================*"
	@echo ""
	@echo " Compilation list $(LIST_COMPILE)"

pdmt.$(SO_SUFFIX): pdmt.cpp pdmtFunctions.hpp convexHull.hpp polyMeshWrite.hpp
	$(FFCXX) -auto pdmt.cpp $(PDMT_DEBUG)
	@echo ""

medCoupling.$(SO_SUFFIX): medCoupling.cpp
	export LD_LIBRARY_PATH=$(MEDFILE_ROOT_DIR)/lib:$(HDF5_ROOT_DIR)/lib:$(MEDCOUPLING_ROOT_DIR)/lib
	$(FFCXX) -auto -mpi medCoupling.cpp -I$(MEDCOUPLING_ROOT_DIR)/include -I$(MEDFILE_ROOT_DIR)/include -I$(HDF5_ROOT_DIR)/include -L$(MEDCOUPLING_ROOT_DIR)/lib -lmedloader -lmedcoupling -linterpkernel -lmedcouplingremapper
	@echo ""

install-exec-local:
	@mkdir -p ${prefix}/include
	@mkdir -p ${prefix}/include/psd
	@cp $(DEP1) ${prefix}/include/psd/$(DEP1)
	
clean-local:
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning files in :"; pwd
	@echo "*============================================================*"
	@echo ""
	rm  -f *~  *.$(SO_SUFFIX) *.o


maintainer-clean-local: clean-local
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning automake generated files"
	@echo "*============================================================*"
	@echo ""	
	rm  -rf Makefile Makefile.in	
