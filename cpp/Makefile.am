LIST_COMPILE= pdmt.cpp

HPP_DEP = pdmtFunctions.hpp \
convexHull.hpp \
polyMeshWrite.hpp \
vtkWriter.hpp  \
medWriter.hpp


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

pdmt.$(SO_SUFFIX): pdmt.cpp  $(HPP_DEP)
if HAVEMED
	$(FFCXX) -auto -mpi pdmt.cpp $(CXXFLAGS_MED) -I$(HDF5_ROOT_DIR)/include $(LDFLAGS_MED) $(LIBS_MED) $(PDMT_DEBUG) -DMEDCOUPLING
else !HAVEMED
	$(FFCXX) -auto pdmt.cpp $(PDMT_DEBUG)
	@echo ""
endif

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
