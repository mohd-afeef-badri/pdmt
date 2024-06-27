LIST_COMPILE= pdmt.cpp

HPP_DEP = pdmtFunctions.hpp \
convexHull.hpp \
polyMeshWrite.hpp \
vtkWriter.hpp  \
typWriter.hpp  \
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
if WITHMEDIO
	$(FFCXX) -auto -mpi pdmt.cpp $(CXXFLAGS_MED) $(LDFLAGS_MED) $(LIBS_MED) $(PDMT_DEBUG) -DMEDCOUPLING
endif
if !WITHMEDIO
	$(FFCXX) -auto -mpi pdmt.cpp $(PDMT_DEBUG)
	@echo ""
	@echo " PDMT will not support med meshes"
endif

install-exec-local:
	@mkdir -p ${prefix}/lib
	@cp $(DEP1) ${prefix}/lib/$(DEP1)

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
