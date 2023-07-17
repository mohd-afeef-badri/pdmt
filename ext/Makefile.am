#--------------------------------------------------------------------------------
# Version info
#--------------------------------------------------------------------------------
MEDCOUPLING_VERSION=9.11.0
SALOME_VERSION=9.11.0


#--------------------------------------------------------------------------------
# Download link #TODO
#--------------------------------------------------------------------------------
MEDCOUPLING_DOWNLOAD      =  wget -c https://files.salome-platform.org/Salome/Salome$(MEDCOUPLING_VERSION)/MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI.tar.gz \
                       -O MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI.tar.gz

#--------------------------------------------------------------------------------
# CMAKE FLAGS FOR MEDCOUPLING
#--------------------------------------------------------------------------------
CMAKE_FLAGS_MEDCOUPLING = \
-DCMAKE_BUILD_TYPE=Release \
-DMEDCOUPLING_ENABLE_PYTHON=OFF \
-DMEDCOUPLING_ENABLE_PARTITIONER=OFF \
-DMEDCOUPLING_ENABLE_RENUMBER=OFF \
-DMEDCOUPLING_PARTITIONER_METIS=OFF \
-DMEDCOUPLING_PARTITIONER_PARMETIS=OFF \
-DMEDCOUPLING_PARTITIONER_SCOTCH=OFF \
-DMEDCOUPLING_PARTITIONER_PTSCOTCH=OFF \
-DMEDCOUPLING_MICROMED=OFF \
-DMEDCOUPLING_BUILD_DOC=OFF \
-DCMAKE_EXE_LINKER_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=1" \
-DMEDCOUPLING_BUILD_TESTS=OFF \
-DSALOME_USE_MPI=ON \
-DCMAKE_INSTALL_PREFIX=./../INSTALL/MEDCOUPLING


#--------------------------------------------------------------------------------
# TRIGGER BUILD OF MEDCOULING BASED ON --build-dependencies flag
#--------------------------------------------------------------------------------
if BUILDDEPENDENCIES
LIST_MAKE_EXECUTE = download-dependencies \
medcoupling-build
else
LIST_MAKE_EXECUTE =
endif


all-local: $(LIST_MAKE_EXECUTE)
	@echo ""



download-dependencies:
	@if [ ! -f MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI.tar.gz ]; then \
		$(MEDCOUPLING_DOWNLOAD); \
	fi

#--------------------------------------------------------------------------------
# MEDCOUPLING installation
#  - Download from wget version MEDCOUPLING_VERSION
#  - User can override download by providing ext/mfront-$(MFRONT_VERSION).tar.gz
#--------------------------------------------------------------------------------

medcoupling-build: MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI/INSTALL/MEDCOUPLING/lib/libmedcoupling.so

MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI/INSTALL/MEDCOUPLING/lib/libmedcoupling.so: MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI.tar.gz
	@echo ""
	@echo "*============================================================*"
	@echo " MEDCOUPLING installation began for you be patient:   "
	@echo "   - Download using wget version $(MEDCOUPLING_VERSION)"
	@echo "   - to overide provide ext/MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI.tar.gz"
	@echo "*============================================================*"
	@echo ""
	@echo ""
	@echo "*---------------------------------*"
	@echo "   untaring MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI.tar.gz"
	@echo "*---------------------------------*"
	@echo ""
	tar -xzf MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI.tar.gz
	@echo ""
	@echo "*---------------------------------*"
	@echo "   creating BUILD and INSTALL directories"
	@echo "*---------------------------------*"
	@echo ""
	cd MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI   && \
	mkdir -p BUILD   && \
	mkdir -p INSTALL
	@echo ""
	@echo "*---------------------------------*"
	@echo "   installing medcoupling"
	@echo "*---------------------------------*"
	@echo ""
	cd MEDCOUPLING-$(MEDCOUPLING_VERSION)-MPI/SOURCES       && \
	tar -xzf ./../ARCHIVES/MEDCOUPLING.tgz && \
	cd ./../BUILD && \
	source ./../env_launch.sh && \
	cmake  $(CMAKE_FLAGS_MEDCOUPLING) ./../SOURCES/MEDCOUPLING && \
	make -j16 && \
	make install


check-local:
	@echo ""
	@echo "*============================================================*"
	@echo " Nothing to check here:"
	@echo "*============================================================*"
	@echo ""

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

install-exec-local:
	@echo ""
	@echo "*============================================================*"
	@echo " Nothing to install here:"
	@echo "*============================================================*"
	@echo ""
