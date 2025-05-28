#--------------------------------------------------------------------------------
# Version info
#--------------------------------------------------------------------------------
MEDCOUPLING_VERSION=9.14.0
MEDIO_VERSION=1.0

#--------------------------------------------------------------------------------
# Download link
#--------------------------------------------------------------------------------
MEDCOUPLING_DOWNLOAD =	wget -c \
  			--header="User-Agent: Mozilla/5.0 (X11; Linux x86_64)" \
  			--header="Referer: https://www.salome-platform.org/" \
  			https://files.salome-platform.org/Salome/Salome$(MEDCOUPLING_VERSION)/MEDCOUPLING-$(MEDCOUPLING_VERSION).tar.gz \
  			-O MEDCOUPLING-$(MEDCOUPLING_VERSION).tar.gz

MEDIO_DOWNLOAD	=	wget -c https://github.com/mohd-afeef-badri/medio/archive/refs/tags/$(MEDIO_VERSION).tar.gz \
			-O medio-$(MEDIO_VERSION).tar.gz

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
medcoupling-build \
medio-build
else
LIST_MAKE_EXECUTE =
endif


if BUILDDEPENDENCIES
LIST_MAKE_INSTALL = medio-install
else
LIST_MAKE_INSTALL =
endif

all-local: $(LIST_MAKE_EXECUTE)
	@echo ""

install-exec-local: $(LIST_MAKE_INSTALL)
	@echo ""
	@echo "*============================================================*"
	@echo "   external depndencies intallation phase "
	@echo "*============================================================*"
	@echo ""
	

download-dependencies:
	@if [ ! -f MEDCOUPLING-$(MEDCOUPLING_VERSION).tar.gz ]; then \
		$(MEDCOUPLING_DOWNLOAD); \
	fi
	@if [ ! -f medio-$(MEDIO_VERSION).tar.gz ]; then \
		$(MEDIO_DOWNLOAD); \
	fi

#--------------------------------------------------------------------------------
# MEDCOUPLING installation
#  - Download from wget version MEDCOUPLING_VERSION
#  - User can override download by providing ext/MEDCOUPLING-$(MEDCOUPLING_VERSION).tar.gz
#--------------------------------------------------------------------------------

medcoupling-build: MEDCOUPLING-$(MEDCOUPLING_VERSION)/INSTALL/MEDCOUPLING/lib/libmedcoupling.so medcoupling-install

MEDCOUPLING-$(MEDCOUPLING_VERSION)/INSTALL/MEDCOUPLING/lib/libmedcoupling.so: MEDCOUPLING-$(MEDCOUPLING_VERSION).tar.gz
	@echo ""
	@echo "*============================================================*"
	@echo " MEDCOUPLING installation began for you be patient:   "
	@echo "   - Download using wget version $(MEDCOUPLING_VERSION)"
	@echo "   - to overide provide ext/MEDCOUPLING-$(MEDCOUPLING_VERSION).tar.gz"
	@echo "*============================================================*"
	@echo ""
	@echo ""
	@echo "*---------------------------------*"
	@echo "   untaring MEDCOUPLING-$(MEDCOUPLING_VERSION).tar.gz"
	@echo "*---------------------------------*"
	@echo ""
	tar -xzf MEDCOUPLING-$(MEDCOUPLING_VERSION).tar.gz
	@echo ""
	@echo "*---------------------------------*"
	@echo "   creating BUILD and INSTALL directories"
	@echo "*---------------------------------*"
	@echo ""
	cd MEDCOUPLING-$(MEDCOUPLING_VERSION)   && \
	mkdir -p BUILD   && \
	mkdir -p INSTALL
	@echo ""
	@echo "*---------------------------------*"
	@echo "   installing medcoupling"
	@echo "*---------------------------------*"
	@echo ""
	cd MEDCOUPLING-$(MEDCOUPLING_VERSION)/SOURCES       && \
	tar -xzf ./../ARCHIVES/MEDCOUPLING.tgz && \
	cd ./../BUILD && \
	source ./../env_launch.sh && \
	cmake  $(CMAKE_FLAGS_MEDCOUPLING) ./../SOURCES/MEDCOUPLING && \
	make -j16

medcoupling-install: MEDCOUPLING-$(MEDCOUPLING_VERSION)/INSTALL/MEDCOUPLING/lib/libmedcoupling.so
	@echo ""
	@echo "*---------------------------------*"
	@echo "   installing medcoupling"
	@echo "*---------------------------------*"
	@echo ""
	cd MEDCOUPLING-$(MEDCOUPLING_VERSION)/BUILD && source ./../env_launch.sh && make install -j all


#--------------------------------------------------------------------------------
# MEDIO installation
#  - Download from wget version MEDIO_VERSION
#  - User can override download by providing ext/medio-$(MEDIO_VERSION).tar.gz
#--------------------------------------------------------------------------------

medio-build: medio-$(MEDIO_VERSION)/src/medio.so

medio-$(MEDIO_VERSION)/src/medio.so: medio-$(MEDIO_VERSION).tar.gz
	@echo ""
	@echo "*============================================================*"
	@echo " MEDIO installation began for you be patient:   "
	@echo "   - Download using wget version $(MEDIO_VERSION)"
	@echo "   - to overide provide ext/medio-$(MEDIO_VERSION).tar.gz"
	@echo "*============================================================*"
	@echo ""
	@echo ""
	@echo "*---------------------------------*"
	@echo "   untaring medio-$(MEDIO_VERSION).tar.gz"
	@echo "*---------------------------------*"
	@echo ""
	tar -xzf medio-$(MEDIO_VERSION).tar.gz
	@echo ""
	cd medio-$(MEDIO_VERSION)/src && \
	$(FFCXX) -auto medio.cpp $(CXXFLAGS_MED)  $(LDFLAGS_MED) $(LIBS_MED) -DMEDCOUPLING

medio-install:  medio-$(MEDIO_VERSION)/src/medio.so
if BUILDDEPENDENCIES
	@echo ""
	@echo "*---------------------------------*"
	@echo "   installing medio"
	@echo "*---------------------------------*"
	@echo ""
	@mkdir -p ${prefix}/lib
	@cp medio-$(MEDIO_VERSION)/src/medio.so ${prefix}/lib/$(DEP1)	
endif

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
	rm  -f medio-$(MEDIO_VERSION)/src/medio.so

maintainer-clean-local: clean-local
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning automake generated files"
	@echo "*============================================================*"
	@echo ""
	rm  -rf Makefile Makefile.in
