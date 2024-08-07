FF_IS = `which $(FFMPI)`


EXPORT_PDMT  =\nexport FF_LOADPATH="${prefix}/lib";\nexport FF_INCLUDEPATH="${prefix}/lib"
EXPORT_PDMT +=\nexport PATH="${prefix}/bin:$$PATH";\nexport LD_LIBRARY_PATH="${prefix}/lib:$$LD_LIBRARY_PATH"


if HAVEDEPENDENCIES
EXPORT_MED  =\nexport out_dir_Path=${prefix}/bin
EXPORT_MED +=\nexport PATH="$${out_dir_Path}$(PDMT_HDF5_DIR)/bin:$${PATH}"
EXPORT_MED +=\nexport LD_LIBRARY_PATH="$(PDMT_HDF5_DIR)/lib:$${LD_LIBRARY_PATH}"
EXPORT_MED +=\nexport PATH="$${out_dir_Path}$(PDMT_MEDFILE_DIR)/bin:$${PATH}"
EXPORT_MED +=\nexport LD_LIBRARY_PATH="$(PDMT_MEDFILE_DIR)/lib:$${LD_LIBRARY_PATH}"
EXPORT_MED +=\nexport MEDCOUPLING_ROOT_DIR="$(PDMT_MEDCOUPLING_DIR)"
EXPORT_MED +=\nexport PATH="$${MEDCOUPLING_ROOT_DIR}/bin:$${PATH}"
EXPORT_MED +=\nexport LD_LIBRARY_PATH="$${MEDCOUPLING_ROOT_DIR}/lib:$${LD_LIBRARY_PATH}"
endif

if HAVEMED
EXPORT_MED  =\nexport out_dir_Path=${prefix}/bin
EXPORT_MED +=\nexport PATH="$${out_dir_Path}$(hdf5)/bin:$${PATH}"
EXPORT_MED +=\nexport LD_LIBRARY_PATH="$(hdf5)/lib:$${LD_LIBRARY_PATH}"
EXPORT_MED +=\nexport PATH="$${out_dir_Path}/$(medfile)/bin:$${PATH}"
EXPORT_MED +=\nexport LD_LIBRARY_PATH="$(medfile)lib:$${LD_LIBRARY_PATH}"
EXPORT_MED +=\nexport MEDCOUPLING_ROOT_DIR="$(medcoupling)"
EXPORT_MED +=\nexport PATH="$${MEDCOUPLING_ROOT_DIR}/bin:$${PATH}"
EXPORT_MED +=\nexport LD_LIBRARY_PATH="$${MEDCOUPLING_ROOT_DIR}/lib:$${LD_LIBRARY_PATH}"
endif

if !WITHMEDIO
EXPORT_MED  =
endif

all-local: PDMT
	@echo ""
	@echo "*============================================================*"
	@echo " Compilation done in :"; pwd
	@echo "*============================================================*"
	@echo ""


install-exec-local:
	@mkdir -p ${prefix}
	@mkdir -p ${prefix}/lib
	@mkdir -p ${prefix}/bin
	cp meshIO.idp meshIO.idp.tmp
if !WITHMEDIO
	find . -name meshIO.idp.tmp -exec sed -i -e "s#load \"medio\"##g" {} \;
	find . -name meshIO.idp.tmp -exec sed -i -e "s#Th2D = loadmedmesh(meshName, meshname = medMeshName);#cout<<\"No med support PDMT not compiled with medio\"<<endl;#g" {} \;
endif
	cp pdmt.edp pdmt.edp.tmp && find . -name pdmt.edp.tmp -exec sed -i -e "s#./../src/pdmt#pdmt#g" {} \;
	mv meshIO.idp.tmp ${prefix}/lib/meshIO.idp
	mv pdmt.edp.tmp ${prefix}/lib/pdmt.edp
	cp ./PDMT ${prefix}/bin/.

PDMT:
	cat $(FF_IS) > PDMT
	find . -name PDMT -exec sed -i -e 's#mpirun=mpirun#mpirun=mpirun;\n\n$(EXPORT_PDMT)$(EXPORT_MED)\n#g' {} \;
	find . -name PDMT -exec sed -i -e 's#ffmpi=FreeFem++-mpi#\n\n<<COMMENT_BLOCK\nffmpi=FreeFem++-mpi#g' {} \;
	find . -name PDMT -exec sed -i -e 's#if \[ -d #COMMENT_BLOCK\n\nif \[ -d #g' {} \;
	find . -name PDMT -exec sed -i -e 's#a\[0\]#\n\n\# Default arguments\ndefault_args=\"-np 1 FreeFem++-mpi ${prefix}/lib/pdmt.edp\"\n\na\[0\]#g' {} \;
	find . -name PDMT -exec sed -i -e 's#---#---\n\n\# Always include default arguments first\nset -- $$default_args \"$$@\"\n\n\n#g' {} \;
#	find . -name PDMT -exec sed -i -e 's#mpi_run_option=\"--oversubscribe\"#mpi_run_option=\"--oversubscribe -np 1 FreeFem++-mpi ${prefix}/lib/pdmt.edp\"#g' {} \;
#	find . -name PDMT -exec sed -i -e 's#mpi_run_option=\"\"#mpi_run_option=\"-np 1 FreeFem++-mpi ${prefix}/lib/pdmt.edp\"#g' {} \;
#	find . -name PDMT -exec sed -i -e 's#a\[((j++))\]=\"$$mpi_run_option\"#a\[((j++))\]=\"$$mpi_run_option -np 1 FreeFem++-mpi ${prefix}/lib/pdmt.edp\"#g' {} \;
	chmod a+x ./PDMT
	less ./PDMT

a[((j++))]="$mpi_run_option"

check-local:
	@echo ""
	@echo "*============================================================*"
	@echo " performing check for pdmt"
	@echo "*============================================================*"
	@echo ""
	@echo "*------------------------------------------------------------*"
	@echo " test 1: square with debug"
	@echo "*------------------------------------------------------------*"
	@echo ""
	./PDMT --debug --mesh square
	@echo ""
	@echo "*------------------------------------------------------------*"
	@echo " test 2: square no debug"
	@echo "*------------------------------------------------------------*"
	@echo ""
	./PDMT --mesh square
	@echo ""
	@echo "*------------------------------------------------------------*"
	@echo " test 3: circle with debug"
	@echo "*------------------------------------------------------------*"
	@echo ""
	./PDMT --debug --mesh circle
	@echo ""
	@echo "*------------------------------------------------------------*"
	@echo " test 4: circle no debug"
	@echo "*------------------------------------------------------------*"
	@echo ""
	./PDMT  --mesh circle
	@echo ""
	@echo "*------------------------------------------------------------*"
	@echo " test 5: circle no debug pdmt_mesh.typ2"
	@echo "*------------------------------------------------------------*"
	@echo ""
	./PDMT  --mesh circle --out_mesh pdmt_mesh.typ2


clean-local:
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning files in :"; pwd
	@echo "*============================================================*"
	@echo ""
	rm  -f *~ *.tmp *.vtk *.mesh *.vtu *.med *.txt *.typ2 ./PDMT

maintainer-clean-local: clean-local
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning automake generated files"
	@echo "*============================================================*"
	@echo ""
	rm  -rf Makefile Makefile.in

