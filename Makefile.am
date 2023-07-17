SUBDIRS = ext src run test


all-local:
	@echo ""
	@echo "*============================================================*"
	@echo " pdmt version 0.1 compilation successfully done"
	@echo ""
	@echo " In order to use medio please ensure to edit LD_LIBRARY_PATH"
	@echo " use the following"
	@echo ""
if HAVEMED
	@echo ' export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(hdf5)/lib:$(medfile)/lib:$(medcoupling)/lib'
else !HAVEMED
	@echo ' export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(PDMT_HDF5_DIR)/lib:$(PDMT_MEDFILE_DIR)/lib:$(PDMT_MEDCOUPLING_DIR)/lib'
endif
	@echo ""
	@echo " after <make install> pdmt can be used with FreeFEM"
	@echo " you will find it in ${prefix} folder."
	@echo ""
	@echo "*============================================================*"
	@echo ""

clean-local:
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning files in :"; pwd
	@echo "*============================================================*"
	@echo ""
	rm  -f *~ *.vtk *.vtu *.med *.txt

install-exec-local:
	@mkdir -p ${prefix}
	@cp $(PDMT_MEDIO_DIR).so ${prefix}/lib/.

maintainer-clean-local: clean-local
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning automake generated files"
	@echo "*============================================================*"
	@echo ""	
	rm  -rf Makefile Makefile.in
