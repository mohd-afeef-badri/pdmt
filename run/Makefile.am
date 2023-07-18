
all-local:
	@echo ""
	@echo "*============================================================*"
	@echo " Compilation done in :"; pwd
	@echo "*============================================================*"
	@echo ""


install-exec-local:
	@mkdir -p ${prefix}
	cp meshIO.edp meshIO.edp.tmp && find . -name meshIO.edp.tmp -exec sed -i -e "s#$./../../medio/src/medio#/lib/medio#g" {} \;
if !HAVEMED
	find . -name meshIO.edp.tmp -exec sed -i -e "s#load \"./lib/medio\"##g" {} \;
	find . -name meshIO.edp.tmp -exec sed -i -e "s#Th2D = loadmedmesh(meshName, meshname = \"Mesh_1\");#cout<<\"No med support\"<<endl;#g" {} \;
endif
	cp Main.edp Main.edp.tmp && find . -name Main.edp.tmp -exec sed -i -e "s#./../src/pdmt#./lib/pdmt#g" {} \;
	mv meshIO.edp.tmp ${prefix}/meshIO.edp
	mv Main.edp.tmp ${prefix}/Main.edp

clean-local:
	@echo ""
	@echo "*============================================================*"
	@echo " Cleaning files in :"; pwd
	@echo "*============================================================*"
	@echo ""
	rm  -f *~ *.tmp


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
	@echo " performing check for pdmt"
	@echo "*============================================================*"
	@echo ""
