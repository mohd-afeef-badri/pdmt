AC_INIT([pdmt],[0.1],[mohd-afeef.badri@cea.fr])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_HEADERS([config.h])

dnl ---------------------------------------------------------------
dnl Small macro to print out a nice, pretty section title.
dnl ---------------------------------------------------------------
define(SECTION_TITLE,
[
	echo >& AS_MESSAGE_FD
	echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
	echo '   $1' >& AS_MESSAGE_FD
	echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
])


AU_ALIAS([ACX_FF], [AX_FF])
AC_DEFUN([AX_FF], [
AC_PREREQ([2.50])
AC_ARG_VAR(FF,[FreeFEM compiler command])
AC_CHECK_PROG(FF,[FreeFem++],[FreeFem++],[no])
])

AU_ALIAS([ACX_FFCXX], [AX_FFCXX])
AC_DEFUN([AX_FFCXX], [
AC_PREREQ([2.50])
AC_ARG_VAR(FFCXX,[FreeFEM compiler command])
AC_CHECK_PROG(FFCXX,[ff-c++],[ff-c++],[no])
])

AU_ALIAS([ACX_FFCXXDEP], [AX_FFCXXDEP])
AC_DEFUN([AX_FFCXXDEP], [
AC_PREREQ([2.50])
AC_ARG_VAR(FFCXXDEP,[ff-c++ compiler dependencies])
AC_CHECK_PROG(FFCXXDEP,[ff-get-dep],[ff-get-dep],[no])
])

AU_ALIAS([ACX_FFMPI], [AX_FFMPI])
AC_DEFUN([AX_FFMPI], [
AC_PREREQ([2.50])
AC_ARG_VAR(FFMPI,[FreeFEM compiler command])
AC_CHECK_PROG(FFMPI,[ff-mpirun],[ff-mpirun],[no])
])

AC_ARG_VAR(MPIRUN,[MPIRUN command])
AC_CHECK_PROGS(MPIRUN, mpirun mpiexec srun ccc_mprun, $MPIRUN)


AC_ARG_WITH(FreeFEM,
	[  --with-FreeFEM=         FreeFEM  installation path : (e.g. --with-FreeFEM=/opt/FreeFem/bin)],
	   FF="${withval}/FreeFem++"
	   FFCXX="${withval}/ff-c++"
	   FFCXXDEP="${withval}/ff-get-dep"
	   FFMPI="${withval}/ff-mpirun")

AC_ARG_WITH([hdf5],
	[  --with-hdf5=               hdf5 root path : (e.g. --with-hdf5=/home/salome/INSTALL/hdf5)],
        [hdf5="$withval"], [hdf5="no"])

AC_ARG_WITH([medfile],
	[  --with-medfile=            medfile root path : (e.g. --with-medfile=/home/salome/INSTALL/medfile)],
        [medfile="$withval"], [medfile="no"])

AC_ARG_WITH([medcoupling],
	[  --with-medcoupling=        medcoupling root path : (e.g. --with-salome=/home/salome/INSTALL/MEDCOUPLING)],
        [medcoupling="$withval"], [medcoupling="no"])


AC_ARG_WITH([dependencies],
	[  --with-dependencies=         Dependencies installation : (e.g. --with-dependencies=yes)],
        [dependencies="$withval"], [dependencies="no"])


AC_SUBST([dependencies])
AM_CONDITIONAL([BUILDDEPENDENCIES], [test $dependencies != no])
test "$dependencies" == "no" && AC_MSG_WARN([

 PDMT can automatically compile its dependencies
 MedCoupling, Medfile, Hdf5 for you.
 In order to do so

 --with-dependencies=yes

 The wget command will be used  to download these,
 if you wish to provided already downloaded files,
 simply add the MEDCOUPLING-$$VERSION-MPI.tar.gz
 in /ext folder.
 ])



AU_ALIAS([ACX_MPI], [AX_MPI])
AC_DEFUN([AX_MPI], [
AC_PREREQ([2.50])

AC_LANG_CASE(
[C], [
	AC_REQUIRE([AC_PROG_CC])
	AC_ARG_VAR(MPICC,[MPI C compiler command])
	AC_CHECK_PROGS(MPICC, mpicc hcc mpxlc_r mpxlc mpcc cmpicc, $CC)
	ax_mpi_save_CC="$CC"
	CC="$MPICC"
	AC_SUBST(MPICC)
],
[C++], [
	AC_REQUIRE([AC_PROG_CXX])
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	AC_CHECK_PROGS(MPICXX, mpic++ mpicxx mpiCC hcp mpxlC_r mpxlC mpCC cmpic++, $CXX)
	ax_mpi_save_CXX="$CXX"
	CXX="$MPICXX"
	AC_SUBST(MPICXX)
],
[Fortran 77], [
	AC_REQUIRE([AC_PROG_F77])
	AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
	AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf_r mpxlf mpf77 cmpifc, $F77)
	ax_mpi_save_F77="$F77"
	F77="$MPIF77"
	AC_SUBST(MPIF77)
],
[Fortran], [
	AC_REQUIRE([AC_PROG_FC])
	AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
	AC_CHECK_PROGS(MPIFC, mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, $FC)
	ax_mpi_save_FC="$FC"
	FC="$MPIFC"
	AC_SUBST(MPIFC)
])

if test x = x"$MPILIBS"; then
	AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
		[Fortran 77], [AC_MSG_CHECKING([for MPI_Init])
			AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])],
		[Fortran], [AC_MSG_CHECKING([for MPI_Init])
			AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
				AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])])
fi
AC_LANG_CASE([Fortran 77], [
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
	fi
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpich, MPI_Init, [MPILIBS="-lfmpich"])
	fi
],
[Fortran], [
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
	fi
	if test x = x"$MPILIBS"; then
		AC_CHECK_LIB(mpichf90, MPI_Init, [MPILIBS="-lmpichf90"])
	fi
])
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
fi
if test x = x"$MPILIBS"; then
	AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl We have to use AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[]])],[],[]) and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpi.h])
	AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[Fortran 77], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpif.h])
	AC_TRY_COMPILE([#include <mpif.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi],
[Fortran], [if test x != x"$MPILIBS"; then
	AC_MSG_CHECKING([for mpif.h])
	AC_TRY_COMPILE([#include <mpif.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
		AC_MSG_RESULT(no)])
fi])

AC_LANG_CASE([C], [CC="$ax_mpi_save_CC"],
	[C++], [CXX="$ax_mpi_save_CXX"],
	[Fortran 77], [F77="$ax_mpi_save_F77"],
	[Fortran], [FC="$ax_mpi_save_FC"])

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])

SECTION_TITLE([ checking for a C++/MPIC++ compilers ])
AC_PROG_CXX
AC_LANG(C++)
AX_MPI

SECTION_TITLE([ checking for a C/MPIC compilers ])
AC_PROG_CC
AC_LANG(C)
AX_MPI


AC_ARG_WITH(MPICXX,
	[  --with-MPICXX=          path of the MPI C++ compiler to be used
	                  (e.g. --with-MPICXX=your/MPIC++/Install/path/bin/mpic++)],
	   MPICXX="${withval}")

AC_ARG_WITH(MPICC,
	[  --with-MPICC=           path of the MPI C compiler to be used
	                  (e.g. --with-MPICC=your/MPIC++/Install/path/bin/mpicc)],
	   MPICC="${withval}")


AC_ARG_WITH(MPIRUN,
	[  --with-MPIRUN=          path of the MPIRUN to be used
	                  (e.g. --with-MPIRUN=your/Install/path/bin/mpirun)],
	   MPIRUN="${withval}")


test "$MPICXX" == "$CXX" && AC_MSG_ERROR([

  Required  MPI C++ compilers were not  found. If you are
  sure that MPI C++ compiler is  present on you operating
  system you might use the ./configure flag:

    '--with-MPICXX=your/MPIC++/Install/path/bin/mpic++'

  to provide the location of MPI C++ compiler.
                                   ])

test "$MPICC" == "$CC" && AC_MSG_ERROR([

  Required  MPI C compilers were not  found. If you are
  sure that MPI C compiler is  present on you operating
  system you might use the ./configure flag:

    '--with-MPICC=your/MPIC/Install/path/bin/mpicc'

  to provide the location of MPI C compiler.
                                   ])



SECTION_TITLE([      checking for a FF++ compilers     ])
AX_FF
test "$FF" == "no" && AC_MSG_ERROR([
  Required  program  'FreeFem++'  not  found.
  If 'FreeFem++' is not  present in  /usr/bin
  or in  /usr/local/bin  consider  using  the
  '--with-FreeFEM=your/FreeFEM/path/bin' flag
  with ./configure command.
                                   ])

AX_FFCXX
test "$FFCXX" == "no" && AC_MSG_ERROR([
  Required  program   'ff-++'   not    found.
  If  'ff-++'  is  not  present  in  /usr/bin
  or  in  /usr/local/bin consider  using  the
  '--with-FreeFEM=your/FreeFEM/path/bin' flag
  with ./configure command.
                                   ])

AX_FFCXXDEP
test "$FFCXXDEP" == "no" && AC_MSG_ERROR([
  Required  program   'ff-get-dep'   not    found.
  If  'ff-get-dep'  is  not  present  in  /usr/bin
  or  in  /usr/local/bin consider  using  the
  '--with-FreeFEM=your/FreeFEM/path/bin' flag
  with ./configure command.
                                   ])

AX_FFMPI
test "$FFMPI" == "no" && AC_MSG_ERROR([
  Required  program  'ff-mpirun'  not  found.
  If 'ff-mpirun' is not present  in  /usr/bin
  or in  /usr/local/bin  consider  using  the
  '--with-FreeFEM=your/FreeFEM/path/bin' flag
  with ./configure command.
                                  ])



dnl ---------------------------------------------------------------
dnl Did the user provide salome/medcoupling
dnl ---------------------------------------------------------------
SECTION_TITLE([      checking for a med support     ])

AC_SUBST([hdf5])
AC_SUBST([medfile])
AC_SUBST([medcoupling])

AM_CONDITIONAL([HAVEMED], [test $hdf5 != no -a $medfile != no -a $medcoupling != no])
AC_MSG_RESULT([])
test $hdf5 != no && AC_MSG_RESULT([checking for hdf5... yes])
test $hdf5 == no && AC_MSG_RESULT([checking for hdf5... no])
test $medfile != no && AC_MSG_RESULT([checking for medfile... yes])
test $medfile == no && AC_MSG_RESULT([checking for medfile... no])
test $medcoupling != no && AC_MSG_RESULT([checking for medcoupling... yes])
test $medcoupling == no && AC_MSG_RESULT([checking for medcoupling... no])
test $hdf5 != no -a $medfile != no -a $medcoupling != no && AC_MSG_RESULT([med support... yes])
AC_MSG_RESULT([])

dnl CXXFLAGS_MED="-I${hdf5}/include -I${medfile}/include -I${medcoupling}/include  ${CXXFLAGS}"
dnl LDFLAGS_MED="-L${hdf5}/lib -L${medfile}/lib -L${medcoupling}/lib ${LDFLAGS}"
dnl LIBS_MED="-lhdf5 -lmedC -lmedloader -lmedcoupling -linterpkernel -lmedcouplingremapper   ${LIBS}"

test $hdf5 != no -a $medfile != no -a $medcoupling != no && AC_SUBST(CXXFLAGS_MED,"-I${hdf5}/include -I${medfile}/include -I${medcoupling}/include  ${CXXFLAGS}"  )
test $hdf5 != no -a $medfile != no -a $medcoupling != no && AC_SUBST(LDFLAGS_MED,"-L${hdf5}/lib -L${medfile}/lib -L${medcoupling}/lib ${LDFLAGS}"  )
test $hdf5 != no -a $medfile != no -a $medcoupling != no && AC_SUBST(LIBS_MED,"-lhdf5 -lmedC -lmedloader -lmedcoupling -linterpkernel -lmedcouplingremapper   ${LIBS}"  )



test $dependencies != no && AC_SUBST(CXXFLAGS_MED,"-I$(pwd)/ext/MEDCOUPLING-9.11.0-MPI/BINARIES-CO7/hdf5/include -I$(pwd)/ext/MEDCOUPLING-9.11.0-MPI/BINARIES-CO7/medfile/include -I$(pwd)/ext/MEDCOUPLING-9.11.0-MPI/INSTALL/MEDCOUPLING/include ${CXXFLAGS}"  )
test $dependencies != no && AC_SUBST(LDFLAGS_MED,"-L$(pwd)/ext/MEDCOUPLING-9.11.0-MPI/BINARIES-CO7/hdf5/lib -L$(pwd)/ext/MEDCOUPLING-9.11.0-MPI/BINARIES-CO7/medfile/lib -L$(pwd)/ext/MEDCOUPLING-9.11.0-MPI/INSTALL/MEDCOUPLING/lib ${LDFLAGS}"  )
test $dependencies != no && AC_SUBST(LIBS_MED,"-lhdf5 -lmedC -lmedloader -lmedcoupling -linterpkernel -lmedcouplingremapper   ${LIBS}"  )

test $dependencies != no && AC_SUBST(PDMT_MEDCOUPLING_DIR,"$(pwd)/ext/MEDCOUPLING-9.11.0-MPI/INSTALL/MEDCOUPLING"  )
test $dependencies != no && AC_SUBST(PDMT_HDF5_DIR,"$(pwd)/ext/MEDCOUPLING-9.11.0-MPI/BINARIES-CO7/hdf5"  )
test $dependencies != no && AC_SUBST(PDMT_MEDFILE_DIR,"$(pwd)/ext/MEDCOUPLING-9.11.0-MPI/BINARIES-CO7/medfile"  )

AC_MSG_RESULT([])

AC_SUBST(CXXFLAGS_MED,[$CXXFLAGS_MED])
AC_SUBST(LDFLAGS_MED,[$LDFLAGS_MED])
AC_SUBST(LIBS_MED,[$LIBS_MED])


AC_MSG_RESULT([ CXXFLAGS_MED FLAGS    : '$CXXFLAGS_MED'])
AC_MSG_RESULT([ LDFLAGS_MED  FLAGS    : '$LDFLAGS_MED'])
AC_MSG_RESULT([ LIBS_MED     FLAGS    : '$LIBS_MED'])

AC_MSG_RESULT([])

test $hdf5 == no -o $medfile == no -o $medcoupling == no && AC_MSG_WARN([
PDMT will is not configured with med support,
you did not provide the flags:
    --with-hdf5
    --with-medfile
    --with-medcoupling
 ])
dnl ---------------------------------------------------------------


dnl ---------------------------------------------------------------
dnl Associate appropriate suffix for different os
dnl ---------------------------------------------------------------
SECTION_TITLE([      checking for operating system type    ])

os_uname=`uname`
case $os_uname in
     CYGWIN*|MINGW*|MSYS_NT*)
	os_type="windows";;
     Darwin)
	os_type="mac";;
     Linux)
	os_type="linux";;
     *)
	os_type="unknown";;
esac

AC_MSG_RESULT([])
test "$os_type" == "windows" && AC_MSG_RESULT([ OS identified: Windows])
test "$os_type" == "mac"     && AC_MSG_RESULT([ OS identified: MacOs])
test "$os_type" == "linux"   && AC_MSG_RESULT([ OS identified: Linux])

test "$os_type" == "unknown" && AC_MSG_ERROR([
  For the moment PDMT is only built on Linux or MacOS
  please contact  mohd-afeef.badri@cea.fr  for  help
  regarding to building PDMT on your OS.
                                   ])

AM_CONDITIONAL([WINDOWS], [test "$os_type" == "windows"])
AM_CONDITIONAL([MAC],     [test "$os_type" == "mac"])
AM_CONDITIONAL([LINUX],   [test "$os_type" == "linux"])

test "$os_type" == "windows" && AC_SUBST(SO_SUFFIX,"dll"  )
test "$os_type" == "mac"     && AC_SUBST(SO_SUFFIX,"dylib")
test "$os_type" == "linux"   && AC_SUBST(SO_SUFFIX,"so"   )

AC_MSG_RESULT([ Shared object sufix   : '$SO_SUFFIX'])
dnl ---------------------------------------------------------------

SECTION_TITLE([    generating the required makefiles   ])


AC_CONFIG_FILES([
 Makefile
 ext/Makefile
 src/Makefile
 test/Makefile
])
AC_OUTPUT

SECTION_TITLE([               status               ])

dnl Print results
AC_MSG_RESULT([])
AC_MSG_RESULT([ PDMT version $VERSION configured successfully.])
AC_MSG_RESULT([])
AC_MSG_RESULT([ install  path      : '$prefix'])
AC_MSG_RESULT([])
AC_MSG_RESULT([ C++ compiler       : '$CXX'])
AC_MSG_RESULT([ C compiler         : '$CC'])
AC_MSG_RESULT([ C/C++  flags       : '$CXXFLAGS'])
AC_MSG_RESULT([])
AC_MSG_RESULT([ MPIRUN             : '$MPIRUN'])
AC_MSG_RESULT([ MPIC++ compiler    : '$MPICXX'])
AC_MSG_RESULT([ MPIC   compiler    : '$MPICC'])
AC_MSG_RESULT([])
AC_MSG_RESULT([ FreeFEM compiler   : '$FF'])
test "$LDFLAGS" != "" && AC_MSG_RESULT([ Building with '$LDFLAGS' for linker flags.])
test "$LIBS"    != "" && AC_MSG_RESULT([ Building with '$LIBS' libraries.])
AC_MSG_RESULT([])
test "$dependencies"          == "no"     &&  AC_MSG_RESULT([ Dependencies  compilation  : DISABLE])
test "$dependencies"          != "no"     &&  AC_MSG_RESULT([ Dependencies  compilation  : ENABLED])
AC_MSG_RESULT([])
test "$hdf5"         != "no"                         						&&  AC_MSG_RESULT([ Hdf5        : YES $hdf5])
test "$hdf5"         == "no" -a "$dependencies" == "no"					&&  AC_MSG_RESULT([ Hdf5        : NO])
test "$dependencies" != "no"					                          &&  AC_MSG_RESULT([ Hdf5        : YES $PDMT_HDF5_DIR])
test "$medfile"      != "no" 						                        &&  AC_MSG_RESULT([ medfile     : YES $medfile])
test "$medfile"      == "no" -a "$dependencies" == "no"					&&  AC_MSG_RESULT([ medfile     : NO])
test "$dependencies" != "no" 						                        &&  AC_MSG_RESULT([ medfile     : YES $PDMT_MEDFILE_DIR])
test "$medcoupling"  != "no" 						                        &&  AC_MSG_RESULT([ medcoupling : YES $medcoupling])
test "$medcoupling"  == "no" -a "$dependencies" == "no"				  &&  AC_MSG_RESULT([ medcoupling : NO])
test "$dependencies" != "no" 						                        &&  AC_MSG_RESULT([ medcoupling : YES $PDMT_MEDCOUPLING_DIR])
AC_MSG_RESULT([])
SECTION_TITLE([               END               ])
