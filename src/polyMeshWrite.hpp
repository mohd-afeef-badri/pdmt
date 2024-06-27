/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.com
     Date     : 17/01/2022
     Comment  : The program is used for writing polyhedral mesh

     -------------------------------------------------------------------

     PDMT a parallel  dual meshing tool uses   finite  element framework
     to convert a triangular / tetrahedral mesh into a  polyhedral  mesh.
     PDMT is distributed  in  the  hope that it  will be useful, HOWEVER
     WITHOUT ANY WARRANTY; or without  even  implied warranty of FITNESS
     FOR A PARTICULAR PURPOSE.

*******************************************************************************/

#include <iostream>
#include <stack>
#include <stdlib.h>
#include <set>

#include <mpi.h>

#include "typWriter.hpp"
#include "vtkWriter.hpp"

#ifdef MEDCOUPLING
#include "MEDLoader.hxx"
#include "MEDFileData.hxx"
/*
#include "MEDLoaderBase.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingMemArray.hxx"

*/
using namespace MEDCoupling;
#include "medWriter.hpp"
#endif

using namespace std;


template<class K>
class polyMeshWrite_Op : public E_F0mps
{
public:
    Expression filename			                ;

    static const int n_name_param = 4		        ;
    static basicAC_F0::name_and_type name_param[]	;
    Expression nargs[n_name_param]			;

    polyMeshWrite_Op(const basicAC_F0& args		,
                      Expression param1
                     ) :
        filename     (param1)
    {
        args.SetNameParam(n_name_param	,
                          name_param	,
                          nargs
                         )		;
    }

    AnyType operator()(Stack stack) const		;
};

template<class K>
basicAC_F0::name_and_type polyMeshWrite_Op<K>::name_param[] =
{
    {"nodes"   , &typeid(KNM<double>*)},
    {"cells"   , &typeid(KN<KN<long>>*)},
    {"edges"   , &typeid(KN<KN<long>>*)},
    {"labels"  , &typeid(KN<long>*)}
};


template<class K>
class polyMeshWrite : public OneOperator
{
public:
    polyMeshWrite() : OneOperator(atype<long>()	,
                                       atype<string*>()
                                      ) {}

    E_F0* code(const basicAC_F0& args) const
    {
        return new polyMeshWrite_Op<K>(args,
                                        t[0]->CastTo(args[0])
                                       );
    }
};

template<class K>
AnyType polyMeshWrite_Op<K>::operator()(Stack stack) const
{
    string* inputfile  = GetAny<string*>((*filename)(stack))	;

    KNM < double > *nodesPoly = 0;
    KN<KN<long>>   *CellsPoly = 0;
    KN<KN<long>>   *EdgesPoly = 0;
    KN<long>       *LabelsPoly= 0;

    bool withEdges = false;
    bool withLabel = false;


    if(nargs[0])
      nodesPoly = GetAny< KNM< double > * >((*nargs[0])(stack));

    if(nargs[1])
     CellsPoly = GetAny< KN<KN<long>> * >((*nargs[1])(stack));

    if(nargs[2]){
     EdgesPoly = GetAny< KN<KN<long>> * >((*nargs[2])(stack));
     withEdges = true;
    }

    if(nargs[3]){
     LabelsPoly = GetAny< KN<long> * >((*nargs[3])(stack));
     withLabel = true;
    }

    if(verbosity){
     std::cout << "------------------------------------------------------ " <<  std::endl;
     std::cout << "PMDT will now write the polyhedral mesh " << *inputfile  << std::endl;
     std::cout << "------------------------------------------------------ " <<  std::endl;
    }


    std::string fullFileName(*inputfile);

    //std::size_t writeVtkFile = (fullFileName).find(".vtk");
    //std::size_t writeMedFile = (fullFileName).find(".vtk");

    if ( (fullFileName).find(".vtk") != std::string::npos){
     writePolyVtk(inputfile,nodesPoly,CellsPoly,EdgesPoly, LabelsPoly);
    }

    if ( (fullFileName).find(".vtu") != std::string::npos)
    {

     int timePvd, mpiSize;
     string basVtuFileName = "";

     parallelIO(inputfile, 0, &timePvd, &mpiSize, &basVtuFileName);
     writePolyVtu(inputfile,nodesPoly,CellsPoly,EdgesPoly, LabelsPoly);

     if( mpiSize > 1 )
       PvtuWriter(inputfile,mpiSize,timePvd,basVtuFileName);
    }


#ifdef MEDCOUPLING
    if ((fullFileName).find(".med") != std::string::npos) {
     writePolyMed(inputfile,nodesPoly,CellsPoly,EdgesPoly, LabelsPoly);
    }
#endif

    return 0L;
};
