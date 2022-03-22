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

#ifdef DMEDCOUPLING
#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingMemArray.hxx"
using namespace MEDCoupling;
#endif

using namespace std;


template<class K>
class polyMeshWrite_Op : public E_F0mps
{
public:
    Expression filename			                ;

    static const int n_name_param = 2		        ;
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
    {"cells"   , &typeid(KN<KN<long>>*)}
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

    if(nargs[0])
      nodesPoly = GetAny< KNM< double > * >((*nargs[0])(stack));

    if(nargs[1])
     CellsPoly = GetAny< KN<KN<long>> * >((*nargs[1])(stack));

    if(verbosity){
     std::cout << "------------------------------------------------------ " <<  std::endl;
     std::cout << "PMDT will now write the polyhedral mesh " << *inputfile  << std::endl;
     std::cout << "------------------------------------------------------ " <<  std::endl;
    }


    std::string fullFileName(*inputfile);

    //std::size_t writeVtkFile = (fullFileName).find(".vtk");
    //std::size_t writeMedFile = (fullFileName).find(".vtk");

    if ( (fullFileName).find(".vtk") != std::string::npos){

      ofstream polyWrite;
      polyWrite.open(*inputfile);

      //------------ Write header ----------------//

      polyWrite << "# vtk DataFile Version 2.0\n"
                << "Unstructured Grid PDMT\n"
                << "ASCII\n"
                << "DATASET UNSTRUCTURED_GRID\n\n";


      //------------ Write nodes ----------------//

      if(verbosity){
       std::cout << "------------------------------------------------------ " <<  std::endl;
       std::cout << "PMDT NODES IN POLYMESH " << nodesPoly->N()  << std::endl;
       std::cout << "------------------------------------------------------ " <<  std::endl;
      }

      {
      int TotalNodes = nodesPoly->N();

      polyWrite << "POINTS " << TotalNodes << " float\n";
      for(int i=0; i < TotalNodes; i++)
        polyWrite << (*nodesPoly)(i,0) << "\t" << (*nodesPoly)(i,1) << "\t 0\n";
      }

      //------------ Write cells ----------------//
      if(verbosity){
       std::cout << "------------------------------------------------------ " <<  std::endl;
       std::cout << "PMDT CELS IN POLYMESH " << CellsPoly->N()  << std::endl;
       //std::cout << "PMDT CELS IN POLYMESH " << (*CellsPoly)(0).N()  << std::endl;
       //std::cout << "PMDT CELS IN POLYMESH " << (*CellsPoly)(0)(0)  << std::endl;
       std::cout << "------------------------------------------------------ " <<  std::endl;
      }

      {
        int TotalCells = CellsPoly->N();
        int TotalCellConnectivity = 0;

        for(int i=0; i < TotalCells; i++)
           TotalCellConnectivity += (*CellsPoly)(i).N() + 1;

        polyWrite << "\n";
        polyWrite << "CELLS " << TotalCells << "\t" << TotalCellConnectivity;
        polyWrite << "\n";

        for(int i=0; i < TotalCells; i++){
          polyWrite << (*CellsPoly)(i).N() << " ";
          for(int j=0; j < (*CellsPoly)(i).N(); j++){
            polyWrite << (*CellsPoly)(i)(j) << " ";
          }
          polyWrite << "\n" ;
        }
      }

      //------------ Write cell types -----------------------//
      {
        int TotalCells = CellsPoly->N();
        polyWrite << "\n";
        polyWrite << "\n";
        polyWrite << "CELL_TYPES " << CellsPoly->N() << " " ;
        polyWrite << "\n";
        for(int i=0; i < TotalCells; i++)
          polyWrite << "7" << "\n";
      }

    }


#ifdef DMEDCOUPLING
    if ( (fullFileName).find(".med") != std::string::npos)
    {

      int TotalNodes = nodesPoly->N();
      double targetCoords[TotalNodes*2];
      for(int i=0; i < TotalNodes; i++){
        targetCoords[i*2]=(*nodesPoly)(i,0) ;
        targetCoords[i*2+1]=(*nodesPoly)(i,1);
      }


      int TotalCells = CellsPoly->N();
      int TotalCellConnectivity = 0;

      for(int i=0; i < TotalCells; i++)
        TotalCellConnectivity += (*CellsPoly)(i).N() + 1;

      mcIdType targetConn[TotalCellConnectivity-TotalCells];

      int count=0;
      for(int i=0; i < TotalCells; i++){
        for(int j=0; j < (*CellsPoly)(i).N(); j++){
          targetConn[count]= (*CellsPoly)(i)(j);
          count++;
        }
      }

      MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
      targetMesh->setMeshDimension(2);
      targetMesh->allocateCells(TotalCells);
      targetMesh->setName("2DPolyMesh_2");

      count=0;
      for(int i=0; i < TotalCells; i++){
        targetMesh->insertNextCell(INTERP_KERNEL::NORM_POLYGON,(*CellsPoly)(i).N(),targetConn+count);
        count += (*CellsPoly)(i).N();
      }

      targetMesh->finishInsertingCells();
      DataArrayDouble *myCoords=DataArrayDouble::New();
      myCoords->alloc(TotalNodes,2);
      myCoords->setInfoOnComponent(0,"x [m]");
      myCoords->setInfoOnComponent(1,"y [m]");
      std::copy(targetCoords,targetCoords+(TotalNodes*2),myCoords->getPointer());
      targetMesh->setCoords(myCoords);
      myCoords->decrRef();
      WriteUMesh(*inputfile,targetMesh,true);

      }
#endif

    return 0L;
};
