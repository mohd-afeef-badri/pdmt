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

#ifdef MEDCOUPLING
#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDFileData.hxx"
using namespace MEDCoupling;
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

        int TotalVTKConnectionList = TotalCells;

        for(int i=0; i < TotalCells; i++)
           TotalCellConnectivity += (*CellsPoly)(i).N() + 1;

        if(withEdges){
        int TotalEdges = EdgesPoly->N();
        for(int i=0; i < TotalEdges; i++)
           TotalCellConnectivity += (*EdgesPoly)(i).N() + 1;

        TotalVTKConnectionList +=  TotalEdges;
        }

        polyWrite << "\n";
        polyWrite << "CELLS " << TotalVTKConnectionList << "\t" << TotalCellConnectivity;
        polyWrite << "\n";

        for(int i=0; i < TotalCells; i++){
          polyWrite << (*CellsPoly)(i).N() << " ";
          for(int j=0; j < (*CellsPoly)(i).N(); j++){
            polyWrite << (*CellsPoly)(i)(j) << " ";
          }
          polyWrite << "\n" ;
        }

        if(withEdges){
        for(int i=0; i < EdgesPoly->N(); i++){
          polyWrite << (*EdgesPoly)(i).N() << " ";
          for(int j=0; j < (*EdgesPoly)(i).N(); j++){
            polyWrite << (*EdgesPoly)(i)(j) << " ";
          }
          polyWrite << "\n" ;
        }
        }
      }

      //------------ Write cell types -----------------------//
      {
        int TotalCells = CellsPoly->N();

        if(withEdges)
           TotalCells += EdgesPoly->N();

        polyWrite << "\n";
        polyWrite << "\n";
        polyWrite << "CELL_TYPES " << TotalCells << " " ;
        polyWrite << "\n";
        for(int i=0; i < CellsPoly->N(); i++)
          polyWrite << "7" << "\n";

        if(withEdges)
        for(int i=0; i < EdgesPoly->N(); i++)
          polyWrite << "3" << "\n";
      }


      //------------ Write cell labels -----------------------//
      {
        int TotalCells = CellsPoly->N();

        if(withEdges)
           TotalCells += EdgesPoly->N();

        polyWrite << "\n";
        polyWrite << "\n";
        polyWrite << "CELL_DATA " << TotalCells << " \n" ;
        polyWrite << "SCALARS label float 1 \n";
        polyWrite << "LOOKUP_TABLE CellColors ";
        polyWrite << "\n";

        if(!withLabel)
          for(int i=0; i < CellsPoly->N(); i++)
            polyWrite << "0" << "\n";
        else
          for(int i=0; i < CellsPoly->N(); i++)
            polyWrite << (*LabelsPoly)(i) << "\n";

        if(withEdges && !withLabel)
          for(int i=0; i < EdgesPoly->N(); i++)
            polyWrite << "1" << "\n";

        if(withEdges && withLabel)
          for(int i = CellsPoly->N(); i < EdgesPoly->N() + CellsPoly->N(); i++)
            polyWrite << (*LabelsPoly)(i)  << "\n";
      }
    }

#ifdef MEDCOUPLING
    if ((fullFileName).find(".med") != std::string::npos) {

      //  get nodes of the mesh  //
      int TotalNodes = nodesPoly -> N();

      double medNodeCoords[TotalNodes * 2];

      for (int i = 0; i < TotalNodes; i++) {
        medNodeCoords[i * 2] = ( * nodesPoly)(i, 0);
        medNodeCoords[i * 2 + 1] = ( * nodesPoly)(i, 1);
      }

      //  get cells | cells + edges //
      int TotalCells = CellsPoly -> N();
      int TotalCellConnectivity = 0;
      int TotalConnectionList = TotalCells;

      for (int i = 0; i < TotalCells; i++)
        TotalCellConnectivity += ( * CellsPoly)(i).N() + 1;

      if (withEdges) {
        int TotalEdges = EdgesPoly -> N();

        for (int i = 0; i < TotalEdges; i++)
          TotalCellConnectivity += ( * EdgesPoly)(i).N() + 1;

        TotalConnectionList += TotalEdges;
      }

      mcIdType medCellConn[TotalCellConnectivity - TotalConnectionList];

      int count = 0;
      for (int i = 0; i < CellsPoly -> N(); i++) {
        for (int j = 0; j < ( * CellsPoly)(i).N(); j++) {
          medCellConn[count] = ( * CellsPoly)(i)(j);
          count++;
        }
      }

      if (withEdges) {
        for (int i = 0; i < EdgesPoly -> N(); i++) {
          for (int j = 0; j < ( * EdgesPoly)(i).N(); j++) {
            medCellConn[count] = ( * EdgesPoly)(i)(j);
            count++;
          }
        }
      }

      // --- fill med meshes ---	//
      MEDCouplingUMesh * medMesh2d = MEDCouplingUMesh::New();
      MEDCouplingUMesh * medMesh1d = MEDCouplingUMesh::New();

      medMesh2d -> setMeshDimension(2);
      medMesh2d -> allocateCells(CellsPoly -> N());
      medMesh2d -> setName("PolyMesh");

      count = 0;
      for (int i = 0; i < CellsPoly -> N(); i++) {
        medMesh2d -> insertNextCell(INTERP_KERNEL::NORM_POLYGON, ( * CellsPoly)(i).N(), medCellConn + count);
        count += ( * CellsPoly)(i).N();
      }

      medMesh2d -> finishInsertingCells();

      if (withEdges) {
        medMesh1d -> setMeshDimension(1);
        medMesh1d -> allocateCells(EdgesPoly -> N());
        medMesh1d -> setName("PolyMesh");

        for (int i = 0; i < EdgesPoly -> N(); i++) {
          medMesh1d -> insertNextCell(INTERP_KERNEL::NORM_SEG2, 2, medCellConn + count);
          count += ( * EdgesPoly)(i).N();
        }
      }

      DataArrayDouble * myCoords = DataArrayDouble::New();
      myCoords -> alloc(TotalNodes, 2);
      myCoords -> setInfoOnComponent(0, "x");
      myCoords -> setInfoOnComponent(1, "y");
      std::copy(medNodeCoords, medNodeCoords + (TotalNodes * 2), myCoords -> getPointer());

      medMesh2d -> setCoords(myCoords);

      if (withEdges)
        medMesh1d -> setCoords(myCoords);

      myCoords -> decrRef();

      if (!withLabel) {
        std::vector < const MEDCouplingUMesh * > finalMesh;
        finalMesh.push_back(medMesh2d);
        if (withEdges)
          finalMesh.push_back(medMesh1d);
        WriteUMeshes( * inputfile, finalMesh, true);
      }

      if (withLabel) {

        MCAuto < MEDFileUMesh > finalMeshWithLabel = MEDFileUMesh::New();

        finalMeshWithLabel -> setMeshAtLevel(0, medMesh2d);
        finalMeshWithLabel -> setMeshAtLevel(-1, medMesh1d);

        MCAuto < DataArrayIdType > fam2d = DataArrayIdType::New();
        MCAuto < DataArrayIdType > fam1d = DataArrayIdType::New();

        fam2d -> alloc(CellsPoly -> N(), 1);
        fam1d -> alloc(EdgesPoly -> N(), 1);

        mcIdType elemsFams[CellsPoly -> N() + EdgesPoly -> N()];

        std::set < int > poly2DUniqueLabels;
        std::set < int > poly1DUniqueLabels;

        for (int i = 0; i < CellsPoly -> N(); i++) {
          poly2DUniqueLabels.insert(( * LabelsPoly)(i));
          elemsFams[i] = ( * LabelsPoly)(i) + 1000; // Adding 1000 because med does not like tag zero
        }

        for (int i = CellsPoly -> N(); i < EdgesPoly -> N() + CellsPoly -> N(); i++) {
          poly1DUniqueLabels.insert(( * LabelsPoly)(i));
          elemsFams[i] = ( * LabelsPoly)(i);
        }

        // Iterate through all the elements in a set and display the value.
        for (std::set < int > ::iterator it = poly2DUniqueLabels.begin(); it != poly2DUniqueLabels.end(); ++it)
          std::cout << " Volume " << * it;

        for (std::set < int > ::iterator it = poly1DUniqueLabels.begin(); it != poly1DUniqueLabels.end(); ++it)
          std::cout << " Surface " << * it;

        std::copy(elemsFams, elemsFams + CellsPoly -> N(), fam2d -> getPointer());
        std::copy(elemsFams + CellsPoly -> N(), elemsFams + int(EdgesPoly -> N() + CellsPoly -> N()), fam1d -> getPointer());

        finalMeshWithLabel -> setFamilyFieldArr(-1, fam1d);
        finalMeshWithLabel -> setFamilyFieldArr(0, fam2d);

        std::map < std::string, std::vector < std::string >> theGroups;
        std::map < std::string, mcIdType > theFamilies;

        for (std::set < int > ::iterator it = poly2DUniqueLabels.begin(); it != poly2DUniqueLabels.end(); ++it) {
          theFamilies["cell_family_" + to_string( * it) + ""] = * it + 1000;
          theGroups["cell_group_" + to_string( * it) + ""].push_back("cell_family_" + to_string( * it) + "");
        }

        for (std::set < int > ::iterator it = poly1DUniqueLabels.begin(); it != poly1DUniqueLabels.end(); ++it) {
          theFamilies["boundary_family_" + to_string( * it) + ""] = * it;
          theGroups["boundary_group_" + to_string( * it) + ""].push_back("boundary_family_" + to_string( * it) + "");
        }

    /*
            theFamilies["cells" ]=0;
            theFamilies["border1"]=1;
            theFamilies["border2"]=2;
            theFamilies["border3"]=3;
            theFamilies["border4"]=4;

            theGroups["Face_group"].push_back("cells");
            theGroups["boundary1"].push_back("border1");
            theGroups["boundary2"].push_back("border2");
            theGroups["boundary3"].push_back("border3");
            theGroups["boundary4"].push_back("border4");
    */

        finalMeshWithLabel -> setFamilyInfo(theFamilies);
        finalMeshWithLabel -> setGroupInfo(theGroups);

        finalMeshWithLabel -> write( * inputfile, 2); // med
      }
    }
#endif

    return 0L;
};
