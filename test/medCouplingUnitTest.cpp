/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.com
     Date     : 21/03/2022
     Comment  : The program uses MedCoupling for loading and writing
                med files (meshes) native to salome

     -------------------------------------------------------------------

     PDMT a parallel  dual meshing tool uses   finite  element framework
     to convert a triangular / tetrahedral mesh into a  polyhedral  mesh.
     PDMT is distributed  in  the  hope that it  will be useful, HOWEVER
     WITHOUT ANY WARRANTY; or without  even  implied warranty of FITNESS
     FOR A PARTICULAR PURPOSE.

*******************************************************************************/

#include <iostream>

#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingMemArray.hxx"

using namespace std;
using namespace MEDCoupling;

/*
int testRun()
{
  double targetCoords[24]={
    -0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7,
    -0.05,0.95, 0.2,1.2, 0.45,0.95
  };
  mcIdType targetConn[24]={1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(5);
  targetMesh->setName("2DMesh_2");
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+12);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+16);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,targetConn+6);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(12,2);
  myCoords->setInfoOnComponent(0,"toto [m]");
  myCoords->setInfoOnComponent(1,"energie [kW]");
  std::copy(targetCoords,targetCoords+24,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  WriteUMesh("file2.med",targetMesh,true);
  return 1;
}
*/

int testRun()
{
double targetCoords[48]={
0.333333	,0.166667,
0.166667	,0.333333,
0.833333	,0.166667,
0.666667	,0.333333,
0.333333	,0.666667,
0.166667	,0.833333,
0.833333	,0.666667,
0.666667	,0.833333,
0	,0,
0.5	,0,
1	,0,
0	,0.5,
1	,0.5,
0	,1,
0.5	,1,
1	,1,
0.25	,0,
0.75	,0,
1	,0.25,
1	,0.75,
0.25	,1,
0.75	,1,
0	,0.25,
0	,0.75
  };


mcIdType targetConn[56]={
/*
5 16 0 1 22 8
6 16 9 17 2 3 0
4 17 10 18 2
6 22 1 4 5 23 11
6 0 3 6 7 4 1
6 2 18 12 19 6 3
4 23 5 20 13
6 4 7 21 14 20 5
5 6 19 15 21 7
*/
16, 0,  1,  22, 8,
16, 9,  17, 2,  3,  0,
17, 10, 18, 2,
22, 1,  4,  5,  23, 11,
0 ,3 , 6  ,7 , 4, 1,
2 ,18, 12 ,19, 6, 3,
23,5,  20 ,13,
4 ,7 , 21 ,14 ,20, 5,
6 ,19, 15 ,21, 7,
22, 8,
8 ,16,
17, 9,
9 ,16,
                 };
  MEDCouplingUMesh *targetMesh2d=MEDCouplingUMesh::New();
  targetMesh2d->setMeshDimension(2);
  targetMesh2d->allocateCells(10);   // total number of cells
  targetMesh2d->setName("2DPolyMesh_2");
  targetMesh2d->insertNextCell(INTERP_KERNEL::NORM_POLYGON,5,targetConn);      // comes from connectivity
  targetMesh2d->insertNextCell(INTERP_KERNEL::NORM_POLYGON,6,targetConn+5);
  targetMesh2d->insertNextCell(INTERP_KERNEL::NORM_POLYGON,4,targetConn+11);
  targetMesh2d->insertNextCell(INTERP_KERNEL::NORM_POLYGON,6,targetConn+15);
  targetMesh2d->insertNextCell(INTERP_KERNEL::NORM_POLYGON,6,targetConn+21);
  targetMesh2d->insertNextCell(INTERP_KERNEL::NORM_POLYGON,6,targetConn+27);
  targetMesh2d->insertNextCell(INTERP_KERNEL::NORM_POLYGON,4,targetConn+33);
  targetMesh2d->insertNextCell(INTERP_KERNEL::NORM_POLYGON,6,targetConn+37);
  targetMesh2d->insertNextCell(INTERP_KERNEL::NORM_POLYGON,5,targetConn+43);
  targetMesh2d->finishInsertingCells();


  MEDCouplingUMesh *targetMesh1d=MEDCouplingUMesh::New();
  targetMesh1d->setMeshDimension(1);
  targetMesh1d->allocateCells(4);
  targetMesh1d->setName("2DPolyMesh_1");
  targetMesh1d->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,targetConn+48);
  targetMesh1d->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,targetConn+50);
  targetMesh1d->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,targetConn+52);
  targetMesh1d->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,targetConn+54);
  targetMesh1d->finishInsertingCells();

  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(24,2); // tottal number of points
  myCoords->setInfoOnComponent(0,"x [m]");
  myCoords->setInfoOnComponent(1,"y [m]");
  std::copy(targetCoords,targetCoords+48,myCoords->getPointer());
  targetMesh2d->setCoords(myCoords);
  targetMesh1d->setCoords(myCoords);

  std::vector<const MEDCouplingUMesh *> finalMesh;
  finalMesh.push_back(targetMesh2d);
  finalMesh.push_back(targetMesh1d);

  WriteUMeshes("finalMesh.med",finalMesh,true);

  myCoords->decrRef();
  WriteUMesh("mesh2d.med",targetMesh2d,true);
  WriteUMesh("mesh1d.med",targetMesh1d,true);
  return 1;
}

int main(){

int err = testRun();
return 1;
}


