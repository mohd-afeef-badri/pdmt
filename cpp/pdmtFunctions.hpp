/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.com
     Date     : 22/02/2022
     Comment  : The program collects functions needed by the pdmt mesher

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

using namespace std;

//-------------------------------------------------------------------------
// PdmtMarkBorderNodes : This function takes 2D mesh Th and a vector 
//                       isBorderNode as input. At the end :
//                         \forall i=[0:Th.nv]
//                           isBorderNode(i) = 0 (if i is not on border)
//                           isBorderNode(i) > 0 (if i is on border)
//-------------------------------------------------------------------------
long PdmtMarkBorderNodes(pmesh pTh, KN< long > *isBorderNode) { 

  const Mesh &Th = *pTh;

#ifdef DEBUG
  cout << "\n"
          "--------------------------------------\n"
          " PdmtMarkBorderNodes function called  \n"
          "--------------------------------------\n";
      
  cout << "  Th.nv  - mesh vertices -    "  << Th.nv  << "\n" 
          "  Th.nbe - mesh edges    -    "  << Th.neb << "\n";
          
  cout << "--------------------------------------\n" << endl;  
#endif 
  
  isBorderNode->resize(Th.nv);                  // resize isBorderNode  
  
  for(int i = 0; i < Th.nv; ++i)                // isBorderNode = 0
    *(isBorderNode[0]+i) = 0;
  
  for(int k = 0; k < Th.neb; ++k)               // isBorderNode > 0
   for(int i = 0; i < 2; ++i)
     *(isBorderNode[0]+Th(Th.be(k)[i])) += 1; 
     
  return 0; 

}



//-------------------------------------------------------------------------
// PdmtFillSearchTableEdges : 
//           This function takes 2D mesh Th and two vectors headVertexBorder
//           and nextVertexBorder as input. At the end these two vectors can 
//           be used to search and navigate through number of edges for each 
//           mesh point i; 
//-------------------------------------------------------------------------
int PdmtFillSearchTableEdges(const Fem2D::Mesh* const &pTh, KN< long > *const &headVertexBorder,  KN< long > *const &nextVertexBorder) {

  const Mesh &Th = *pTh;
  
#ifdef DEBUG
  cout << "\n"
          "--------------------------------------\n"
          " PdmtFillSearchTableEdges function called  \n"
          "--------------------------------------\n";
      
  cout << "  Th.nv  - mesh vertices -    "  << Th.nv  << "\n" 
          "  Th.nbe - mesh edges    -    "  << Th.neb << "\n";
          
  cout << "--------------------------------------\n" << endl;  
#endif


  headVertexBorder->resize(Th.nv);              // resize headVertexBorder
  nextVertexBorder->resize(Th.neb*2);           // resize nextVertexBorder   
  
  for(int i = 0; i < Th.nv; ++i)                // set headVertexBorder = -1
    *(headVertexBorder[0]+i) = -1; 

  for(int k = 0; k < Th.neb; ++k){
    for(int i = 0; i < 2; ++i){
      int v = Th(Th.be(k)[i]);                                        // current vertex number 
      *(nextVertexBorder[0]+(2*k+i))  =  *(headVertexBorder[0]+v);    // next vertex           
      *(headVertexBorder[0]+v)     = 2*k+i;                           // update head vertex    
    }
  }
 
  return 0;  
}


//-------------------------------------------------------------------------
// PdmtFillSearchTableTriangles : 
//           This function takes 2D mesh Th and two vectors  headVertex  and 
//           and  nextVertex as input. At  the  end these two vectors can be  
//           used to search &  navigate through number of triangles for each 
//           mesh point i; 
//-------------------------------------------------------------------------
int PdmtFillSearchTableTriangles(const Fem2D::Mesh* const &pTh, KN< long > *const &headVertex,  KN< long > *const &nextVertex) {

  const Mesh &Th = *pTh;
  
#ifdef DEBUG
  cout << "\n"
          "--------------------------------------\n"
          " PdmtFillSearchTableTriangles function called  \n"
          "--------------------------------------\n";
      
  cout << "  Th.nv  - mesh vertices -    "  << Th.nv  << "\n" 
          "  Th.nbe - mesh edges    -    "  << Th.neb << "\n"
          "  Th.nt  - mesh triangles-    "  << Th.nt  << "\n";          
          
  cout << "--------------------------------------\n" << endl;  
#endif


  headVertex->resize(Th.nv);              // resize headVertex
  nextVertex->resize(Th.nt*3);            // resize nextVertex   
  
  for(int i = 0; i < Th.nv; ++i)          // set headVertex = -1
    *(headVertex[0]+i) = -1; 

  for(int k = 0; k < Th.nt; ++k){
    for(int i = 0; i < 3; ++i){
      int v = Th(Th[k][i]);                                     // current vertex number 
      *(nextVertex[0]+(3*k+i))  =  *(headVertex[0]+v);          // next vertex           
      *(headVertex[0]+v)        = 3*k+i;                        // update head vertex    
    }
  }
 
  return 0;  
}


double PdmtMaxinTwoP1(KN<double> *const & f, KN<double> *const & f1)
{
  long int nn = f->N();

  for(long int i=0; i<nn; i++)
    *(f[0]+i)=max(*(f1[0]+i),*(f[0]+i));

  return 0.0;
}
