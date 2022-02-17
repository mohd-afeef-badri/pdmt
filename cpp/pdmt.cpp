/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.com
     Date     : 15/01/2022
     Comment  : The program finds convex hull of a set of  points.
                Original  program  found  online  was  adapted for
                PDMT.  For explanation of the orientation()
                www.geeksforgeeks.org/orientation-3-ordered-points

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

#include "ff++.hpp"
#include "pdmt.hpp"

#include "polyMeshWrite.hpp"

using namespace std;

static void InitFF()
{
  Global.Add("pdmtConvexHull" ,"(",new OneOperator3_<int,KN<double>*,KN<double>*,KN<double>*>(PdmtConvexHull));
  Global.Add("PdmtMaxinTwoP1" ,"(",new OneOperator2_<double,KN<double>*, KN<double>*>(PdmtMaxinTwoP1));
  Global.Add("PdmtUnitTest"   ,"(",new OneOperator0<int>(PdmtTest));

  Global.Add("PdmtPolyMeshWrite"  ,"(", new polyMeshWrite<double>);
}
LOADFUNC(InitFF)
