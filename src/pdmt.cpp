/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.com
     Date     : 15/01/2022
     Comment  : PDMT interface for FreeFEM

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
#include "AFunction_ext.hpp"

#include "convexHull.hpp"
#include "polyMeshWrite.hpp"
#include "pdmtFunctions.hpp"
#include "pdmtVersion.hpp"
#include "pdmtHelp.hpp"

using namespace std;
using namespace Fem2D;

static void InitFF()
{

  Global.Add("PdmtVersion"   ,"(",new OneOperator0<int>(PdmtVersion));
  Global.Add("PdmtHelp"   ,"(",new OneOperator0<int>(PdmtHelp));
  Global.Add("pdmtConvexHull" ,"(",new OneOperator3_<int,KN<double>*,KN<double>*,KN<long>*>(PdmtConvexHull));
  Global.Add("PdmtMaxinTwoP1" ,"(",new OneOperator2_<double,KN<double>*, KN<double>*>(PdmtMaxinTwoP1));
  Global.Add("PdmtUnitTest"   ,"(",new OneOperator0<int>(PdmtTest));
  Global.Add("PdmtFillSearchTableTriangles", "(", new OneOperator3_< int,  pmesh, KN< long > *, KN< long > * >(PdmtFillSearchTableTriangles));
  Global.Add("PdmtFillSearchTableEdges", "(", new OneOperator3_< int,  pmesh, KN< long > *, KN< long > * >(PdmtFillSearchTableEdges));
  Global.Add("PdmtGetMeshInfo", "(", new OneOperator3_< int,  pmesh, KN< long > *, KN< long > * >(PdmtGetMeshInfo));
  Global.Add("PdmtPolyMeshWrite"  ,"(", new polyMeshWrite<double>);
}
LOADFUNC(InitFF)
