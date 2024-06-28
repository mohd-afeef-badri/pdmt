/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     -------------------------------------------------------------------

     Author(s): Mohd Afeef Badri
     Email    : mohd-afeef.badri@cea.com
     Date     : 28/06/2024
     Comment  : The program collects functions needed by the PDMT mesher

     -------------------------------------------------------------------

     PDMT a parallel  dual meshing tool uses   finite  element framework
     to convert a triangular / tetrahedral mesh into a  polyhedral  mesh.
     PDMT is distributed  in  the  hope that it  will be useful, HOWEVER
     WITHOUT ANY WARRANTY; or without  even  implied warranty of FITNESS
     FOR A PARTICULAR PURPOSE.

*******************************************************************************/

#include <string>

using namespace std;

int PdmtVersion()
{
  cout << R"(




===================================================================

  PDMT Version 0.2
    Copyright (C) CEA 2022 - 2024


    This is free software; see the source for copying conditions.
    There is NO warranty; not even for MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.


    Report bugs/issues ::    mohd-afeef.badri@cea.fr

===================================================================



)" << endl;
  return 0;
}
