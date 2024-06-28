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

using namespace std;

int PdmtHelp()
{
  PdmtVersion();
    cout << R"(
===================================================================
                            Command Line Parameters
===================================================================

  --version       : Print version info
  --help          : Print help info
  --debug         : Print verbose info
  --mesh          : Provide input mesh
                    (accepts .mesh, .msh, .vtk, .med)
                    Also accepts ("square" or "circle")
  --out_mesh      : Provide name for saved mesh
                    (accepts .med, .vtu, .vtk, .typ2)
  --med_mesh_name : Provide name of mesh in MED file

===================================================================
                            Usage Examples
===================================================================

  # Mesh unit square and print debug info
  PDMT --debug --mesh square

  # Mesh unit circle and print debug info
  PDMT --debug --mesh circle

  # Mesh unit square and save it as out.typ2
  PDMT --mesh square --out_mesh out.typ2

  # Mesh input triangulation in.mesh and save as out.med
  PDMT --mesh ./in.mesh --out_mesh out.med

  # Mesh input triangulation in.msh and save as out.vtu
  PDMT --mesh ./in.msh --out_mesh out.vtu

  # Mesh input .vtk file and save as out.vtk with verbose info
  PDMT --debug --mesh ./in.vtk --out_mesh out.vtk

  # Mesh input .med file and save as out.typ2 with verbose info
  PDMT --debug --mesh ./in.med --out_mesh out.typ2

  # Provide MED mesh name
  PDMT --mesh ./in.med --med_mesh_name my_mesh --out_mesh out.vtu

===================================================================
)" << endl;
  return 0;
}
