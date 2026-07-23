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

  --version          : Print version info
  --help             : Print help info
  --debug            : Print verbose info
  --mesh             : Provide input mesh
                         accepts: .mesh, .msh, .vtk, .med
                         Also accepts : "square" or "circle"
  --dimension        : Input type: 2 (default), 3, or 3S surface
  --feature_angle    : 3D/3S feature angle in degrees (default: 45)
  --conserve_edge    : Comma-separated Gmsh physical edge-group names
  --mode             : 3D/3S construction: subdivided_dual or smooth_dual
                         defaults: smooth_dual for 3D, subdivided_dual for 3S
  --out_mesh         : Provide name for saved mesh
                         accepts: .med, .vtu, .vtk, .typ2
  --square_mesh_size : Provide mesh size for square mesh
                       works only with '--mesh square'
  --circle_mesh_size : Provide mesh size for circle mesh
                       works only with '--mesh circle'
  --med_mesh_name    : Provide name of mesh in MED file
                       works with MED input '--mesh file.med'
===================================================================
                            Usage Examples
===================================================================

  # Mesh unit square and print debug info
  PDMT --debug --mesh square

  # Mesh unit square with size
  PDMT --mesh square --square_mesh_size 15

  # Mesh unit circle and print debug info
  PDMT --debug --mesh circle

  # Mesh unit circle with size
  PDMT --mesh circle --circle_mesh_size 120

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

  # Convert a tetrahedral mesh to a 3D polyhedral VTU mesh
  PDMT --dimension 3 --mesh ./tetra.mesh --out_mesh polyhedra.vtu

  # Read a named tetrahedral MED mesh and write native MED polyhedra
  PDMT --dimension 3 --mesh ./tetra.med --med_mesh_name TetrahedralMesh \
       --out_mesh polyhedra.med

  # Retain all barycentric subdivision points in a tetrahedral dual
  PDMT --dimension 3 --mode subdivided_dual --mesh ./tetra.mesh \
       --out_mesh subdivided-polyhedra.vtu

  # Convert a triangular surface mesh embedded in 3D
  PDMT --dimension 3S --mesh ./surface.msh --out_mesh polygons.vtu

  # Read a named triangular MED surface and write MED polygons
  PDMT --dimension 3S --mesh ./surface.med --med_mesh_name TriangularMesh \
       --out_mesh polygons.med

  # Connect triangle barycentres directly except at protected edges
  PDMT --dimension 3S --mode smooth_dual --mesh ./surface.msh \
       --conserve_edge ridge --out_mesh smooth-polygons.vtu

  # Preserve named Gmsh physical edges while merging boundary faces
  PDMT --dimension 3 --mesh ./tetra.msh --conserve_edge ridge,corner \
       --out_mesh polyhedra.vtu

===================================================================
)" << endl;
  return 0;
}
