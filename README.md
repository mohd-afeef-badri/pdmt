![image](https://github.com/user-attachments/assets/74f0a22c-2310-41fd-8ca6-302c8d349eaf)

PDMT, acryonm for Parallel Dual Meshing Tool, is a polyhedral meshing/remsehing tool. It harnesses the power of the finite element framework (FreeFEM) to facilitate the seamless transformation of a triangular or tetrahedral mesh into a more versatile and efficient polyhedral mesh. PDMT adeptly identifies the dual structure of the original triangular mesh, thereby establishing a solid foundation for the subsequent creation of the polyhedral mesh. The underlying Voronoi frame is utilized to form this new polyhedral mesh, ensuring optimal utilization of computational resources and enhancing the mesh's adaptability for a diverse range of engineering simulations and scientific analyses.

**2D mesh Example**
![image](https://github.com/mohd-afeef-badri/pdmt/assets/52162083/bc7f98a6-7631-439d-934f-7daa49250721)

**3D mesh Example**
<img width="1335" height="511" alt="flight" src="https://github.com/user-attachments/assets/13170768-71bc-449d-8321-7e02eea5a09d" />

## Dependencies

To compile and use PDMT, you will need the following dependencies:

- FreeFEM
- MedCoupling (optional for .med mesh support)

## Compilation

### Compilation with precompiled MedCoupling (install procedure 1)

Below, we present a step-by-step guide on how to compile and install PDMT with precompiled MedCoupling support for FreeFEM:

##### Prepare the Build Configuration:

Run the following command to prepare the build configuration:

```bash
autoreconf -i
```

##### Configure the Build:

Now, let's configure the build by specifying the necessary options using the `configure` script:

```bash
./configure \
--prefix=/home/Work/tmp/pdmt \
--with-medcoupling=/home/Install/TarPackages/SALOME-9.15.0-native-UB24.04/BINARIES-UB24.04/MEDCOUPLING \
--with-medfile=/home/Install/TarPackages/SALOME-9.15.0-native-UB24.04/BINARIES-UB24.04/medfile         \
--with-hdf5=/home/Install/TarPackages/SALOME-9.15.0-native-UB24.04/BINARIES-UB24.04/hdf5
```

In this configuration:

- PDMT will be installed in the  `/home/Work/tmp/pdmt` directory. You can choose a directory that you wish to install PDMT.
- The root directories for `medcoupling` , `medfile`, and `hdf5` which come as precompiled with SALOME are provided `/home/Install/TarPackages/SALOME-9.15.0-native-UB24.04/BINARIES-UB24.04/`.  To get your precompiled SALOME click [here](https://www.salome-platform.org/?page_id=2433).
- To proceed with the build, ensure that FreeFEM is already installed and available in your  `$PATH`. If not use flag `--with-FreeFEM` to configure with FreeFEM installed elsewhere.

Please note you will need to adapt each flag to your specific system.

##### Compile PDMT:

With the configuration set, it's time to compile PDMT and make sure it is ready for use:

```bash
make
```

##### Install PDMT:

Once the compilation process is successful, proceed to install PDMT using the following command:

```bash
make install
```

By following these steps, you should have successfully compiled and installed PDMT linked with precompiled MedCoupling support for FreeFEM.

### Compilation with MedCoupling (install procedure 2)

Here is a step-by-step guide for the typical compilation process:

##### Prepare the Build Configuration:

Begin by running the following command to set up the build configuration:

```bash
autoreconf -i
```

##### Configure the Build:

Next, use the `configure` script to configure the build with the required options:

```bash
./configure \
--prefix=/home/Work/tmp/pdmt \
--with-dependencies
```

In this configuration:

- PDMT will be installed in the  `/home/Work/tmp/pdmt` directory. You can choose a directory that you wish to install PDMT.
- After successfull install, the root directories for medcoupling, medfile, and hdf5 will be automatically located in `ext/MEDCOUPLING-9.15.0/INSTALL`.
- Ensure that FreeFEM is already installed and available in your `$PATH`. If not use flag `--with-FreeFEM` to configure with FreeFEM installed elsewhere.

##### Compile PDMT:

With the configuration set, it's time to compile PDMT and make sure it is ready for use:

```bash
make
```

##### Install PDMT:

Once the compilation process is successful, proceed to install PDMT using the following command:

```bash
make install
```

By following these steps, you should have successfully compiled and installed PDMT linked with MedCoupling support for FreeFEM.

### Check the compilation

To ensure a successful compilation, you can run the following command to perform checks:

```
make check
```

### Note

to use and run PDMT with med support please make sure that your `$LD_LIBRARY_PATH` variable contains  `medcoupling`, `medfile`, and `hdf5`  paths. For example for the install above

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/Install/TarPackages/SALOME-9.15.0-native-UB24.04/BINARIES-UB24.04/MEDCOUPLING/lib:/home/Install/TarPackages/SALOME-9.15.0-native-UB24.04/BINARIES-UB24.04/medfile/lib:/home/Install/TarPackages/SALOME-9.15.0-native-UB24.04/BINARIES-UB24.04/hdf5/lib
```

By following these steps, you will have successfully compiled and installed PDMT with MedCoupling support for FreeFEM. Enjoy the poly meshes for your computational simulations!

## Usage

After installation is done you can simply launch the PDMT mesh conversion via a TUI in any folder of choice. The list of command line flags it accepts

- `--debug`    : to print some verbos info about the meshing process
- `--mesh`     : to provide mesh for conversion, it accepts .mesh, .msh, .vtk, .med(conditional) formats. Also accepts ("square" or "circle").
- `--dimension`: input mesh type (`2` by default, `3` for tetrahedra, or `3S` for a triangular surface embedded in 3D).
- `--feature_angle`: preserve 3D boundary edges sharper than this angle (`45` degrees by default).
- `--conserve_edge`: comma-separated Gmsh physical edge-group names that must remain as feature-edge chains in the 3D polyhedral boundary.
- `--mode`: dual construction for every dimension, either `subdivided_dual` or `smooth_dual`. The defaults are `smooth_dual` for 2D/3D and `subdivided_dual` for 3S.

![image](https://github.com/mohd-afeef-badri/pdmt/assets/52162083/8ae5798d-5a4f-474d-ae39-c7207085f7bd)
![image](https://github.com/mohd-afeef-badri/pdmt/assets/52162083/03f0e8ae-75dd-4823-870b-4c65fab363fe)
![image](https://github.com/mohd-afeef-badri/pdmt/assets/52162083/9052499a-3993-425e-a111-2f94c4ca8798)

### Example 1:

```
PDMT --debug --mesh /your/mesh/file.mesh
```

### 2D polygonal meshes

Use `--dimension 2` to convert a planar triangular mesh into a polygonal dual mesh. Each input vertex produces one output polygon.

The default `--mode smooth_dual` connects neighboring triangle barycentres directly. This is the original PDMT 2D construction and produces the smallest number of polygon edges:

```bash
PDMT --dimension 2 \
  --mesh square \
  --square_mesh_size 4 \
  --mode smooth_dual \
  --out_mesh square_smooth.vtu
```

With `--mode subdivided_dual`, every dual line passes through the corresponding primal-edge midpoint. Around an interior primal vertex, the polygon therefore alternates between triangle barycentres and edge midpoints. Boundary polygons also retain the original boundary vertex, so the output remains conforming and has the same exterior boundary as the input mesh:

```bash
PDMT --dimension 2 \
  --mesh square \
  --square_mesh_size 4 \
  --mode subdivided_dual \
  --out_mesh square_subdivided.vtu
```

The two modes support the same 2D input and output formats. Use `subdivided_dual` when the dual should explicitly follow the barycentric subdivision, and `smooth_dual` when straighter, less subdivided polygon boundaries are preferred.

### 3D Polyhedral meshes

PDMT can convert a  tetrahedral mesh into a conforming polyhedral mesh. Every input vertex is treated as a seed/generator for polytopal mesh hence produces one polyhedron. Triangle fans around each primal edge are merged into a single polygonal face, matching the normal dual topology instead of exposing the barycentric sub-tetrahedra. Boundary labels and sharp feature edges are retained.

Use `--dimension 3` and write legacy VTK (`.vtk`), XML VTK (`.vtu`), or—when PDMT is built with MEDCoupling—MED (`.med`):

```bash
PDMT --dimension 3 \
	--mesh /your/mesh/tetrahedra.mesh \
  	--out_mesh polyhedra.vtu
```

The default `--feature_angle 45` keeps sharp corners while merging smooth surface regions. Use a larger value for more aggressive merging, for example`--feature_angle 80`; use a smaller value to preserve more boundary edges. Edges separating different boundary labels are always preserved.

Tetrahedral 3D mode supports both dual representations. `smooth_dual` is the default and connects neighboring tetrahedron barycentres directly across ordinary internal subdivisions. `subdivided_dual` retains the intervening tetra-face barycentres and smooth boundary-edge midpoints:

```bash
PDMT --dimension 3 \
  --mesh tetrahedra.msh \
  --mode subdivided_dual \
  --out_mesh subdivided_polyhedra.vtu
```

In either mode, domain-boundary anchors, sharp features, region interfaces, and curves selected by `--conserve_edge` take priority and remain in the polyhedron connectivity. This keeps every polyhedron closed and conforming.

Named physical curve groups in an ASCII Gmsh 2.x input can be protected independently of the angle criterion. For example, to retain the group `bla` from `cylinder_edge.msh`:

```bash
PDMT --dimension 3 \
  --mesh msh/cylinder_edge.msh \
  --conserve_edge bla \
  --out_mesh cylinder_poly.vtu
```

Multiple group names are comma-separated, for example `--conserve_edge inlet_rim,outlet_rim`. Each original curve segment remains on the output boundary as a geometrically identical chain split at the dual edge midpoint.

The 3D loader accepts tetrahedral `.mesh`, `.msh`, `.vtk`, and—when PDMT is built with MED support—`.med` files. Select a non-default MED input mesh with `--med_mesh_name`. MED output uses native `NORM_POLYHED` cells and stores exterior polygonal faces at level `-1`, including boundary family labels.
For example:

```bash
PDMT --dimension 3 --mesh tetrahedra.med \
  --med_mesh_name TetrahedralMesh --out_mesh polyhedra.med
```

Three-dimensional `.typ2` output is not supported. The VTK files use cell type 42 (`VTK_POLYHEDRON`). In VTU output, `pdmt_face_connectivity`,`pdmt_face_offsets`, and `pdmt_face_labels` are stored as field data so polygon topology and boundary tags remain available to downstream tools.

The plugin API is also available directly in FreeFEM:

```freefem
real[int,int] nodes(0,0);
int[int][int] faces, cells;
int[int] cellLabels, faceLabels;

PdmtBuildDual3D(Th3,
  nodes=nodes, faces=faces, cells=cells,
  labels=cellLabels, faceLabels=faceLabels,
  featureAngle=45.0,
  meshFile="tetrahedra.msh", conserveEdge="ridge,corner",
  mode="smooth_dual");

PdmtPolyMeshWrite("polyhedra.vtu",
  nodes=nodes, cells=cells, faces=faces,
  labels=cellLabels, faceLabels=faceLabels);
```

Cell face references use signed, one-based face IDs: the sign records whether the face orientation agrees with the global face connectivity.

### 3D  surface polytopal meshes

Use `--dimension 3S` for a triangular surface embedded in 3D. Each input vertex produces one dual polygon on a smooth part of the surface. A sharp, boundary, region-interface, or explicitly conserved edge splits the polygon fan so the feature remains part of the polygon connectivity.

```bash
PDMT --dimension 3S \
  --mesh surface.msh \
  --feature_angle 45 \
  --conserve_edge ridge,rim \
  --out_mesh surface_dual.vtu
```

The default `--mode subdivided_dual` connects triangle barycentres through primal-edge midpoints and therefore follows the piecewise-triangular surface. With `--mode smooth_dual`, ordinary midpoint vertices are removed and adjacent triangle barycentres are connected by straight polygon edges:

```bash
PDMT --dimension 3S \
  --mesh surface.msh \
  --mode smooth_dual \
  --conserve_edge ridge,rim \
  --out_mesh smooth_surface_dual.vtu
```

Protected geometry has priority in both modes. Boundary edges, edges selected by `--feature_angle`, region interfaces, and `--conserve_edge` curves retain their primal vertices and edge midpoints so those feature segments remain in the output connectivity.

The 3S loader accepts `.mesh`, `.msh`, `.vtk`, and—when MED support is enabled—`.med` triangular surfaces. It writes `.vtk`, `.vtu`, or native MED polygon cells. Use `--med_mesh_name` for the surface mesh name inside an input MED file:

```bash
PDMT --dimension 3S --mesh surface.med \
  --med_mesh_name TriangularMesh --out_mesh surface_dual.med
```

Named edge groups selected with `--conserve_edge` still require an ASCII Gmsh 2.x input, as for the volumetric 3D mode.

The Gmsh file must contain type-2 triangle elements; exporting only the physical curves is not a surface mesh. From a `.geo` file, an ASCII 2.x surface export can be generated with:

```bash
gmsh -2 model.geo -format msh2 -o surface.msh
```

For generating meshes that are presented above you can use the mesh files provided in `msh` folder.  From top left clockwise:

- for the disk mesh with five holes mesh is provided in `.vtk` format

```bash
PDMT --debug --mesh /msh/disk5holes.vtk
```

- for the pentagon mesh is provided in `.mesh` fromat

```bash
PDMT --debug --mesh /msh/pentagon.mesh
```

- for the triangle with one hole at the center mesh is provided in `.med` format

```bash
PDMT --debug --mesh /msh/triangle1hole.med --med_mesh_name Mesh_1
```

- for the disk with large decentered hole  mesh is provided in `.msh` format

```bash
PDMT --debug --mesh /msh/disk1hole.msh
```
