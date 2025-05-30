![image](https://github.com/mohd-afeef-badri/pdmt/assets/52162083/bc7f98a6-7631-439d-934f-7daa49250721)

PDMT, acryonm for Parallel Dual Meshing Tool, is a polyhedral meshing/remsehing tool. It harnesses the power of the finite element framework (FreeFEM) to facilitate the seamless transformation of a triangular or tetrahedral mesh into a more versatile and efficient polyhedral mesh. PDMT adeptly identifies the dual structure of the original triangular mesh, thereby establishing a solid foundation for the subsequent creation of the polyhedral mesh. The underlying Voronoi frame is utilized to form this new polyhedral mesh, ensuring optimal utilization of computational resources and enhancing the mesh's adaptability for a diverse range of engineering simulations and scientific analyses.

## Dependencies ##
To compile and use PDMT, you will need the following dependencies:
- FreeFEM
- MedCoupling (optional for .med mesh support)

## Compilation ##

### Compilation with precompiled MedCoupling (install procedure 1)

Below, we present a step-by-step guide on how to compile and install PDMT with precompiled MedCoupling support for FreeFEM:

##### Prepare the Build Configuration: #####
Run the following command to prepare the build configuration:
```bash 
autoreconf -i
```
##### Configure the Build: #####
Now, let's configure the build by specifying the necessary options using the `configure` script:
```bash
./configure \
--prefix=/home/Work/tmp/pdmt \
--with-medcoupling=/home/Install/TarPackages/SALOME-9.14.0-native-UB24.04/BINARIES-UB24.04/MEDCOUPLING \
--with-medfile=/home/Install/TarPackages/SALOME-9.14.0-native-UB24.04/BINARIES-UB24.04/medfile         \
--with-hdf5=/home/Install/TarPackages/SALOME-9.14.0-native-UB24.04/BINARIES-UB24.04/hdf5
```

In this configuration:
-  PDMT will be installed in the  `/home/Work/tmp/pdmt` directory. You can choose a directory that you wish to install PDMT.  
- The root directories for `medcoupling` , `medfile`, and `hdf5` which come as precompiled with SALOME are provided `/home/Install/TarPackages/SALOME-9.14.0-native-UB24.04/BINARIES-UB24.04/`.  To get your precompiled SALOME click [here](https://www.salome-platform.org/?page_id=2433).
- To proceed with the build, ensure that FreeFEM is already installed and available in your  `$PATH`. If not use flag `--with-FreeFEM` to configure with FreeFEM installed elsewhere.

Please note you will need to adapt each flag to your specific system. 

##### Compile PDMT: #####
With the configuration set, it's time to compile PDMT and make sure it is ready for use:
```bash
make
```

##### Install PDMT: #####
Once the compilation process is successful, proceed to install PDMT using the following command:
```bash
make install
```
By following these steps, you should have successfully compiled and installed PDMT linked with precompiled MedCoupling support for FreeFEM.

### Compilation with MedCoupling (install procedure 2)

Here is a step-by-step guide for the typical compilation process:

##### Prepare the Build Configuration: #####

Begin by running the following command to set up the build configuration:
```bash 
autoreconf -i
```

##### Configure the Build: #####
Next, use the `configure` script to configure the build with the required options:
```bash
./configure \
--prefix=/home/Work/tmp/pdmt \
--with-dependencies
```


In this configuration:

-  PDMT will be installed in the  `/home/Work/tmp/pdmt` directory. You can choose a directory that you wish to install PDMT.  
- After successfull install, the root directories for medcoupling, medfile, and hdf5 will be automatically located in `ext/MEDCOUPLING-9.14.0/INSTALL`.
-  Ensure that FreeFEM is already installed and available in your `$PATH`. If not use flag `--with-FreeFEM` to configure with FreeFEM installed elsewhere. 

##### Compile PDMT: #####
With the configuration set, it's time to compile PDMT and make sure it is ready for use:
```bash
make
```

##### Install PDMT: #####
Once the compilation process is successful, proceed to install PDMT using the following command:
```bash
make install
```
By following these steps, you should have successfully compiled and installed PDMT linked with MedCoupling support for FreeFEM.


### Check the compilation ###
To ensure a successful compilation, you can run the following command to perform checks:
```
make check
```

### Note ###

to use and run PDMT with med support please make sure that your `$LD_LIBRARY_PATH` variable contains  `medcoupling`, `medfile`, and `hdf5`  paths. For example for the install above 

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/Install/TarPackages/SALOME-9.14.0-native-UB24.04/BINARIES-UB24.04/MEDCOUPLING/lib:/home/Install/TarPackages/SALOME-9.14.0-native-UB24.04/BINARIES-UB24.04/medfile/lib:/home/Install/TarPackages/SALOME-9.14.0-native-UB24.04/BINARIES-UB24.04/hdf5/lib
```
By following these steps, you will have successfully compiled and installed PDMT with MedCoupling support for FreeFEM. Enjoy the poly meshes for your computational simulations!

## Usage ##

After installation is done you can simply launch the PDMT mesh conversion via a TUI in any folder of choice. The list of command line flags it accepts
- `--debug`    : to print some verbos info about the meshing process
- `--mesh`     : to provide mesh for conversion, it accepts .mesh, .msh, .vtk, .med(conditional) formats. Also accepts ("square" or "circle").



![image](https://github.com/mohd-afeef-badri/pdmt/assets/52162083/8ae5798d-5a4f-474d-ae39-c7207085f7bd)
![image](https://github.com/mohd-afeef-badri/pdmt/assets/52162083/03f0e8ae-75dd-4823-870b-4c65fab363fe)
![image](https://github.com/mohd-afeef-badri/pdmt/assets/52162083/9052499a-3993-425e-a111-2f94c4ca8798)

### Example 1:
```
PDMT --debug --mesh /your/mesh/file.mesh
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



