/*****************************************************************************

     Robust MED import for tetrahedral volume and triangular surface meshes.

*****************************************************************************/

#ifndef PDMT_MED_MESH_LOADER_HPP
#define PDMT_MED_MESH_LOADER_HPP

#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"
#include "MEDCouplingUMesh.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace PdmtMedLoader {

inline bool hasLevel(const MEDCoupling::MEDFileUMesh &mesh, int requested) {
  const std::vector<int> levels = mesh.getNonEmptyLevels();
  return std::find(levels.begin(), levels.end(), requested) != levels.end();
}

inline long familyLabel(const MEDCoupling::DataArrayIdType *families,
                        long cell) {
  return families ? static_cast<long>(families->getIJ(cell, 0)) : 1L;
}

inline void ensureMesh(const MEDCoupling::MEDCouplingUMesh *mesh,
                       const char *description) {
  if (!mesh)
    ExecError(description);
  mesh->checkConsistencyLight();
}

inline Mesh3 *loadVolume(const std::string &fileName,
                         const std::string &meshName) {
  using namespace MEDCoupling;
  MCAuto<MEDFileUMesh> fileMesh = MEDFileUMesh::New(fileName, meshName);
  MCAuto<MEDCouplingUMesh> volume =
      dynamic_cast<MEDCouplingUMesh *>(fileMesh->getMeshAtLevel(0));
  ensureMesh(volume, "PdmtLoadMedMesh3: MED level 0 is not an unstructured mesh");
  if (volume->getMeshDimension() != 3 || volume->getSpaceDimension() != 3)
    ExecError("PdmtLoadMedMesh3: MED level 0 must be a volume mesh in 3D");

  const bool suppliedBoundary = hasLevel(*fileMesh, -1);
  MCAuto<MEDCouplingUMesh> boundary;
  if (suppliedBoundary) {
    boundary = dynamic_cast<MEDCouplingUMesh *>(fileMesh->getMeshAtLevel(-1));
    ensureMesh(boundary, "PdmtLoadMedMesh3: MED level -1 is not an unstructured mesh");
    if (boundary->getMeshDimension() != 2)
      ExecError("PdmtLoadMedMesh3: MED level -1 must contain boundary triangles");
  }

  const long nv = volume->getNumberOfNodes();
  const long nt = volume->getNumberOfCells();
  long nbe = suppliedBoundary ? boundary->getNumberOfCells() : 0;
  const DataArrayIdType *cellFamilies = fileMesh->getFamilyFieldAtLevel(0);
  const DataArrayIdType *boundaryFamilies =
      suppliedBoundary ? fileMesh->getFamilyFieldAtLevel(-1) : 0;

  Vertex3 *vertices = new Vertex3[nv];
  const double *coordinates = volume->getCoords()->getConstPointer();
  for (long node = 0; node < nv; ++node) {
    vertices[node].x = coordinates[3 * node];
    vertices[node].y = coordinates[3 * node + 1];
    vertices[node].z = coordinates[3 * node + 2];
    vertices[node].lab = 1;
  }

  Tet *tetrahedra = new Tet[nt];
  std::vector<std::array<int, 4> > tetNodes(nt);
  const mcIdType *connectivity = volume->getNodalConnectivity()->getConstPointer();
  const mcIdType *offsets = volume->getNodalConnectivityIndex()->getConstPointer();
  for (long cell = 0; cell < nt; ++cell) {
    if (volume->getTypeOfCell(cell) != INTERP_KERNEL::NORM_TETRA4)
      ExecError("PdmtLoadMedMesh3: MED level 0 must contain only linear TETRA4 cells");
    const mcIdType begin = offsets[cell] + 1;
    int nodes[4] = {
        static_cast<int>(connectivity[begin]),
        static_cast<int>(connectivity[begin + 1]),
        static_cast<int>(connectivity[begin + 2]),
        static_cast<int>(connectivity[begin + 3])};

    const double abx = vertices[nodes[1]].x - vertices[nodes[0]].x;
    const double aby = vertices[nodes[1]].y - vertices[nodes[0]].y;
    const double abz = vertices[nodes[1]].z - vertices[nodes[0]].z;
    const double acx = vertices[nodes[2]].x - vertices[nodes[0]].x;
    const double acy = vertices[nodes[2]].y - vertices[nodes[0]].y;
    const double acz = vertices[nodes[2]].z - vertices[nodes[0]].z;
    const double adx = vertices[nodes[3]].x - vertices[nodes[0]].x;
    const double ady = vertices[nodes[3]].y - vertices[nodes[0]].y;
    const double adz = vertices[nodes[3]].z - vertices[nodes[0]].z;
    const double determinant =
        abx * (acy * adz - acz * ady) -
        aby * (acx * adz - acz * adx) +
        abz * (acx * ady - acy * adx);
    if (std::abs(determinant) <= 1.e-30)
      ExecError("PdmtLoadMedMesh3: MED input contains a degenerate TETRA4 cell");
    if (determinant < 0.0)
      std::swap(nodes[1], nodes[2]);
    tetNodes[cell] = {{nodes[0], nodes[1], nodes[2], nodes[3]}};
    tetrahedra[cell].set(vertices, nodes, familyLabel(cellFamilies, cell));
  }

  Triangle3 *boundaryTriangles = 0;
  if (suppliedBoundary) {
    boundaryTriangles = new Triangle3[nbe];
    const mcIdType *boundaryConnectivity =
        boundary->getNodalConnectivity()->getConstPointer();
    const mcIdType *boundaryOffsets =
        boundary->getNodalConnectivityIndex()->getConstPointer();
    for (long face = 0; face < nbe; ++face) {
      if (boundary->getTypeOfCell(face) != INTERP_KERNEL::NORM_TRI3)
        ExecError("PdmtLoadMedMesh3: MED level -1 must contain only linear TRI3 cells");
      const mcIdType begin = boundaryOffsets[face] + 1;
      int nodes[3] = {
          static_cast<int>(boundaryConnectivity[begin]),
          static_cast<int>(boundaryConnectivity[begin + 1]),
          static_cast<int>(boundaryConnectivity[begin + 2])};
      boundaryTriangles[face].set(
          vertices, nodes, familyLabel(boundaryFamilies, face));
    }
  } else {
    typedef std::array<int, 3> Face;
    struct BoundaryFaceCandidate {
      Face oriented;
      int count;
      BoundaryFaceCandidate() : oriented({{0, 0, 0}}), count(0) {}
    };
    std::map<Face, BoundaryFaceCandidate> faces;
    static const int localFaces[4][3] = {
        {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    for (long cell = 0; cell < nt; ++cell) {
      for (int local = 0; local < 4; ++local) {
        Face oriented = {{
            tetNodes[cell][localFaces[local][0]],
            tetNodes[cell][localFaces[local][1]],
            tetNodes[cell][localFaces[local][2]]}};
        Face key = oriented;
        std::sort(key.begin(), key.end());
        BoundaryFaceCandidate &candidate = faces[key];
        if (candidate.count == 0)
          candidate.oriented = oriented;
        ++candidate.count;
      }
    }
    for (std::map<Face, BoundaryFaceCandidate>::const_iterator face = faces.begin();
         face != faces.end(); ++face)
      if (face->second.count == 1)
        ++nbe;
    boundaryTriangles = new Triangle3[nbe];
    long boundaryId = 0;
    for (std::map<Face, BoundaryFaceCandidate>::const_iterator face = faces.begin();
         face != faces.end(); ++face) {
      if (face->second.count != 1)
        continue;
      int nodes[3] = {
          face->second.oriented[0],
          face->second.oriented[1],
          face->second.oriented[2]};
      boundaryTriangles[boundaryId++].set(vertices, nodes, 1);
    }
  }

  Mesh3 *result = new Mesh3(nv, nt, nbe, vertices, tetrahedra,
                            boundaryTriangles, false, false,
                            false);
  if (!suppliedBoundary && verbosity)
    std::cout << "PdmtLoadMedMesh3: MED mesh has no level -1; reconstructed "
              << result->nbe << " boundary triangles" << std::endl;
  return result;
}

inline MeshS *loadSurface(const std::string &fileName,
                          const std::string &meshName) {
  using namespace MEDCoupling;
  MCAuto<MEDFileUMesh> fileMesh = MEDFileUMesh::New(fileName, meshName);
  MCAuto<MEDCouplingUMesh> surface =
      dynamic_cast<MEDCouplingUMesh *>(fileMesh->getMeshAtLevel(0));
  ensureMesh(surface, "PdmtLoadMedMeshS: MED level 0 is not an unstructured mesh");
  if (surface->getMeshDimension() != 2 || surface->getSpaceDimension() != 3)
    ExecError("PdmtLoadMedMeshS: MED level 0 must be a surface mesh embedded in 3D");

  const bool suppliedBoundary = hasLevel(*fileMesh, -1);
  MCAuto<MEDCouplingUMesh> boundary;
  if (suppliedBoundary) {
    boundary = dynamic_cast<MEDCouplingUMesh *>(fileMesh->getMeshAtLevel(-1));
    ensureMesh(boundary, "PdmtLoadMedMeshS: MED level -1 is not an unstructured mesh");
    if (boundary->getMeshDimension() != 1)
      ExecError("PdmtLoadMedMeshS: MED level -1 must contain boundary edges");
  }

  const long nv = surface->getNumberOfNodes();
  const long nt = surface->getNumberOfCells();
  const long nbe = suppliedBoundary ? boundary->getNumberOfCells() : 0;
  const DataArrayIdType *cellFamilies = fileMesh->getFamilyFieldAtLevel(0);
  const DataArrayIdType *boundaryFamilies =
      suppliedBoundary ? fileMesh->getFamilyFieldAtLevel(-1) : 0;

  Vertex3 *vertices = new Vertex3[nv];
  const double *coordinates = surface->getCoords()->getConstPointer();
  for (long node = 0; node < nv; ++node) {
    vertices[node].x = coordinates[3 * node];
    vertices[node].y = coordinates[3 * node + 1];
    vertices[node].z = coordinates[3 * node + 2];
    vertices[node].lab = 1;
  }

  TriangleS *triangles = new TriangleS[nt];
  const mcIdType *connectivity = surface->getNodalConnectivity()->getConstPointer();
  const mcIdType *offsets = surface->getNodalConnectivityIndex()->getConstPointer();
  for (long cell = 0; cell < nt; ++cell) {
    if (surface->getTypeOfCell(cell) != INTERP_KERNEL::NORM_TRI3)
      ExecError("PdmtLoadMedMeshS: MED level 0 must contain only linear TRI3 cells");
    const mcIdType begin = offsets[cell] + 1;
    int nodes[3] = {
        static_cast<int>(connectivity[begin]),
        static_cast<int>(connectivity[begin + 1]),
        static_cast<int>(connectivity[begin + 2])};
    triangles[cell].set(vertices, nodes, familyLabel(cellFamilies, cell));
  }

  BoundaryEdgeS *boundaryEdges = suppliedBoundary ? new BoundaryEdgeS[nbe] : 0;
  if (suppliedBoundary) {
    const mcIdType *boundaryConnectivity =
        boundary->getNodalConnectivity()->getConstPointer();
    const mcIdType *boundaryOffsets =
        boundary->getNodalConnectivityIndex()->getConstPointer();
    for (long edge = 0; edge < nbe; ++edge) {
      if (boundary->getTypeOfCell(edge) != INTERP_KERNEL::NORM_SEG2)
        ExecError("PdmtLoadMedMeshS: MED level -1 must contain only linear SEG2 cells");
      const mcIdType begin = boundaryOffsets[edge] + 1;
      int nodes[2] = {
          static_cast<int>(boundaryConnectivity[begin]),
          static_cast<int>(boundaryConnectivity[begin + 1])};
      boundaryEdges[edge].set(
          vertices, nodes, familyLabel(boundaryFamilies, edge));
    }
  }

  MeshS *result = new MeshS(nv, nt, nbe, vertices, triangles, boundaryEdges,
                            false, false, !suppliedBoundary);
  if (!suppliedBoundary && verbosity)
    std::cout << "PdmtLoadMedMeshS: MED mesh has no level -1; reconstructed "
              << result->nbe << " boundary/feature edges" << std::endl;
  return result;
}

inline std::string exceptionMessage(const char *operation,
                                    const std::exception &error) {
  return std::string(operation) + ": " + error.what();
}

} // namespace PdmtMedLoader

class pdmtMedLoader3_Op : public E_F0mps {
public:
  Expression filename;
  static const int n_name_param = 1;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];

  pdmtMedLoader3_Op(const basicAC_F0 &args, Expression input)
      : filename(input) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator()(Stack stack) const {
    std::string *file = GetAny<std::string *>((*filename)(stack));
    if (!nargs[0])
      ExecError("PdmtLoadMedMesh3: meshname is required");
    std::string *meshName = GetAny<std::string *>((*nargs[0])(stack));
    try {
      Mesh3 *mesh = PdmtMedLoader::loadVolume(*file, *meshName);
      Add2StackOfPtr2FreeRC(stack, mesh);
      return mesh;
    } catch (const std::exception &error) {
      const std::string message =
          PdmtMedLoader::exceptionMessage("PdmtLoadMedMesh3", error);
      ExecError(message.c_str());
    }
    return pmesh3(0);
  }
};

basicAC_F0::name_and_type pdmtMedLoader3_Op::name_param[] = {
    {"meshname", &typeid(std::string *)}};

class pdmtMedLoader3 : public OneOperator {
public:
  pdmtMedLoader3() : OneOperator(atype<pmesh3>(), atype<std::string *>()) {}
  E_F0 *code(const basicAC_F0 &args) const {
    return new pdmtMedLoader3_Op(args, t[0]->CastTo(args[0]));
  }
};

class pdmtMedLoaderS_Op : public E_F0mps {
public:
  Expression filename;
  static const int n_name_param = 1;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];

  pdmtMedLoaderS_Op(const basicAC_F0 &args, Expression input)
      : filename(input) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator()(Stack stack) const {
    std::string *file = GetAny<std::string *>((*filename)(stack));
    if (!nargs[0])
      ExecError("PdmtLoadMedMeshS: meshname is required");
    std::string *meshName = GetAny<std::string *>((*nargs[0])(stack));
    try {
      MeshS *mesh = PdmtMedLoader::loadSurface(*file, *meshName);
      Add2StackOfPtr2FreeRC(stack, mesh);
      return mesh;
    } catch (const std::exception &error) {
      const std::string message =
          PdmtMedLoader::exceptionMessage("PdmtLoadMedMeshS", error);
      ExecError(message.c_str());
    }
    return pmeshS(0);
  }
};

basicAC_F0::name_and_type pdmtMedLoaderS_Op::name_param[] = {
    {"meshname", &typeid(std::string *)}};

class pdmtMedLoaderS : public OneOperator {
public:
  pdmtMedLoaderS() : OneOperator(atype<pmeshS>(), atype<std::string *>()) {}
  E_F0 *code(const basicAC_F0 &args) const {
    return new pdmtMedLoaderS_Op(args, t[0]->CastTo(args[0]));
  }
};

#endif
