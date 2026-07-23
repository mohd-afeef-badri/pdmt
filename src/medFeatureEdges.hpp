/*****************************************************************************

     Read named edge groups from MED meshes and map them to FreeFEM edges.

*****************************************************************************/

#ifndef PDMT_MED_FEATURE_EDGES_HPP
#define PDMT_MED_FEATURE_EDGES_HPP

#ifdef MEDCOUPLING

#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"
#include "MEDCouplingUMesh.hxx"

#include <limits>
#include <sstream>

namespace Pdmt3D {

template<class MeshType>
inline std::map<long, long> mapMedNodesToMeshVertices(
    const MEDCoupling::MEDCouplingUMesh &medMesh,
    const MeshType &mesh,
    const char *context) {
  double minCoordinate[3] = {
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max()};
  double maxCoordinate[3] = {
      -std::numeric_limits<double>::max(),
      -std::numeric_limits<double>::max(),
      -std::numeric_limits<double>::max()};
  for (long vertex = 0; vertex < mesh.nv; ++vertex) {
    minCoordinate[0] = std::min(minCoordinate[0], mesh(vertex).x);
    minCoordinate[1] = std::min(minCoordinate[1], mesh(vertex).y);
    minCoordinate[2] = std::min(minCoordinate[2], mesh(vertex).z);
    maxCoordinate[0] = std::max(maxCoordinate[0], mesh(vertex).x);
    maxCoordinate[1] = std::max(maxCoordinate[1], mesh(vertex).y);
    maxCoordinate[2] = std::max(maxCoordinate[2], mesh(vertex).z);
  }
  const double extent = std::max(maxCoordinate[0] - minCoordinate[0],
      std::max(maxCoordinate[1] - minCoordinate[1],
               maxCoordinate[2] - minCoordinate[2]));
  const double tolerance = std::max(1.e-12, extent * 1.e-10);

  typedef std::array<long long, 3> GridKey;
  std::map<GridKey, std::vector<long> > vertexGrid;
  for (long vertex = 0; vertex < mesh.nv; ++vertex) {
    const GridKey key = {{
        static_cast<long long>(std::floor((mesh(vertex).x - minCoordinate[0]) / tolerance)),
        static_cast<long long>(std::floor((mesh(vertex).y - minCoordinate[1]) / tolerance)),
        static_cast<long long>(std::floor((mesh(vertex).z - minCoordinate[2]) / tolerance))}};
    vertexGrid[key].push_back(vertex);
  }

  const double *coordinates = medMesh.getCoords()->getConstPointer();
  std::map<long, long> mappedNodes;
  for (long node = 0; node < medMesh.getNumberOfNodes(); ++node) {
    const Point medPoint = {
        coordinates[3 * node],
        coordinates[3 * node + 1],
        coordinates[3 * node + 2]};
    const GridKey key = {{
        static_cast<long long>(std::floor((medPoint.x - minCoordinate[0]) / tolerance)),
        static_cast<long long>(std::floor((medPoint.y - minCoordinate[1]) / tolerance)),
        static_cast<long long>(std::floor((medPoint.z - minCoordinate[2]) / tolerance))}};
    long mapped = -1;
    double bestDistance = tolerance * tolerance;
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz) {
          const GridKey nearby = {{key[0] + dx, key[1] + dy, key[2] + dz}};
          std::map<GridKey, std::vector<long> >::const_iterator bucket =
              vertexGrid.find(nearby);
          if (bucket == vertexGrid.end())
            continue;
          for (std::vector<long>::const_iterator candidate = bucket->second.begin();
               candidate != bucket->second.end(); ++candidate) {
            const double x = mesh(*candidate).x - medPoint.x;
            const double y = mesh(*candidate).y - medPoint.y;
            const double z = mesh(*candidate).z - medPoint.z;
            const double distance = x * x + y * y + z * z;
            if (distance <= bestDistance) {
              bestDistance = distance;
              mapped = *candidate;
            }
          }
        }
    if (mapped >= 0)
      mappedNodes[node] = mapped;
  }
  if (mappedNodes.empty())
    ExecError(std::string(context) + ": cannot map MED nodes to FreeFEM mesh vertices");
  return mappedNodes;
}

template<class MeshType>
inline std::set<EdgeKey> readMedFeatureEdges(
    const std::string &fileName,
    const std::string &meshName,
    const std::string &requestedCsv,
    const MeshType &mesh,
    const std::set<EdgeKey> &meshEdges,
    int edgeLevel,
    const char *context) {
  using namespace MEDCoupling;
  const std::set<std::string> requested = splitPhysicalNames(requestedCsv);
  if (requested.empty())
    return std::set<EdgeKey>();

  std::string selectedMeshName = meshName;
  if (selectedMeshName.empty()) {
    const std::vector<std::string> names = GetMeshNames(fileName);
    if (names.empty())
      ExecError(std::string(context) + ": MED file contains no mesh");
    if (names.size() > 1)
      ExecError(std::string(context) + ": medMeshName is required for MED files with multiple meshes");
    selectedMeshName = names[0];
  }

  MCAuto<MEDFileUMesh> fileMesh = MEDFileUMesh::New(fileName, selectedMeshName);
  MCAuto<MEDCouplingUMesh> edgeMesh =
      dynamic_cast<MEDCouplingUMesh *>(fileMesh->getMeshAtLevel(edgeLevel));
  if (!edgeMesh)
    ExecError(std::string(context) + ": MED edge level is not an unstructured mesh");
  if (edgeMesh->getMeshDimension() != 1)
    ExecError(std::string(context) + ": requested MED groups must be 1D edge groups");

  const std::vector<std::string> levelGroups =
      fileMesh->getGroupsOnSpecifiedLev(edgeLevel);
  std::set<std::string> foundNames;
  for (std::vector<std::string>::const_iterator group = levelGroups.begin();
       group != levelGroups.end(); ++group)
    if (requested.find(*group) != requested.end())
      foundNames.insert(*group);
  if (foundNames.size() != requested.size()) {
    std::ostringstream message;
    message << context << ": unknown MED edge group(s):";
    for (std::set<std::string>::const_iterator name = requested.begin();
         name != requested.end(); ++name)
      if (foundNames.find(*name) == foundNames.end())
        message << " " << *name;
    ExecError(message.str());
  }

  const std::map<long, long> mappedNodes =
      mapMedNodesToMeshVertices(*edgeMesh, mesh, context);
  std::set<EdgeKey> result;
  for (std::set<std::string>::const_iterator name = requested.begin();
       name != requested.end(); ++name) {
    MCAuto<DataArrayIdType> groupCells =
        fileMesh->getGroupArr(edgeLevel, *name, false);
    if (!groupCells || groupCells->getNumberOfTuples() == 0)
      ExecError(std::string(context) + ": a requested MED edge group contains no SEG2 cells");
    for (mcIdType i = 0; i < groupCells->getNumberOfTuples(); ++i) {
      const mcIdType cell = groupCells->getIJ(i, 0);
      if (cell < 0 || cell >= edgeMesh->getNumberOfCells())
        ExecError(std::string(context) + ": MED edge group refers to an invalid cell");
      if (edgeMesh->getTypeOfCell(cell) != INTERP_KERNEL::NORM_SEG2)
        ExecError(std::string(context) + ": conserveEdge MED groups must contain only SEG2 cells");
      const mcIdType *connectivity =
          edgeMesh->getNodalConnectivity()->getConstPointer();
      const mcIdType *offsets =
          edgeMesh->getNodalConnectivityIndex()->getConstPointer();
      const mcIdType begin = offsets[cell] + 1;
      std::map<long, long>::const_iterator first =
          mappedNodes.find(connectivity[begin]);
      std::map<long, long>::const_iterator second =
          mappedNodes.find(connectivity[begin + 1]);
      if (first == mappedNodes.end() || second == mappedNodes.end())
        ExecError(std::string(context) + ": cannot map a MED edge endpoint to the input mesh");
      const EdgeKey edge = edgeKey(first->second, second->second);
      if (meshEdges.find(edge) == meshEdges.end())
        ExecError(std::string(context) + ": a conserved MED edge is not an edge of the input mesh");
      result.insert(edge);
    }
  }
  return result;
}

} // namespace Pdmt3D

#endif

#endif
