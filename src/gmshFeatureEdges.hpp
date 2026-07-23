/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     Read named physical line groups from an ASCII Gmsh 2.x mesh and map
     their endpoints to a FreeFEM mesh3.

*****************************************************************************/

#ifndef PDMT_GMSH_FEATURE_EDGES_HPP
#define PDMT_GMSH_FEATURE_EDGES_HPP

#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace Pdmt3D {

inline std::string trim(const std::string &value) {
  std::string::size_type begin = 0;
  while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin])))
    ++begin;
  std::string::size_type end = value.size();
  while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1])))
    --end;
  return value.substr(begin, end - begin);
}

inline std::set<std::string> splitPhysicalNames(const std::string &csv) {
  std::set<std::string> names;
  std::string::size_type begin = 0;
  while (begin <= csv.size()) {
    const std::string::size_type comma = csv.find(',', begin);
    const std::string name = trim(csv.substr(
        begin, comma == std::string::npos ? std::string::npos : comma - begin));
    if (!name.empty())
      names.insert(name);
    if (comma == std::string::npos)
      break;
    begin = comma + 1;
  }
  return names;
}

inline bool isGmshLineType(long elementType) {
  // Linear and higher-order line element identifiers used by Gmsh 2.x.
  return elementType == 1 || elementType == 8 || elementType == 26 ||
         elementType == 27 || elementType == 28;
}

struct GmshLine {
  long first;
  long second;
  long physicalTag;
};

template<class MeshType>
inline std::set<EdgeKey> readGmshFeatureEdges(
    const std::string &fileName,
    const std::string &requestedCsv,
    const MeshType &mesh,
    const std::set<EdgeKey> &meshEdges) {
  const std::set<std::string> requested = splitPhysicalNames(requestedCsv);
  if (requested.empty())
    return std::set<EdgeKey>();

  std::ifstream input(fileName.c_str());
  if (!input)
    ExecError("PdmtBuildDual3D: cannot open the Gmsh file used by conserveEdge");

  std::map<long, Point> gmshNodes;
  std::map<long, std::string> lineGroupNames;
  std::vector<GmshLine> lines;
  bool sawMeshFormat = false;
  std::string section;
  while (std::getline(input, section)) {
    section = trim(section);
    if (section == "$MeshFormat") {
      std::string formatLine;
      if (!std::getline(input, formatLine))
        ExecError("PdmtBuildDual3D: invalid Gmsh MeshFormat section");
      std::istringstream format(formatLine);
      double version = 0.0;
      long binary = 0;
      long dataSize = 0;
      format >> version >> binary >> dataSize;
      if (!format || binary != 0 || version >= 4.0)
        ExecError("PdmtBuildDual3D: conserveEdge currently requires an ASCII Gmsh 2.x file");
      sawMeshFormat = true;
    } else if (section == "$PhysicalNames") {
      std::string countLine;
      std::getline(input, countLine);
      const long count = std::strtol(countLine.c_str(), 0, 10);
      for (long i = 0; i < count; ++i) {
        std::string line;
        std::getline(input, line);
        std::istringstream values(line);
        long dimension = -1;
        long tag = -1;
        values >> dimension >> tag;
        const std::string::size_type firstQuote = line.find('"');
        const std::string::size_type lastQuote = line.rfind('"');
        if (dimension == 1 && firstQuote != std::string::npos && lastQuote > firstQuote)
          lineGroupNames[tag] = line.substr(firstQuote + 1, lastQuote - firstQuote - 1);
      }
    } else if (section == "$Nodes") {
      std::string countLine;
      std::getline(input, countLine);
      const long count = std::strtol(countLine.c_str(), 0, 10);
      for (long i = 0; i < count; ++i) {
        std::string line;
        std::getline(input, line);
        std::istringstream values(line);
        long id;
        Point point;
        values >> id >> point.x >> point.y >> point.z;
        if (!values)
          ExecError("PdmtBuildDual3D: invalid node in Gmsh file");
        gmshNodes[id] = point;
      }
    } else if (section == "$Elements") {
      std::string countLine;
      std::getline(input, countLine);
      const long count = std::strtol(countLine.c_str(), 0, 10);
      for (long i = 0; i < count; ++i) {
        std::string line;
        std::getline(input, line);
        std::istringstream values(line);
        long id, elementType, numberOfTags;
        values >> id >> elementType >> numberOfTags;
        if (!values || numberOfTags < 0)
          ExecError("PdmtBuildDual3D: invalid element in Gmsh file");
        long physicalTag = 0;
        for (long tag = 0; tag < numberOfTags; ++tag) {
          long value;
          values >> value;
          if (tag == 0)
            physicalTag = value;
        }
        if (isGmshLineType(elementType)) {
          GmshLine edge;
          edge.physicalTag = physicalTag;
          values >> edge.first >> edge.second;
          if (!values)
            ExecError("PdmtBuildDual3D: invalid line element in Gmsh file");
          lines.push_back(edge);
        }
      }
    }
  }
  if (!sawMeshFormat)
    ExecError("PdmtBuildDual3D: conserveEdge input is not a Gmsh file");

  std::set<long> requestedTags;
  std::set<std::string> foundNames;
  for (std::map<long, std::string>::const_iterator group = lineGroupNames.begin();
       group != lineGroupNames.end(); ++group)
    if (requested.find(group->second) != requested.end()) {
      requestedTags.insert(group->first);
      foundNames.insert(group->second);
    }
  if (foundNames.size() != requested.size()) {
    std::ostringstream message;
    message << "PdmtBuildDual3D: unknown physical edge group(s):";
    for (std::set<std::string>::const_iterator name = requested.begin();
         name != requested.end(); ++name)
      if (foundNames.find(*name) == foundNames.end())
        message << " " << *name;
    ExecError(message.str());
  }

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

  std::map<long, long> mappedNodes;
  std::set<EdgeKey> result;
  std::map<long, long> selectedElementCount;
  for (std::vector<GmshLine>::const_iterator line = lines.begin(); line != lines.end(); ++line) {
    if (requestedTags.find(line->physicalTag) == requestedTags.end())
      continue;
    ++selectedElementCount[line->physicalTag];
    long endpoints[2] = {line->first, line->second};
    long mapped[2] = {-1, -1};
    for (int endpoint = 0; endpoint < 2; ++endpoint) {
      std::map<long, long>::const_iterator known = mappedNodes.find(endpoints[endpoint]);
      if (known != mappedNodes.end()) {
        mapped[endpoint] = known->second;
        continue;
      }
      std::map<long, Point>::const_iterator gmshPoint = gmshNodes.find(endpoints[endpoint]);
      if (gmshPoint == gmshNodes.end())
        ExecError("PdmtBuildDual3D: a conserved Gmsh edge refers to an unknown node");
      const GridKey key = {{
          static_cast<long long>(std::floor((gmshPoint->second.x - minCoordinate[0]) / tolerance)),
          static_cast<long long>(std::floor((gmshPoint->second.y - minCoordinate[1]) / tolerance)),
          static_cast<long long>(std::floor((gmshPoint->second.z - minCoordinate[2]) / tolerance))}};
      double bestDistance = tolerance * tolerance;
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
          for (int dz = -1; dz <= 1; ++dz) {
            const GridKey nearby = {{key[0] + dx, key[1] + dy, key[2] + dz}};
            std::map<GridKey, std::vector<long> >::const_iterator bucket = vertexGrid.find(nearby);
            if (bucket == vertexGrid.end())
              continue;
            for (std::vector<long>::const_iterator candidate = bucket->second.begin();
                 candidate != bucket->second.end(); ++candidate) {
              const double x = mesh(*candidate).x - gmshPoint->second.x;
              const double y = mesh(*candidate).y - gmshPoint->second.y;
              const double z = mesh(*candidate).z - gmshPoint->second.z;
              const double distance = x * x + y * y + z * z;
              if (distance <= bestDistance) {
                bestDistance = distance;
                mapped[endpoint] = *candidate;
              }
            }
          }
      if (mapped[endpoint] < 0)
        ExecError("PdmtBuildDual3D: cannot map a conserved Gmsh edge node to mesh3");
      mappedNodes[endpoints[endpoint]] = mapped[endpoint];
    }
    const EdgeKey edge = edgeKey(mapped[0], mapped[1]);
    if (meshEdges.find(edge) == meshEdges.end())
      ExecError("PdmtBuildDual3D: a conserved physical line is not an edge of the tetrahedral mesh");
    result.insert(edge);
  }
  for (std::set<long>::const_iterator tag = requestedTags.begin(); tag != requestedTags.end(); ++tag)
    if (!selectedElementCount[*tag])
      ExecError("PdmtBuildDual3D: a requested physical edge group contains no line elements");
  return result;
}

} // namespace Pdmt3D

#endif
