/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     Add conforming quadrilateral viscous layers to selected named Gmsh or
     MED boundaries of a two-dimensional dual mesh.

*****************************************************************************/

#ifndef PDMT_VISCOUS_LAYERS_2D_HPP
#define PDMT_VISCOUS_LAYERS_2D_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "gmshFeatureEdges.hpp"

#ifdef MEDCOUPLING
#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"
#endif

namespace Pdmt2D {

struct Point {
  double x;
  double y;
};

typedef std::pair<long, long> EdgeKey;

inline EdgeKey edgeKey(long first, long second) {
  return first < second ? EdgeKey(first, second) : EdgeKey(second, first);
}

inline double signedArea(const std::vector<long> &polygon,
                         const std::vector<Point> &points) {
  double twiceArea = 0.0;
  for (std::size_t i = 0; i < polygon.size(); ++i) {
    const Point &first = points[polygon[i]];
    const Point &second = points[polygon[(i + 1) % polygon.size()]];
    twiceArea += first.x * second.y - second.x * first.y;
  }
  return 0.5 * twiceArea;
}

inline double orientation(const Point &first, const Point &second,
                          const Point &third) {
  return (second.x - first.x) * (third.y - first.y) -
         (second.y - first.y) * (third.x - first.x);
}

inline bool pointOnSegment(const Point &point, const Point &first,
                           const Point &second, double lengthTolerance,
                           double areaTolerance) {
  if (std::fabs(orientation(first, second, point)) > areaTolerance)
    return false;
  return point.x >= std::min(first.x, second.x) - lengthTolerance &&
         point.x <= std::max(first.x, second.x) + lengthTolerance &&
         point.y >= std::min(first.y, second.y) - lengthTolerance &&
         point.y <= std::max(first.y, second.y) + lengthTolerance;
}

inline bool segmentsIntersect(const Point &a, const Point &b,
                              const Point &c, const Point &d,
                              double lengthTolerance,
                              double areaTolerance) {
  const double abc = orientation(a, b, c);
  const double abd = orientation(a, b, d);
  const double cda = orientation(c, d, a);
  const double cdb = orientation(c, d, b);
  if (((abc > areaTolerance && abd < -areaTolerance) ||
       (abc < -areaTolerance && abd > areaTolerance)) &&
      ((cda > areaTolerance && cdb < -areaTolerance) ||
       (cda < -areaTolerance && cdb > areaTolerance)))
    return true;
  return (std::fabs(abc) <= areaTolerance &&
          pointOnSegment(c, a, b, lengthTolerance, areaTolerance)) ||
         (std::fabs(abd) <= areaTolerance &&
          pointOnSegment(d, a, b, lengthTolerance, areaTolerance)) ||
         (std::fabs(cda) <= areaTolerance &&
          pointOnSegment(a, c, d, lengthTolerance, areaTolerance)) ||
         (std::fabs(cdb) <= areaTolerance &&
          pointOnSegment(b, c, d, lengthTolerance, areaTolerance));
}

inline Point lineIntersection(const Point &a, const Point &b,
                              const Point &c, const Point &d,
                              double areaTolerance) {
  const double firstX = b.x - a.x;
  const double firstY = b.y - a.y;
  const double secondX = d.x - c.x;
  const double secondY = d.y - c.y;
  const double denominator =
      firstX * secondY - firstY * secondX;
  if (std::fabs(denominator) <= areaTolerance)
    ExecError("PdmtAddViscousLayers2D: cannot resolve a parallel layer-front intersection");
  const double parameter =
      ((c.x - a.x) * secondY - (c.y - a.y) * secondX) /
      denominator;
  return Point{a.x + parameter * firstX, a.y + parameter * firstY};
}

inline bool simplePolygon(const std::vector<long> &polygon,
                          const std::vector<Point> &points,
                          double lengthTolerance,
                          double areaTolerance) {
  if (polygon.size() < 3)
    return false;
  std::set<long> uniqueNodes(polygon.begin(), polygon.end());
  if (uniqueNodes.size() != polygon.size())
    return false;
  for (std::size_t first = 0; first < polygon.size(); ++first) {
    const std::size_t firstNext = (first + 1) % polygon.size();
    const Point &a = points[polygon[first]];
    const Point &b = points[polygon[firstNext]];
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    if (dx * dx + dy * dy <= lengthTolerance * lengthTolerance)
      return false;
    for (std::size_t second = first + 1; second < polygon.size(); ++second) {
      const std::size_t secondNext = (second + 1) % polygon.size();
      if (first == second || firstNext == second || secondNext == first)
        continue;
      // The first and last edges are adjacent as well.
      if (first == 0 && secondNext == 0)
        continue;
      if (segmentsIntersect(a, b, points[polygon[second]],
                            points[polygon[secondNext]], lengthTolerance,
                            areaTolerance))
        return false;
    }
  }
  return true;
}

inline bool firstSelfIntersection(const std::vector<long> &polygon,
                                  const std::vector<Point> &points,
                                  double lengthTolerance,
                                  double areaTolerance,
                                  std::size_t &firstEdge,
                                  std::size_t &secondEdge) {
  for (std::size_t first = 0; first < polygon.size(); ++first) {
    const std::size_t firstNext = (first + 1) % polygon.size();
    for (std::size_t second = first + 1; second < polygon.size(); ++second) {
      const std::size_t secondNext = (second + 1) % polygon.size();
      if (firstNext == second || secondNext == first ||
          (first == 0 && secondNext == 0))
        continue;
      if (segmentsIntersect(points[polygon[first]],
                            points[polygon[firstNext]],
                            points[polygon[second]],
                            points[polygon[secondNext]], lengthTolerance,
                            areaTolerance)) {
        firstEdge = first;
        secondEdge = second;
        return true;
      }
    }
  }
  return false;
}

inline std::vector<long> mergePolygonsAcrossEdge(
    const std::vector<long> &first, const std::vector<long> &second,
    const EdgeKey &shared) {
  std::size_t firstLocal = first.size();
  for (std::size_t i = 0; i < first.size(); ++i)
    if (edgeKey(first[i], first[(i + 1) % first.size()]) == shared) {
      firstLocal = i;
      break;
    }
  if (firstLocal == first.size())
    ExecError("PdmtAddViscousLayers2D: cannot find the shared core edge");
  std::vector<long> reversedSecond;
  const std::vector<long> *orientedSecond = &second;
  std::size_t secondLocal = second.size();
  for (std::size_t i = 0; i < second.size(); ++i)
    if (second[i] == first[(firstLocal + 1) % first.size()] &&
        second[(i + 1) % second.size()] == first[firstLocal]) {
      secondLocal = i;
      break;
    }
  if (secondLocal == second.size()) {
    reversedSecond.assign(second.rbegin(), second.rend());
    orientedSecond = &reversedSecond;
    for (std::size_t i = 0; i < orientedSecond->size(); ++i)
      if ((*orientedSecond)[i] ==
              first[(firstLocal + 1) % first.size()] &&
          (*orientedSecond)[(i + 1) % orientedSecond->size()] ==
              first[firstLocal]) {
        secondLocal = i;
        break;
      }
  }
  if (secondLocal == orientedSecond->size())
    ExecError("PdmtAddViscousLayers2D: cannot merge polygons across a shared edge");

  // Follow the non-shared path of the first polygon from the shared edge's
  // second endpoint back to its first endpoint, then the corresponding
  // non-shared path through the second polygon.
  std::vector<long> merged;
  for (std::size_t step = 0; step < first.size(); ++step) {
    const long node = first[(firstLocal + 1 + step) % first.size()];
    if (merged.empty() || merged.back() != node)
      merged.push_back(node);
  }
  for (std::size_t step = 2; step < orientedSecond->size(); ++step) {
    const long node =
        (*orientedSecond)[(secondLocal + step) % orientedSecond->size()];
    if (merged.empty() || merged.back() != node)
      merged.push_back(node);
  }
  if (merged.size() > 1 && merged.front() == merged.back())
    merged.pop_back();

  // Two regions that already touch along a second edge can leave a
  // two-vertex backtracking spur in the concatenated walk. It encloses no
  // area and is not part of the union boundary, so remove it before further
  // intersection processing.
  bool removedSpur = true;
  while (removedSpur) {
    removedSpur = false;
    for (std::size_t firstOccurrence = 0;
         firstOccurrence < merged.size() && !removedSpur;
         ++firstOccurrence)
      for (std::size_t secondOccurrence = firstOccurrence + 1;
           secondOccurrence < merged.size(); ++secondOccurrence) {
        if (merged[firstOccurrence] != merged[secondOccurrence])
          continue;
        const std::size_t firstLoopSize =
            secondOccurrence - firstOccurrence;
        const std::size_t secondLoopSize =
            merged.size() - firstLoopSize;
        std::vector<long> kept;
        if (firstLoopSize < 3) {
          for (std::size_t position = secondOccurrence;
               position < merged.size(); ++position)
            kept.push_back(merged[position]);
          for (std::size_t position = 0;
               position < firstOccurrence; ++position)
            kept.push_back(merged[position]);
        } else if (secondLoopSize < 3) {
          for (std::size_t position = firstOccurrence;
               position < secondOccurrence; ++position)
            kept.push_back(merged[position]);
        } else {
          continue;
        }
        merged.swap(kept);
        removedSpur = true;
        break;
      }
  }
  return merged;
}

inline bool convexPolygon(const std::vector<long> &polygon,
                          const std::vector<Point> &points,
                          double areaTolerance) {
  if (polygon.size() < 3)
    return false;
  const double area = signedArea(polygon, points);
  if (std::fabs(area) <= areaTolerance)
    return false;
  const double sign = area > 0.0 ? 1.0 : -1.0;
  for (std::size_t vertex = 0; vertex < polygon.size(); ++vertex) {
    const Point &previous =
        points[polygon[(vertex + polygon.size() - 1) % polygon.size()]];
    const Point &current = points[polygon[vertex]];
    const Point &next = points[polygon[(vertex + 1) % polygon.size()]];
    if (sign * orientation(previous, current, next) < -areaTolerance)
      return false;
  }
  return true;
}

inline bool pointInOrOnTriangle(const Point &point, const Point &first,
                                const Point &second, const Point &third,
                                double sign, double areaTolerance) {
  return sign * orientation(first, second, point) >= -areaTolerance &&
         sign * orientation(second, third, point) >= -areaTolerance &&
         sign * orientation(third, first, point) >= -areaTolerance;
}

// Ear clipping is used only as an intermediate representation. Adjacent ears
// are immediately merged again whenever their union is convex, yielding a
// compact polygon partition instead of exposing a triangle fan.
inline std::vector<std::vector<long> > convexPartition(
    const std::vector<long> &polygon, const std::vector<Point> &points,
    double lengthTolerance, double areaTolerance) {
  std::vector<std::vector<long> > pieces;
  if (convexPolygon(polygon, points, areaTolerance)) {
    pieces.push_back(polygon);
    return pieces;
  }
  const double polygonArea = signedArea(polygon, points);
  const double sign = polygonArea > 0.0 ? 1.0 : -1.0;
  std::vector<long> remaining = polygon;
  while (remaining.size() > 3) {
    bool clipped = false;
    for (std::size_t vertex = 0; vertex < remaining.size(); ++vertex) {
      const long previous =
          remaining[(vertex + remaining.size() - 1) % remaining.size()];
      const long current = remaining[vertex];
      const long next = remaining[(vertex + 1) % remaining.size()];
      if (sign * orientation(points[previous], points[current], points[next]) <=
          areaTolerance)
        continue;
      bool containsVertex = false;
      for (std::size_t candidate = 0; candidate < remaining.size();
           ++candidate) {
        const long node = remaining[candidate];
        if (node == previous || node == current || node == next)
          continue;
        if (pointInOrOnTriangle(points[node], points[previous],
                               points[current], points[next], sign,
                               areaTolerance)) {
          containsVertex = true;
          break;
        }
      }
      if (containsVertex)
        continue;
      std::vector<long> ear(3);
      ear[0] = previous;
      ear[1] = current;
      ear[2] = next;
      pieces.push_back(ear);
      remaining.erase(remaining.begin() + vertex);
      clipped = true;
      break;
    }
    if (!clipped)
      ExecError("PdmtAddViscousLayers2D: cannot partition an absorbed non-convex core patch");
  }
  pieces.push_back(remaining);

  bool mergedPiece = true;
  while (mergedPiece) {
    mergedPiece = false;
    std::map<EdgeKey, std::vector<long> > owners;
    for (std::size_t piece = 0; piece < pieces.size(); ++piece)
      for (std::size_t local = 0; local < pieces[piece].size(); ++local)
        owners[edgeKey(pieces[piece][local],
                       pieces[piece][(local + 1) % pieces[piece].size()])]
            .push_back(static_cast<long>(piece));
    for (std::map<EdgeKey, std::vector<long> >::const_iterator edge =
             owners.begin();
         edge != owners.end() && !mergedPiece; ++edge) {
      if (edge->second.size() != 2)
        continue;
      const long first = edge->second[0];
      const long second = edge->second[1];
      const std::vector<long> merged = mergePolygonsAcrossEdge(
          pieces[first], pieces[second], edge->first);
      if (!simplePolygon(merged, points, lengthTolerance, areaTolerance) ||
          !convexPolygon(merged, points, areaTolerance))
        continue;
      pieces[first] = merged;
      pieces.erase(pieces.begin() + second);
      mergedPiece = true;
    }
  }
  return pieces;
}

inline std::set<long> readPhysicalBoundaryTags(
    const std::string &fileName, const std::string &requestedCsv,
    const Fem2D::Mesh &mesh) {
  using Pdmt3D::allPhysicalEdgeGroupsRequested;
  using Pdmt3D::splitPhysicalNames;
  using Pdmt3D::trim;

  const std::set<std::string> requested = splitPhysicalNames(requestedCsv);
  if (requested.empty())
    ExecError("PdmtAddViscousLayers2D: viscousLayerGroups must not be empty");

  std::ifstream input(fileName.c_str());
  if (!input)
    ExecError("PdmtAddViscousLayers2D: cannot open the Gmsh input file");

  bool sawMeshFormat = false;
  std::map<long, std::string> boundaryNames;
  std::string section;
  while (std::getline(input, section)) {
    section = trim(section);
    if (section == "$MeshFormat") {
      std::string formatLine;
      if (!std::getline(input, formatLine))
        ExecError("PdmtAddViscousLayers2D: invalid Gmsh MeshFormat section");
      std::istringstream format(formatLine);
      double version = 0.0;
      long binary = 0;
      long dataSize = 0;
      format >> version >> binary >> dataSize;
      if (!format || binary != 0 || version < 2.0 || version >= 5.0)
        ExecError("PdmtAddViscousLayers2D: viscous layers require an ASCII Gmsh 2.x or 4.x file");
      sawMeshFormat = true;
    } else if (section == "$PhysicalNames") {
      std::string countLine;
      if (!std::getline(input, countLine))
        ExecError("PdmtAddViscousLayers2D: invalid Gmsh PhysicalNames section");
      const long count = std::strtol(countLine.c_str(), 0, 10);
      for (long i = 0; i < count; ++i) {
        std::string line;
        if (!std::getline(input, line))
          ExecError("PdmtAddViscousLayers2D: truncated Gmsh PhysicalNames section");
        std::istringstream values(line);
        long dimension = -1;
        long tag = -1;
        values >> dimension >> tag;
        const std::string::size_type firstQuote = line.find('"');
        const std::string::size_type lastQuote = line.rfind('"');
        if (dimension == 1 && firstQuote != std::string::npos &&
            lastQuote > firstQuote)
          boundaryNames[tag] =
              line.substr(firstQuote + 1, lastQuote - firstQuote - 1);
      }
    }
  }
  if (!sawMeshFormat)
    ExecError("PdmtAddViscousLayers2D: input is not a Gmsh file");

  std::set<long> loadedBoundaryTags;
  for (long edge = 0; edge < mesh.neb; ++edge)
    loadedBoundaryTags.insert(mesh.be(edge).lab);

  std::set<long> selectedTags;
  std::set<std::string> foundNames;
  for (std::map<long, std::string>::const_iterator group =
           boundaryNames.begin();
       group != boundaryNames.end(); ++group) {
    if (loadedBoundaryTags.find(group->first) == loadedBoundaryTags.end())
      continue;
    if (allPhysicalEdgeGroupsRequested(requested) ||
        requested.find(group->second) != requested.end()) {
      selectedTags.insert(group->first);
      foundNames.insert(group->second);
    }
  }

  if (allPhysicalEdgeGroupsRequested(requested)) {
    if (selectedTags.empty())
      ExecError("PdmtAddViscousLayers2D: Gmsh file contains no populated physical boundary groups");
  } else if (foundNames.size() != requested.size()) {
    std::ostringstream message;
    message << "PdmtAddViscousLayers2D: unknown or empty physical boundary group(s):";
    for (std::set<std::string>::const_iterator name = requested.begin();
         name != requested.end(); ++name)
      if (foundNames.find(*name) == foundNames.end())
        message << " " << *name;
    ExecError(message.str());
  }
  return selectedTags;
}

#ifdef MEDCOUPLING
inline std::set<long> readMedBoundaryTags(
    const std::string &fileName, const std::string &meshName,
    const std::string &requestedCsv, const Fem2D::Mesh &mesh) {
  using Pdmt3D::allPhysicalEdgeGroupsRequested;
  using Pdmt3D::splitPhysicalNames;
  using namespace MEDCoupling;

  const std::set<std::string> requested = splitPhysicalNames(requestedCsv);
  if (requested.empty())
    ExecError("PdmtAddViscousLayers2D: viscousLayerGroups must not be empty");

  std::string selectedMeshName = meshName;
  if (selectedMeshName.empty()) {
    const std::vector<std::string> names = GetMeshNames(fileName);
    if (names.empty())
      ExecError("PdmtAddViscousLayers2D: MED file contains no mesh");
    if (names.size() > 1)
      ExecError("PdmtAddViscousLayers2D: medMeshName is required for MED files with multiple meshes");
    selectedMeshName = names[0];
  }

  MCAuto<MEDFileUMesh> fileMesh =
      MEDFileUMesh::New(fileName, selectedMeshName);
  const std::vector<std::string> levelGroups =
      fileMesh->getGroupsOnSpecifiedLev(-1);
  std::set<std::string> selectedNames;
  if (allPhysicalEdgeGroupsRequested(requested)) {
    selectedNames.insert(levelGroups.begin(), levelGroups.end());
    if (selectedNames.empty())
      ExecError("PdmtAddViscousLayers2D: MED boundary level contains no groups");
  } else {
    for (std::vector<std::string>::const_iterator group =
             levelGroups.begin();
         group != levelGroups.end(); ++group)
      if (requested.find(*group) != requested.end())
        selectedNames.insert(*group);
  }

  if (!allPhysicalEdgeGroupsRequested(requested) &&
      selectedNames.size() != requested.size()) {
    std::ostringstream message;
    message << "PdmtAddViscousLayers2D: unknown MED boundary group(s):";
    for (std::set<std::string>::const_iterator name = requested.begin();
         name != requested.end(); ++name)
      if (selectedNames.find(*name) == selectedNames.end())
        message << " " << *name;
    ExecError(message.str());
  }

  std::set<long> loadedBoundaryTags;
  for (long edge = 0; edge < mesh.neb; ++edge)
    loadedBoundaryTags.insert(mesh.be(edge).lab);

  std::set<long> selectedTags;
  std::set<std::string> emptyNames;
  for (std::set<std::string>::const_iterator name = selectedNames.begin();
       name != selectedNames.end(); ++name) {
    const std::vector<mcIdType> familyIds =
        fileMesh->getFamiliesIdsOnGroup(*name);
    bool populated = false;
    for (std::vector<mcIdType>::const_iterator family = familyIds.begin();
         family != familyIds.end(); ++family)
      if (loadedBoundaryTags.find(static_cast<long>(*family)) !=
          loadedBoundaryTags.end()) {
        selectedTags.insert(static_cast<long>(*family));
        populated = true;
      }
    if (!populated)
      emptyNames.insert(*name);
  }

  if (!emptyNames.empty() && !allPhysicalEdgeGroupsRequested(requested)) {
    std::ostringstream message;
    message << "PdmtAddViscousLayers2D: empty MED boundary group(s):";
    for (std::set<std::string>::const_iterator name = emptyNames.begin();
         name != emptyNames.end(); ++name)
      message << " " << *name;
    ExecError(message.str());
  }
  if (selectedTags.empty())
    ExecError("PdmtAddViscousLayers2D: selected MED groups contain no populated boundary families");
  return selectedTags;
}
#endif

struct PolygonEdge {
  long owner;
  long first;
  long second;
};

struct SelectedEdge {
  long owner;
  long first;
  long second;
  long explicitEdge;
};

struct BoundaryLink {
  long neighbour;
  long explicitEdge;
};

inline void appendUnique(std::vector<long> &values, long value) {
  if (values.empty() || values.back() != value)
    values.push_back(value);
}

inline long addViscousLayers(
    const Fem2D::Mesh &mesh, const std::string &meshFile,
    const std::string &medMeshName, const std::string &requestedGroups,
    long layerCount, double totalThickness, KNM<double> &nodes,
    KN<KN<long> > &cells, KN<KN<long> > &edges, KN<long> &labels) {
  if (layerCount <= 0)
    ExecError("PdmtAddViscousLayers2D: viscousLayerCount must be positive");
  if (!(totalThickness > 0.0) ||
      !std::isfinite(totalThickness))
    ExecError("PdmtAddViscousLayers2D: viscousLayerThickness must be a finite positive number");
  if (nodes.M() < 2)
    ExecError("PdmtAddViscousLayers2D: nodes must have two coordinates");
  if (cells.N() <= 0 || edges.N() <= 0)
    ExecError("PdmtAddViscousLayers2D: the dual mesh has no cells or boundary edges");
  if (labels.N() != cells.N() + edges.N())
    ExecError("PdmtAddViscousLayers2D: labels must contain all cell and boundary-edge labels");

  std::set<long> selectedTags;
  if (meshFile.size() >= 4 &&
      meshFile.substr(meshFile.size() - 4) == ".med") {
#ifdef MEDCOUPLING
    selectedTags =
        readMedBoundaryTags(meshFile, medMeshName, requestedGroups, mesh);
#else
    ExecError("PdmtAddViscousLayers2D: this PDMT build has no MED support");
#endif
  } else {
    selectedTags =
        readPhysicalBoundaryTags(meshFile, requestedGroups, mesh);
  }
  const long oldCellCount = cells.N();
  const long oldEdgeCount = edges.N();

  std::vector<Point> pointList(nodes.N());
  double minX = std::numeric_limits<double>::max();
  double minY = std::numeric_limits<double>::max();
  double maxX = -std::numeric_limits<double>::max();
  double maxY = -std::numeric_limits<double>::max();
  for (long node = 0; node < nodes.N(); ++node) {
    pointList[node].x = nodes(node, 0L);
    pointList[node].y = nodes(node, 1L);
    minX = std::min(minX, pointList[node].x);
    minY = std::min(minY, pointList[node].y);
    maxX = std::max(maxX, pointList[node].x);
    maxY = std::max(maxY, pointList[node].y);
  }
  const double extent = std::max(maxX - minX, maxY - minY);
  const double lengthTolerance = std::max(1.e-12, extent * 1.e-11);
  const double areaTolerance =
      std::max(1.e-24, extent * extent * 1.e-13);

  std::vector<std::vector<long> > polygonList(oldCellCount);
  std::vector<double> originalArea(oldCellCount, 0.0);
  std::map<EdgeKey, std::vector<PolygonEdge> > polygonEdges;
  for (long cell = 0; cell < oldCellCount; ++cell) {
    polygonList[cell].resize(cells[cell].N());
    for (long local = 0; local < cells[cell].N(); ++local) {
      const long node = cells[cell][local];
      if (node < 0 || node >= nodes.N())
        ExecError("PdmtAddViscousLayers2D: cell connectivity contains an invalid node");
      polygonList[cell][local] = node;
    }
    originalArea[cell] = signedArea(polygonList[cell], pointList);
    if (std::fabs(originalArea[cell]) <= areaTolerance)
      ExecError("PdmtAddViscousLayers2D: input contains a zero-area dual polygon");
    for (std::size_t local = 0; local < polygonList[cell].size(); ++local) {
      PolygonEdge occurrence;
      occurrence.owner = cell;
      occurrence.first = polygonList[cell][local];
      occurrence.second =
          polygonList[cell][(local + 1) % polygonList[cell].size()];
      polygonEdges[edgeKey(occurrence.first, occurrence.second)]
          .push_back(occurrence);
    }
  }

  std::vector<long> oldCellLabels(oldCellCount);
  std::vector<long> oldEdgeLabels(oldEdgeCount);
  for (long cell = 0; cell < oldCellCount; ++cell)
    oldCellLabels[cell] = labels[cell];
  for (long edge = 0; edge < oldEdgeCount; ++edge)
    oldEdgeLabels[edge] = labels[oldCellCount + edge];

  const long oldNodeCount = nodes.N();
  std::vector<SelectedEdge> selectedEdges;
  std::vector<char> selectedExplicitEdge(oldEdgeCount, 0);
  std::vector<std::vector<BoundaryLink> > boundaryLinks(oldNodeCount);
  for (long edge = 0; edge < oldEdgeCount; ++edge) {
    if (edges[edge].N() != 2)
      ExecError("PdmtAddViscousLayers2D: explicit boundary connectivity must contain line segments");
    const long first = edges[edge][0];
    const long second = edges[edge][1];
    if (first < 0 || first >= oldNodeCount ||
        second < 0 || second >= oldNodeCount)
      ExecError("PdmtAddViscousLayers2D: boundary connectivity contains an invalid node");
    BoundaryLink firstLink = {second, edge};
    BoundaryLink secondLink = {first, edge};
    boundaryLinks[first].push_back(firstLink);
    boundaryLinks[second].push_back(secondLink);
    if (selectedTags.find(oldEdgeLabels[edge]) == selectedTags.end())
      continue;
    selectedExplicitEdge[edge] = 1;
    const EdgeKey key = edgeKey(first, second);
    std::map<EdgeKey, std::vector<PolygonEdge> >::const_iterator occurrence =
        polygonEdges.find(key);
    if (occurrence == polygonEdges.end() || occurrence->second.size() != 1)
      ExecError("PdmtAddViscousLayers2D: selected boundary segment is not owned by exactly one polygon");
    SelectedEdge selected;
    selected.owner = occurrence->second[0].owner;
    selected.first = occurrence->second[0].first;
    selected.second = occurrence->second[0].second;
    selected.explicitEdge = edge;
    selectedEdges.push_back(selected);
  }
  if (selectedEdges.empty())
    ExecError("PdmtAddViscousLayers2D: selected groups contain no output boundary segments");

  std::map<long, std::vector<Point> > incidentNormals;
  std::vector<std::set<EdgeKey> > selectedCellEdges(oldCellCount);
  std::vector<std::map<EdgeKey, long> > selectedCellEdgeIndex(oldCellCount);
  std::vector<Point> selectedEdgeNormals(selectedEdges.size());
  std::map<long, std::vector<long> > incidentSelectedEdges;
  for (std::size_t edge = 0; edge < selectedEdges.size(); ++edge) {
    const SelectedEdge &selected = selectedEdges[edge];
    const Point &first = pointList[selected.first];
    const Point &second = pointList[selected.second];
    const double dx = second.x - first.x;
    const double dy = second.y - first.y;
    const double length = std::sqrt(dx * dx + dy * dy);
    if (length <= lengthTolerance)
      ExecError("PdmtAddViscousLayers2D: selected boundary contains a zero-length segment");
    Point normal;
    if (originalArea[selected.owner] > 0.0) {
      normal.x = -dy / length;
      normal.y = dx / length;
    } else {
      normal.x = dy / length;
      normal.y = -dx / length;
    }
    incidentNormals[selected.first].push_back(normal);
    incidentNormals[selected.second].push_back(normal);
    selectedEdgeNormals[edge] = normal;
    incidentSelectedEdges[selected.first].push_back(
        static_cast<long>(edge));
    incidentSelectedEdges[selected.second].push_back(
        static_cast<long>(edge));
    selectedCellEdges[selected.owner].insert(
        edgeKey(selected.first, selected.second));
    selectedCellEdgeIndex[selected.owner][
        edgeKey(selected.first, selected.second)] =
        static_cast<long>(edge);
  }

  // Compute the full-thickness displacement of every selected wall node.
  // At the end of a wall, prefer the direction of the adjacent unselected
  // boundary. For an orthogonal corner this slides the inner endpoint down
  // the side wall and keeps every layer full-width.
  std::map<long, Point> selectedDisplacement;
  std::map<long, long> slidingEndpointEdge;
  std::set<long> internalEndpoint;
  std::set<long> beveledCorner;
  for (std::map<long, std::vector<Point> >::const_iterator nodeNormals =
           incidentNormals.begin();
       nodeNormals != incidentNormals.end(); ++nodeNormals) {
    const long node = nodeNormals->first;
    Point direction = {0.0, 0.0};
    double scale = totalThickness;
    long adjacentUnselectedEdge = -1;
    long adjacentUnselectedNode = -1;
    for (std::size_t link = 0; link < boundaryLinks[node].size(); ++link)
      if (!selectedExplicitEdge[boundaryLinks[node][link].explicitEdge]) {
        adjacentUnselectedEdge =
            boundaryLinks[node][link].explicitEdge;
        adjacentUnselectedNode = boundaryLinks[node][link].neighbour;
      }

    bool slidesOnBoundary = false;
    if (nodeNormals->second.size() == 1 &&
        adjacentUnselectedNode >= 0) {
      direction.x =
          pointList[adjacentUnselectedNode].x - pointList[node].x;
      direction.y =
          pointList[adjacentUnselectedNode].y - pointList[node].y;
      const double length =
          std::sqrt(direction.x * direction.x + direction.y * direction.y);
      if (length > lengthTolerance) {
        direction.x /= length;
        direction.y /= length;
        const double projection =
            direction.x * nodeNormals->second[0].x +
            direction.y * nodeNormals->second[0].y;
        if (projection > 1.e-8) {
          scale = totalThickness / projection;
          slidesOnBoundary = true;
          slidingEndpointEdge[node] = adjacentUnselectedEdge;
        }
      }
    }

    if (!slidesOnBoundary) {
      direction.x = 0.0;
      direction.y = 0.0;
      for (std::size_t i = 0; i < nodeNormals->second.size(); ++i) {
        direction.x += nodeNormals->second[i].x;
        direction.y += nodeNormals->second[i].y;
      }
      const double directionLength =
          std::sqrt(direction.x * direction.x + direction.y * direction.y);
      if (directionLength <= 1.e-12)
        ExecError("PdmtAddViscousLayers2D: cannot determine an inward normal at a selected boundary corner");
      direction.x /= directionLength;
      direction.y /= directionLength;
      double minimumProjection = 1.0;
      for (std::size_t i = 0; i < nodeNormals->second.size(); ++i)
        minimumProjection = std::min(
            minimumProjection,
            direction.x * nodeNormals->second[i].x +
                direction.y * nodeNormals->second[i].y);
      if (minimumProjection <= 1.e-8)
        ExecError("PdmtAddViscousLayers2D: selected boundary corner does not admit an inward layer");
      // The intersection of two offset lines escapes to infinity when a
      // boundary cusp approaches 180 degrees. A large but finite miter can
      // also run past the next boundary vertex when the local edge is short,
      // making two otherwise valid offset segments cross. Replace either
      // case by a bevel fan with an independent normal-offset chain for each
      // incident edge.
      const double requestedMiter = totalThickness / minimumProjection;
      const double maximumMiter = 2.0 * totalThickness;
      bool miterOverrunsEdge = false;
      if (incidentSelectedEdges[node].size() == 2) {
        const double tangentialReach = std::sqrt(std::max(
            0.0, requestedMiter * requestedMiter -
                     totalThickness * totalThickness));
        double minimumIncidentLength =
            std::numeric_limits<double>::max();
        for (std::size_t incident = 0;
             incident < incidentSelectedEdges[node].size(); ++incident) {
          const SelectedEdge &selected =
              selectedEdges[incidentSelectedEdges[node][incident]];
          const Point &first = pointList[selected.first];
          const Point &second = pointList[selected.second];
          const double dx = second.x - first.x;
          const double dy = second.y - first.y;
          minimumIncidentLength = std::min(
              minimumIncidentLength, std::sqrt(dx * dx + dy * dy));
        }
        miterOverrunsEdge =
            tangentialReach >= minimumIncidentLength - lengthTolerance;
      }
      if ((requestedMiter > maximumMiter || miterOverrunsEdge) &&
          incidentSelectedEdges[node].size() == 2) {
        beveledCorner.insert(node);
        // The harmonic core deformation sees the midpoint of the bevel. The
        // two exact edge-normal endpoints are created separately below.
        scale = totalThickness * directionLength /
                static_cast<double>(nodeNormals->second.size());
      } else {
        scale = requestedMiter;
      }
      if (adjacentUnselectedEdge >= 0)
        internalEndpoint.insert(node);
    }
    Point displacement = {scale * direction.x, scale * direction.y};
    selectedDisplacement[node] = displacement;
  }

  // Build the dual-node graph used by the deformation solve.
  std::vector<std::set<long> > nodeNeighbours(oldNodeCount);
  for (long cell = 0; cell < oldCellCount; ++cell)
    for (std::size_t local = 0; local < polygonList[cell].size(); ++local) {
      const long first = polygonList[cell][local];
      const long second =
          polygonList[cell][(local + 1) % polygonList[cell].size()];
      nodeNeighbours[first].insert(second);
      nodeNeighbours[second].insert(first);
    }

  std::vector<Point> displacement(oldNodeCount, Point{0.0, 0.0});
  std::vector<char> boundaryFixed(oldNodeCount, 0);
  std::vector<char> zeroAnchor(oldNodeCount, 0);
  for (std::map<long, Point>::const_iterator selected =
           selectedDisplacement.begin();
       selected != selectedDisplacement.end(); ++selected) {
    displacement[selected->first] = selected->second;
    boundaryFixed[selected->first] = 1;
  }

  // Preserve unselected corners and physical-label junctions as zero
  // displacement anchors. Intermediate side-wall nodes remain free to slide.
  const double straightLimit = -std::cos(20.0 * std::acos(-1.0) / 180.0);
  for (long node = 0; node < oldNodeCount; ++node) {
    if (boundaryLinks[node].size() != 2 ||
        selectedDisplacement.find(node) != selectedDisplacement.end())
      continue;
    const BoundaryLink &firstLink = boundaryLinks[node][0];
    const BoundaryLink &secondLink = boundaryLinks[node][1];
    const Point firstVector = {
        pointList[firstLink.neighbour].x - pointList[node].x,
        pointList[firstLink.neighbour].y - pointList[node].y};
    const Point secondVector = {
        pointList[secondLink.neighbour].x - pointList[node].x,
        pointList[secondLink.neighbour].y - pointList[node].y};
    const double firstLength = std::sqrt(
        firstVector.x * firstVector.x + firstVector.y * firstVector.y);
    const double secondLength = std::sqrt(
        secondVector.x * secondVector.x + secondVector.y * secondVector.y);
    const double cosine =
        (firstVector.x * secondVector.x +
         firstVector.y * secondVector.y) /
        (firstLength * secondLength);
    if (oldEdgeLabels[firstLink.explicitEdge] !=
            oldEdgeLabels[secondLink.explicitEdge] ||
        cosine > straightLimit) {
      boundaryFixed[node] = 1;
      zeroAnchor[node] = 1;
    }
  }

  // Treat boundary components without a selected wall as fixed. On a smooth
  // component with a selected arc but no feature anchor, fix the boundary
  // node farthest from the selected arc.
  std::vector<char> boundaryVisited(oldNodeCount, 0);
  for (long seed = 0; seed < oldNodeCount; ++seed) {
    if (boundaryLinks[seed].empty() || boundaryVisited[seed])
      continue;
    std::vector<long> component;
    std::vector<long> pending(1, seed);
    std::vector<long> componentSelected;
    boundaryVisited[seed] = 1;
    bool hasSelected = false;
    bool hasZeroAnchor = false;
    for (std::size_t position = 0; position < pending.size(); ++position) {
      const long node = pending[position];
      component.push_back(node);
      hasSelected =
          hasSelected ||
          selectedDisplacement.find(node) != selectedDisplacement.end();
      if (selectedDisplacement.find(node) !=
          selectedDisplacement.end())
        componentSelected.push_back(node);
      hasZeroAnchor = hasZeroAnchor || zeroAnchor[node];
      for (std::size_t link = 0; link < boundaryLinks[node].size(); ++link) {
        const long neighbour = boundaryLinks[node][link].neighbour;
        if (!boundaryVisited[neighbour]) {
          boundaryVisited[neighbour] = 1;
          pending.push_back(neighbour);
        }
      }
    }
    if (!hasSelected) {
      for (std::size_t i = 0; i < component.size(); ++i)
        boundaryFixed[component[i]] = 1;
    } else if (!hasZeroAnchor &&
               component.size() > componentSelected.size()) {
      long farthest = -1;
      double farthestDistance = -1.0;
      for (std::size_t i = 0; i < component.size(); ++i) {
        const long candidate = component[i];
        if (selectedDisplacement.find(candidate) !=
            selectedDisplacement.end())
          continue;
        double nearestSelected = std::numeric_limits<double>::max();
        for (std::size_t selected = 0;
             selected < componentSelected.size(); ++selected) {
          const double dx =
              pointList[candidate].x -
              pointList[componentSelected[selected]].x;
          const double dy =
              pointList[candidate].y -
              pointList[componentSelected[selected]].y;
          nearestSelected =
              std::min(nearestSelected, dx * dx + dy * dy);
        }
        if (nearestSelected > farthestDistance) {
          farthestDistance = nearestSelected;
          farthest = candidate;
        }
      }
      if (farthest >= 0)
        boundaryFixed[farthest] = 1;
    }
  }

  const double displacementTolerance =
      std::max(1.e-13, totalThickness * 1.e-12);
  // One-dimensional harmonic extension along unselected boundary chains.
  for (long iteration = 0; iteration < 20000; ++iteration) {
    double maximumChange = 0.0;
    for (long node = 0; node < oldNodeCount; ++node) {
      if (boundaryLinks[node].empty() || boundaryFixed[node])
        continue;
      Point updated = {0.0, 0.0};
      double weightSum = 0.0;
      for (std::size_t link = 0; link < boundaryLinks[node].size(); ++link) {
        const long neighbour = boundaryLinks[node][link].neighbour;
        const double dx = pointList[node].x - pointList[neighbour].x;
        const double dy = pointList[node].y - pointList[neighbour].y;
        const double weight = 1.0 / std::sqrt(dx * dx + dy * dy);
        updated.x += weight * displacement[neighbour].x;
        updated.y += weight * displacement[neighbour].y;
        weightSum += weight;
      }
      updated.x /= weightSum;
      updated.y /= weightSum;
      maximumChange =
          std::max(maximumChange,
                   std::max(std::fabs(updated.x - displacement[node].x),
                            std::fabs(updated.y - displacement[node].y)));
      displacement[node] = updated;
    }
    if (maximumChange < displacementTolerance)
      break;
  }
  for (long node = 0; node < oldNodeCount; ++node)
    if (!boundaryLinks[node].empty())
      boundaryFixed[node] = 1;

  // Harmonic deformation of the interior. Inverse-square edge weights make
  // the coarse dual graph follow affine compression much more closely than
  // an unweighted graph, allowing the layer to push several cell rows.
  for (long iteration = 0; iteration < 30000; ++iteration) {
    double maximumChange = 0.0;
    for (long node = 0; node < oldNodeCount; ++node) {
      if (boundaryFixed[node] || nodeNeighbours[node].empty())
        continue;
      Point updated = {0.0, 0.0};
      double weightSum = 0.0;
      for (std::set<long>::const_iterator neighbour =
               nodeNeighbours[node].begin();
           neighbour != nodeNeighbours[node].end(); ++neighbour) {
        const double dx = pointList[node].x - pointList[*neighbour].x;
        const double dy = pointList[node].y - pointList[*neighbour].y;
        const double lengthSquared = dx * dx + dy * dy;
        const double weight = 1.0 / lengthSquared;
        updated.x += weight * displacement[*neighbour].x;
        updated.y += weight * displacement[*neighbour].y;
        weightSum += weight;
      }
      updated.x /= weightSum;
      updated.y /= weightSum;
      maximumChange =
          std::max(maximumChange,
                   std::max(std::fabs(updated.x - displacement[node].x),
                            std::fabs(updated.y - displacement[node].y)));
      displacement[node] = updated;
    }
    if (maximumChange < displacementTolerance)
      break;
  }

  // Move the reusable core nodes. Selected wall nodes remain on the outer
  // boundary and get separate displaced copies for the core interface.
  for (long node = 0; node < oldNodeCount; ++node)
    if (selectedDisplacement.find(node) == selectedDisplacement.end()) {
      pointList[node].x += displacement[node].x;
      pointList[node].y += displacement[node].y;
    }

  std::map<long, std::vector<long> > layerNodes;
  for (std::map<long, Point>::const_iterator selected =
           selectedDisplacement.begin();
       selected != selectedDisplacement.end(); ++selected) {
    if (beveledCorner.find(selected->first) != beveledCorner.end())
      continue;
    std::vector<long> levels(layerCount + 1);
    levels[0] = selected->first;
    const Point original = pointList[selected->first];
    for (long layer = 1; layer <= layerCount; ++layer) {
      const double fraction =
          static_cast<double>(layer) / static_cast<double>(layerCount);
      const Point point = {
          original.x + fraction * selected->second.x,
          original.y + fraction * selected->second.y};
      levels[layer] = static_cast<long>(pointList.size());
      pointList.push_back(point);
    }
    layerNodes[selected->first] = levels;
  }

  // A beveled corner has one offset chain per incident selected edge. Other
  // selected nodes share their ordinary miter/normal chain between edges.
  typedef std::pair<long, long> SelectedEndpointKey;
  std::map<SelectedEndpointKey, std::vector<long> > endpointLayerNodes;
  for (std::size_t edge = 0; edge < selectedEdges.size(); ++edge) {
    const long endpoints[2] = {selectedEdges[edge].first,
                               selectedEdges[edge].second};
    for (long endpoint = 0; endpoint < 2; ++endpoint) {
      const long node = endpoints[endpoint];
      const SelectedEndpointKey key(static_cast<long>(edge), node);
      if (beveledCorner.find(node) == beveledCorner.end()) {
        endpointLayerNodes[key] = layerNodes[node];
        continue;
      }
      std::vector<long> levels(layerCount + 1);
      levels[0] = node;
      const Point original = pointList[node];
      const Point normal = selectedEdgeNormals[edge];
      for (long layer = 1; layer <= layerCount; ++layer) {
        const double distance =
            totalThickness * static_cast<double>(layer) /
            static_cast<double>(layerCount);
        const Point point = {original.x + distance * normal.x,
                             original.y + distance * normal.y};
        levels[layer] = static_cast<long>(pointList.size());
        pointList.push_back(point);
      }
      endpointLayerNodes[key] = levels;
    }
  }

  // Replace the single straight bevel chord by a circular polygonal cap.
  // The handle order follows the owning polygon through the corner, so the
  // cap sweep is on the fluid side and has the same orientation as the core.
  std::map<long, std::vector<SelectedEndpointKey> > roundedCornerHandles;
  long roundedCornerRayCount = 0;
  long artificialHandle = -1;
  const double maximumCapAngle = std::acos(-1.0) / 12.0;
  for (std::set<long>::const_iterator corner = beveledCorner.begin();
       corner != beveledCorner.end(); ++corner) {
    const std::vector<long> &incident = incidentSelectedEdges[*corner];
    if (incident.size() != 2)
      ExecError("PdmtAddViscousLayers2D: a rounded corner must have two incident selected edges");
    const long owner = selectedEdges[incident[0]].owner;
    if (selectedEdges[incident[1]].owner != owner)
      ExecError("PdmtAddViscousLayers2D: rounded corner segments must belong to one dual polygon");
    const std::vector<long> &ownerPolygon = polygonList[owner];
    std::size_t position = ownerPolygon.size();
    for (std::size_t local = 0; local < ownerPolygon.size(); ++local)
      if (ownerPolygon[local] == *corner) {
        position = local;
        break;
      }
    if (position == ownerPolygon.size())
      ExecError("PdmtAddViscousLayers2D: cannot locate a rounded corner in its owner polygon");
    const EdgeKey previousEdge = edgeKey(
        ownerPolygon[(position + ownerPolygon.size() - 1) %
                     ownerPolygon.size()],
        *corner);
    const EdgeKey nextEdge = edgeKey(
        *corner, ownerPolygon[(position + 1) % ownerPolygon.size()]);
    const long firstEdge = selectedCellEdgeIndex[owner][previousEdge];
    const long secondEdge = selectedCellEdgeIndex[owner][nextEdge];
    const Point &firstNormal = selectedEdgeNormals[firstEdge];
    const Point &secondNormal = selectedEdgeNormals[secondEdge];
    const double firstAngle = std::atan2(firstNormal.y, firstNormal.x);
    const double secondAngle = std::atan2(secondNormal.y, secondNormal.x);
    double sweep = secondAngle - firstAngle;
    const double twoPi = 2.0 * std::acos(-1.0);
    // The normal bisector used by the deformation lies on the shorter arc.
    // Follow that same arc here; using the polygon's global orientation can
    // accidentally select the complementary near-360-degree sweep at a
    // re-entrant boundary vertex.
    while (sweep <= -std::acos(-1.0))
      sweep += twoPi;
    while (sweep > std::acos(-1.0))
      sweep -= twoPi;
    const long sectors = std::max(
        3L, static_cast<long>(std::ceil(std::fabs(sweep) /
                                       maximumCapAngle)));
    std::vector<SelectedEndpointKey> handles;
    handles.push_back(SelectedEndpointKey(firstEdge, *corner));
    const Point origin = pointList[*corner];
    for (long sector = 1; sector < sectors; ++sector) {
      const double fraction = static_cast<double>(sector) /
                              static_cast<double>(sectors);
      const double angle = firstAngle + fraction * sweep;
      const SelectedEndpointKey handle(artificialHandle--, *corner);
      std::vector<long> levels(layerCount + 1);
      levels[0] = *corner;
      for (long layer = 1; layer <= layerCount; ++layer) {
        const double distance = totalThickness *
                                static_cast<double>(layer) /
                                static_cast<double>(layerCount);
        levels[layer] = static_cast<long>(pointList.size());
        pointList.push_back(Point{origin.x + distance * std::cos(angle),
                                  origin.y + distance * std::sin(angle)});
      }
      endpointLayerNodes[handle] = levels;
      handles.push_back(handle);
      ++roundedCornerRayCount;
    }
    handles.push_back(SelectedEndpointKey(secondEdge, *corner));
    roundedCornerHandles[*corner] = handles;
  }

  // Build the cyclic order of the selected front. For a thickness larger
  // than the local radius of curvature, a raw offset contour grows a small
  // self-intersecting loop. At each layer level, trim that loop at its first
  // intersection and identify all of its intermediate front vertices with
  // the intersection point. The adjacent layer cells consequently become a
  // conforming triangle fan instead of overlapping quadrilaterals.
  typedef std::map<SelectedEndpointKey,
                   std::vector<SelectedEndpointKey> > FrontAdjacency;
  FrontAdjacency frontAdjacency;
  for (std::size_t edge = 0; edge < selectedEdges.size(); ++edge) {
    const SelectedEndpointKey first(static_cast<long>(edge),
                                    selectedEdges[edge].first);
    const SelectedEndpointKey second(static_cast<long>(edge),
                                     selectedEdges[edge].second);
    frontAdjacency[first].push_back(second);
    frontAdjacency[second].push_back(first);
  }
  for (std::map<long, std::vector<long> >::const_iterator incident =
           incidentSelectedEdges.begin();
       incident != incidentSelectedEdges.end(); ++incident) {
    if (incident->second.size() != 2)
      continue;
    const std::map<long, std::vector<SelectedEndpointKey> >::const_iterator
        rounded = roundedCornerHandles.find(incident->first);
    if (rounded != roundedCornerHandles.end()) {
      for (std::size_t handle = 0; handle + 1 < rounded->second.size();
           ++handle) {
        frontAdjacency[rounded->second[handle]].push_back(
            rounded->second[handle + 1]);
        frontAdjacency[rounded->second[handle + 1]].push_back(
            rounded->second[handle]);
      }
    } else {
      const SelectedEndpointKey first(incident->second[0], incident->first);
      const SelectedEndpointKey second(incident->second[1], incident->first);
      frontAdjacency[first].push_back(second);
      frontAdjacency[second].push_back(first);
    }
  }

  std::vector<std::vector<SelectedEndpointKey> > closedFrontComponents;
  std::set<SelectedEndpointKey> visitedFrontHandles;
  for (FrontAdjacency::const_iterator seed = frontAdjacency.begin();
       seed != frontAdjacency.end(); ++seed) {
    if (visitedFrontHandles.find(seed->first) !=
            visitedFrontHandles.end() ||
        seed->second.size() != 2)
      continue;
    std::vector<SelectedEndpointKey> component;
    SelectedEndpointKey previous(-1, -1);
    SelectedEndpointKey current = seed->first;
    bool closed = false;
    for (std::size_t step = 0; step <= frontAdjacency.size(); ++step) {
      if (step > 0 && current == seed->first) {
        closed = true;
        break;
      }
      const FrontAdjacency::const_iterator neighbours =
          frontAdjacency.find(current);
      if (neighbours == frontAdjacency.end() ||
          neighbours->second.size() != 2)
        break;
      component.push_back(current);
      visitedFrontHandles.insert(current);
      const SelectedEndpointKey next =
          neighbours->second[0] == previous ? neighbours->second[1]
                                            : neighbours->second[0];
      previous = current;
      current = next;
    }
    if (closed && component.size() >= 3)
      closedFrontComponents.push_back(component);
  }

  long trimmedFrontLoopCount = 0;
  std::vector<std::set<SelectedEndpointKey> > persistentFrontGroups;
  std::vector<Point> persistentFrontAnchors;
  std::vector<Point> persistentFrontDirections;
  for (long level = 1; level <= layerCount; ++level) {
    // Once an offset loop has disappeared, keep the same handle group
    // collapsed at every subsequent level. Re-expanding or shifting that
    // group between levels creates crossed spokes and visually disconnected
    // layer rows at sharp ends.
    for (std::size_t group = 0; group < persistentFrontGroups.size();
         ++group) {
      const double fraction =
          static_cast<double>(level) / static_cast<double>(layerCount);
      Point representative = {
          persistentFrontAnchors[group].x +
              fraction * persistentFrontDirections[group].x,
          persistentFrontAnchors[group].y +
              fraction * persistentFrontDirections[group].y};
      std::vector<std::pair<Point, Point> > boundingLines;
      for (std::set<SelectedEndpointKey>::const_iterator handle =
               persistentFrontGroups[group].begin();
           handle != persistentFrontGroups[group].end(); ++handle) {
        const std::vector<SelectedEndpointKey> &neighbours =
            frontAdjacency[*handle];
        for (std::size_t neighbour = 0; neighbour < neighbours.size();
             ++neighbour)
          if (persistentFrontGroups[group].find(neighbours[neighbour]) ==
              persistentFrontGroups[group].end())
            boundingLines.push_back(std::make_pair(
                pointList[endpointLayerNodes[*handle][level]],
                pointList[endpointLayerNodes[neighbours[neighbour]][level]]));
      }
      if (boundingLines.size() == 2 &&
          std::fabs(orientation(boundingLines[0].first,
                                boundingLines[0].second,
                                boundingLines[1].second) -
                    orientation(boundingLines[0].first,
                                boundingLines[0].second,
                                boundingLines[1].first)) > areaTolerance)
        representative = lineIntersection(
            boundingLines[0].first, boundingLines[0].second,
            boundingLines[1].first, boundingLines[1].second,
            areaTolerance);
      std::set<long> oldNodes;
      for (std::set<SelectedEndpointKey>::const_iterator handle =
               persistentFrontGroups[group].begin();
           handle != persistentFrontGroups[group].end(); ++handle) {
        const long node = endpointLayerNodes[*handle][level];
        oldNodes.insert(node);
      }
      if (oldNodes.empty())
        continue;
      const long representativeNode = static_cast<long>(pointList.size());
      pointList.push_back(representative);
      for (std::set<SelectedEndpointKey>::const_iterator handle =
               persistentFrontGroups[group].begin();
           handle != persistentFrontGroups[group].end(); ++handle)
        endpointLayerNodes[*handle][level] = representativeNode;
      for (std::map<long, std::vector<long> >::iterator node =
               layerNodes.begin(); node != layerNodes.end(); ++node)
        if (oldNodes.find(node->second[level]) != oldNodes.end())
          node->second[level] = representativeNode;
    }
    for (std::size_t componentIndex = 0;
         componentIndex < closedFrontComponents.size(); ++componentIndex) {
      const std::vector<SelectedEndpointKey> &handles =
          closedFrontComponents[componentIndex];
      for (std::size_t repair = 0; repair < handles.size(); ++repair) {
        std::vector<long> contour;
        for (std::size_t handle = 0; handle < handles.size(); ++handle) {
          const long node = endpointLayerNodes[handles[handle]][level];
          if (contour.empty() || contour.back() != node)
            contour.push_back(node);
        }
        if (contour.size() > 1 && contour.front() == contour.back())
          contour.pop_back();
        if (contour.size() < 3)
          ExecError("PdmtAddViscousLayers2D: layer-front trimming collapsed a boundary component");

        std::size_t firstCrossing = 0;
        std::size_t secondCrossing = 0;
        if (!firstSelfIntersection(contour, pointList, lengthTolerance,
                                   areaTolerance, firstCrossing,
                                   secondCrossing))
          break;
        const Point intersection = lineIntersection(
            pointList[contour[firstCrossing]],
            pointList[contour[(firstCrossing + 1) % contour.size()]],
            pointList[contour[secondCrossing]],
            pointList[contour[(secondCrossing + 1) % contour.size()]],
            areaTolerance);

        std::vector<long> firstLoop;
        for (std::size_t position = (firstCrossing + 1) % contour.size();;
             position = (position + 1) % contour.size()) {
          firstLoop.push_back(contour[position]);
          if (position == secondCrossing)
            break;
        }
        std::vector<long> secondLoop;
        for (std::size_t position = (secondCrossing + 1) % contour.size();;
             position = (position + 1) % contour.size()) {
          secondLoop.push_back(contour[position]);
          if (position == firstCrossing)
            break;
        }
        const std::vector<long> *loops[2] = {&firstLoop, &secondLoop};
        double loopAreas[2] = {0.0, 0.0};
        for (long loop = 0; loop < 2; ++loop) {
          Point previousPoint = intersection;
          for (std::size_t vertex = 0; vertex < loops[loop]->size();
               ++vertex) {
            const Point &currentPoint =
                pointList[(*loops[loop])[vertex]];
            loopAreas[loop] += previousPoint.x * currentPoint.y -
                               currentPoint.x * previousPoint.y;
            previousPoint = currentPoint;
          }
          loopAreas[loop] += previousPoint.x * intersection.y -
                             intersection.x * previousPoint.y;
          loopAreas[loop] = std::fabs(0.5 * loopAreas[loop]);
        }
        const std::vector<long> &trimmed =
            loopAreas[0] <= loopAreas[1] ? firstLoop : secondLoop;
        const long intersectionNode = static_cast<long>(pointList.size());
        pointList.push_back(intersection);
        std::set<long> collapsedNodes(trimmed.begin(), trimmed.end());
        std::set<SelectedEndpointKey> collapsedHandles;
        for (std::map<SelectedEndpointKey, std::vector<long> >::const_iterator
                 endpoint = endpointLayerNodes.begin();
             endpoint != endpointLayerNodes.end(); ++endpoint)
          if (collapsedNodes.find(endpoint->second[level]) !=
              collapsedNodes.end())
            collapsedHandles.insert(endpoint->first);
        for (std::map<SelectedEndpointKey, std::vector<long> >::iterator
                 endpoint = endpointLayerNodes.begin();
             endpoint != endpointLayerNodes.end(); ++endpoint)
          if (collapsedNodes.find(endpoint->second[level]) !=
              collapsedNodes.end())
            endpoint->second[level] = intersectionNode;
        for (std::map<long, std::vector<long> >::iterator node =
                 layerNodes.begin();
             node != layerNodes.end(); ++node)
          if (collapsedNodes.find(node->second[level]) !=
              collapsedNodes.end())
            node->second[level] = intersectionNode;
        // Merge overlapping persistent groups so each handle belongs to one
        // nested collapse only.
        for (std::size_t group = 0; group < persistentFrontGroups.size();) {
          bool overlaps = false;
          for (std::set<SelectedEndpointKey>::const_iterator handle =
                   collapsedHandles.begin();
               handle != collapsedHandles.end(); ++handle)
            if (persistentFrontGroups[group].find(*handle) !=
                persistentFrontGroups[group].end()) {
              overlaps = true;
              break;
            }
          if (overlaps) {
            collapsedHandles.insert(persistentFrontGroups[group].begin(),
                                    persistentFrontGroups[group].end());
            persistentFrontGroups.erase(
                persistentFrontGroups.begin() + group);
            persistentFrontAnchors.erase(
                persistentFrontAnchors.begin() + group);
            persistentFrontDirections.erase(
                persistentFrontDirections.begin() + group);
          } else {
            ++group;
          }
        }
        if (!collapsedHandles.empty()) {
          Point anchor = {0.0, 0.0};
          std::set<long> wallNodes;
          for (std::set<SelectedEndpointKey>::const_iterator handle =
                   collapsedHandles.begin();
               handle != collapsedHandles.end(); ++handle) {
            const long wallNode = endpointLayerNodes[*handle][0];
            if (wallNodes.insert(wallNode).second) {
              anchor.x += pointList[wallNode].x;
              anchor.y += pointList[wallNode].y;
            }
          }
          anchor.x /= static_cast<double>(wallNodes.size());
          anchor.y /= static_cast<double>(wallNodes.size());
          const double fraction =
              static_cast<double>(level) / static_cast<double>(layerCount);
          const Point direction = {
              (intersection.x - anchor.x) / fraction,
              (intersection.y - anchor.y) / fraction};
          persistentFrontGroups.push_back(collapsedHandles);
          persistentFrontAnchors.push_back(anchor);
          persistentFrontDirections.push_back(direction);
        }
        ++trimmedFrontLoopCount;
      }
    }
  }

  // Persistent collapse groups can exchange geometric order between two
  // coarse levels even though each individual front is simple. Untwist that
  // correspondence by swapping the complete next-level equivalence classes;
  // this preserves both front contours and reconnects the layer rows in their
  // cyclic order at sharp ends.
  long untwistedLayerConnectionCount = 0;
  for (long layer = 0; layer < layerCount; ++layer) {
    for (std::size_t pass = 0; pass < selectedEdges.size(); ++pass) {
      bool changed = false;
      for (std::size_t edge = 0; edge < selectedEdges.size(); ++edge) {
        const SelectedEdge &selected = selectedEdges[edge];
        const SelectedEndpointKey firstKey(static_cast<long>(edge),
                                           selected.first);
        const SelectedEndpointKey secondKey(static_cast<long>(edge),
                                            selected.second);
        const long firstCurrent = endpointLayerNodes[firstKey][layer];
        const long secondCurrent = endpointLayerNodes[secondKey][layer];
        const long firstNext = endpointLayerNodes[firstKey][layer + 1];
        const long secondNext = endpointLayerNodes[secondKey][layer + 1];
        if (firstNext == secondNext)
          continue;
        std::vector<long> currentQuad(4);
        currentQuad[0] = firstCurrent;
        currentQuad[1] = secondCurrent;
        currentQuad[2] = secondNext;
        currentQuad[3] = firstNext;
        const double currentArea = signedArea(currentQuad, pointList);
        if (simplePolygon(currentQuad, pointList, lengthTolerance,
                          areaTolerance) &&
            std::fabs(currentArea) > areaTolerance &&
            currentArea * originalArea[selected.owner] > 0.0)
          continue;
        std::vector<long> swappedQuad(4);
        swappedQuad[0] = firstCurrent;
        swappedQuad[1] = secondCurrent;
        swappedQuad[2] = firstNext;
        swappedQuad[3] = secondNext;
        const double swappedArea = signedArea(swappedQuad, pointList);
        if (!simplePolygon(swappedQuad, pointList, lengthTolerance,
                           areaTolerance) ||
            std::fabs(swappedArea) <= areaTolerance ||
            swappedArea * originalArea[selected.owner] <= 0.0)
          continue;
        for (std::map<SelectedEndpointKey, std::vector<long> >::iterator
                 endpoint = endpointLayerNodes.begin();
             endpoint != endpointLayerNodes.end(); ++endpoint) {
          long &node = endpoint->second[layer + 1];
          if (node == firstNext)
            node = secondNext;
          else if (node == secondNext)
            node = firstNext;
        }
        for (std::map<long, std::vector<long> >::iterator node =
                 layerNodes.begin(); node != layerNodes.end(); ++node) {
          long &levelNode = node->second[layer + 1];
          if (levelNode == firstNext)
            levelNode = secondNext;
          else if (levelNode == secondNext)
            levelNode = firstNext;
        }
        ++untwistedLayerConnectionCount;
        changed = true;
      }
      if (!changed)
        break;
    }
  }
  // Apply the same equivalence-class correction to the angular connections
  // inside every rounded cap. These quads use the reverse contour ordering
  // because they occupy the annulus between two circular fronts.
  for (long layer = 1; layer < layerCount; ++layer)
    for (std::map<long, std::vector<SelectedEndpointKey> >::const_iterator
             cap = roundedCornerHandles.begin();
         cap != roundedCornerHandles.end(); ++cap) {
      const long owner = selectedEdges[cap->second.front().first].owner;
      for (std::size_t pass = 0; pass < cap->second.size(); ++pass) {
        bool changed = false;
        for (std::size_t sector = 0; sector + 1 < cap->second.size();
             ++sector) {
          const SelectedEndpointKey firstKey = cap->second[sector];
          const SelectedEndpointKey secondKey = cap->second[sector + 1];
          const long firstCurrent = endpointLayerNodes[firstKey][layer];
          const long secondCurrent = endpointLayerNodes[secondKey][layer];
          const long firstNext = endpointLayerNodes[firstKey][layer + 1];
          const long secondNext = endpointLayerNodes[secondKey][layer + 1];
          if (firstNext == secondNext)
            continue;
          std::vector<long> currentQuad(4);
          currentQuad[0] = secondCurrent;
          currentQuad[1] = firstCurrent;
          currentQuad[2] = firstNext;
          currentQuad[3] = secondNext;
          const double currentArea = signedArea(currentQuad, pointList);
          if (simplePolygon(currentQuad, pointList, lengthTolerance,
                            areaTolerance) &&
              std::fabs(currentArea) > areaTolerance &&
              currentArea * originalArea[owner] > 0.0)
            continue;
          std::vector<long> swappedQuad(4);
          swappedQuad[0] = secondCurrent;
          swappedQuad[1] = firstCurrent;
          swappedQuad[2] = secondNext;
          swappedQuad[3] = firstNext;
          const double swappedArea = signedArea(swappedQuad, pointList);
          if (!simplePolygon(swappedQuad, pointList, lengthTolerance,
                             areaTolerance) ||
              std::fabs(swappedArea) <= areaTolerance ||
              swappedArea * originalArea[owner] <= 0.0)
            continue;
          for (std::map<SelectedEndpointKey, std::vector<long> >::iterator
                   endpoint = endpointLayerNodes.begin();
               endpoint != endpointLayerNodes.end(); ++endpoint) {
            long &node = endpoint->second[layer + 1];
            if (node == firstNext)
              node = secondNext;
            else if (node == secondNext)
              node = firstNext;
          }
          for (std::map<long, std::vector<long> >::iterator node =
                   layerNodes.begin();
               node != layerNodes.end(); ++node) {
            long &levelNode = node->second[layer + 1];
            if (levelNode == firstNext)
              levelNode = secondNext;
            else if (levelNode == secondNext)
              levelNode = firstNext;
          }
          ++untwistedLayerConnectionCount;
          changed = true;
        }
        if (!changed)
          break;
      }
    }

  // Build the deformed core. Sliding wall endpoints are replaced by their
  // innermost copies. A straight physical-group termination retains a
  // wall-normal internal cap, represented at every layer level.
  std::vector<std::vector<long> > coreCandidates(oldCellCount);
  for (long cell = 0; cell < oldCellCount; ++cell) {
    const std::vector<long> original = polygonList[cell];
    std::vector<long> core;
    for (std::size_t local = 0; local < original.size(); ++local) {
      const long previous =
          original[(local + original.size() - 1) % original.size()];
      const long current = original[local];
      const long next = original[(local + 1) % original.size()];
      if (selectedDisplacement.find(current) ==
          selectedDisplacement.end()) {
        appendUnique(core, current);
        continue;
      }
      const std::map<EdgeKey, long>::const_iterator previousEntry =
          selectedCellEdgeIndex[cell].find(edgeKey(previous, current));
      const std::map<EdgeKey, long>::const_iterator nextEntry =
          selectedCellEdgeIndex[cell].find(edgeKey(current, next));
      const bool previousSelected =
          previousEntry != selectedCellEdgeIndex[cell].end();
      const bool nextSelected =
          nextEntry != selectedCellEdgeIndex[cell].end();
      if (beveledCorner.find(current) != beveledCorner.end()) {
        if (!previousSelected || !nextSelected)
          ExecError("PdmtAddViscousLayers2D: a beveled corner must join two selected boundary segments");
        const SelectedEndpointKey previousKey(previousEntry->second,
                                              current);
        const SelectedEndpointKey nextKey(nextEntry->second, current);
        const std::vector<SelectedEndpointKey> &cap =
            roundedCornerHandles[current];
        if (cap.front() == previousKey && cap.back() == nextKey) {
          for (std::size_t handle = 0; handle < cap.size(); ++handle)
            appendUnique(core,
                         endpointLayerNodes[cap[handle]][layerCount]);
        } else if (cap.front() == nextKey && cap.back() == previousKey) {
          for (std::size_t handle = cap.size(); handle > 0; --handle)
            appendUnique(core,
                         endpointLayerNodes[cap[handle - 1]][layerCount]);
        } else {
          ExecError("PdmtAddViscousLayers2D: rounded cap does not match its core corner");
        }
      } else if (internalEndpoint.find(current) != internalEndpoint.end() &&
          previousSelected && !nextSelected) {
        for (long level = layerCount; level >= 0; --level)
          appendUnique(core, layerNodes[current][level]);
      } else if (internalEndpoint.find(current) != internalEndpoint.end() &&
                 !previousSelected && nextSelected) {
        for (long level = 0; level <= layerCount; ++level)
          appendUnique(core, layerNodes[current][level]);
      } else {
        appendUnique(core, layerNodes[current][layerCount]);
      }
    }
    if (core.size() > 1 && core.front() == core.back())
      core.pop_back();
    coreCandidates[cell].swap(core);
  }

  // A thick offset may move the layer interface completely through one or
  // more small boundary-adjacent dual cells. Absorb such a cell across the
  // internal edge crossed by its deformed outline. This removes that edge
  // from the core topology before the resulting (usually non-convex) region
  // is retained as one polygon. Only cells with the same material label may
  // be merged; a physical interface or another exterior boundary remains a
  // hard limit.
  std::vector<char> activeCore(oldCellCount, 1);
  std::vector<char> absorbedCoreRegion(oldCellCount, 0);
  std::vector<double> coreReferenceArea = originalArea;
  long absorbedCorePolygonCount = 0;
  for (;;) {
    long foldedCell = -1;
    for (long cell = 0; cell < oldCellCount; ++cell) {
      if (!activeCore[cell])
        continue;
      const double area = signedArea(coreCandidates[cell], pointList);
      if (!simplePolygon(coreCandidates[cell], pointList, lengthTolerance,
                         areaTolerance) ||
          std::fabs(area) <= areaTolerance ||
          area * coreReferenceArea[cell] <= 0.0) {
        foldedCell = cell;
        break;
      }
    }
    if (foldedCell < 0)
      break;

    std::size_t firstCrossing = 0;
    std::size_t secondCrossing = 0;
    const bool hasCrossing = firstSelfIntersection(
        coreCandidates[foldedCell], pointList, lengthTolerance,
        areaTolerance, firstCrossing, secondCrossing);

    std::map<EdgeKey, std::vector<long> > coreEdgeOwners;
    for (long cell = 0; cell < oldCellCount; ++cell) {
      if (!activeCore[cell])
        continue;
      for (std::size_t local = 0; local < coreCandidates[cell].size();
           ++local)
        coreEdgeOwners[edgeKey(
            coreCandidates[cell][local],
            coreCandidates[cell][(local + 1) %
                                 coreCandidates[cell].size()])]
            .push_back(cell);
    }

    const std::size_t crossingEdges[2] = {firstCrossing, secondCrossing};
    long neighbour = -1;
    EdgeKey absorbedEdge;
    if (hasCrossing) {
      for (long crossing = 0; crossing < 2 && neighbour < 0; ++crossing) {
        const std::size_t local = crossingEdges[crossing];
        const EdgeKey candidate = edgeKey(
            coreCandidates[foldedCell][local],
            coreCandidates[foldedCell][(local + 1) %
                                       coreCandidates[foldedCell].size()]);
        const std::vector<long> &owners = coreEdgeOwners[candidate];
        if (owners.size() != 2)
          continue;
        const long other = owners[0] == foldedCell ? owners[1] : owners[0];
        if (other != foldedCell && activeCore[other] &&
            oldCellLabels[other] == oldCellLabels[foldedCell]) {
          neighbour = other;
          absorbedEdge = candidate;
        }
      }
    } else {
      // A persistent sharp-end collapse can reduce a core polygon to zero
      // area without producing a crossing. Absorb it across its longest
      // same-material internal edge and retry the combined region.
      double longestInternalEdge = -1.0;
      for (std::size_t local = 0;
           local < coreCandidates[foldedCell].size(); ++local) {
        const EdgeKey candidate = edgeKey(
            coreCandidates[foldedCell][local],
            coreCandidates[foldedCell][(local + 1) %
                                       coreCandidates[foldedCell].size()]);
        const std::vector<long> &owners = coreEdgeOwners[candidate];
        if (owners.size() != 2)
          continue;
        const long other = owners[0] == foldedCell ? owners[1] : owners[0];
        if (other == foldedCell || !activeCore[other] ||
            oldCellLabels[other] != oldCellLabels[foldedCell])
          continue;
        const double dx = pointList[candidate.second].x -
                          pointList[candidate.first].x;
        const double dy = pointList[candidate.second].y -
                          pointList[candidate.first].y;
        const double lengthSquared = dx * dx + dy * dy;
        if (lengthSquared > longestInternalEdge) {
          longestInternalEdge = lengthSquared;
          neighbour = other;
          absorbedEdge = candidate;
        }
      }
    }
    if (neighbour < 0) {
      std::ostringstream message;
      message << "PdmtAddViscousLayers2D: thickness " << totalThickness
              << " makes the layer front collide with a boundary or "
              << "material interface near dual polygon " << foldedCell;
      if (hasCrossing)
        message << "; crossed edges";
      for (long crossing = 0; hasCrossing && crossing < 2; ++crossing) {
        const std::size_t local = crossingEdges[crossing];
        const EdgeKey candidate = edgeKey(
            coreCandidates[foldedCell][local],
            coreCandidates[foldedCell][(local + 1) %
                                       coreCandidates[foldedCell].size()]);
        message << " (" << candidate.first << "," << candidate.second
                << ", owners=" << coreEdgeOwners[candidate].size()
                << ", xy=" << pointList[candidate.first].x << ","
                << pointList[candidate.first].y << "->"
                << pointList[candidate.second].x << ","
                << pointList[candidate.second].y << ")";
      }
      message << "; core";
      for (std::size_t local = 0;
           local < coreCandidates[foldedCell].size(); ++local)
        message << " " << coreCandidates[foldedCell][local];
      ExecError(message.str());
    }

    coreCandidates[foldedCell] = mergePolygonsAcrossEdge(
        coreCandidates[foldedCell], coreCandidates[neighbour], absorbedEdge);
    coreReferenceArea[foldedCell] += coreReferenceArea[neighbour];
    absorbedCoreRegion[foldedCell] = 1;
    activeCore[neighbour] = 0;
    coreCandidates[neighbour].clear();
    ++absorbedCorePolygonCount;
  }

  std::vector<std::vector<long> > rebuiltCore;
  std::vector<long> rebuiltCoreLabels;
  long partitionedCorePolygonCount = 0;
  for (long cell = 0; cell < oldCellCount; ++cell) {
    if (!activeCore[cell])
      continue;
    // Preserve ordinary dual polygons, including their mild non-convexity.
    // An absorbed sharp-tip patch is different: retaining the whole patch as
    // one cell creates a large re-entrant polygon behind the rounded cap.
    // Partition only that patch, then greedily merge the ears into convex
    // polygon cells so the repair remains local and does not expose a fan of
    // unnecessary triangles.
    if (absorbedCoreRegion[cell] &&
        !convexPolygon(coreCandidates[cell], pointList, areaTolerance)) {
      const std::vector<std::vector<long> > pieces = convexPartition(
          coreCandidates[cell], pointList, lengthTolerance, areaTolerance);
      for (std::size_t piece = 0; piece < pieces.size(); ++piece) {
        rebuiltCore.push_back(pieces[piece]);
        rebuiltCoreLabels.push_back(oldCellLabels[cell]);
      }
      partitionedCorePolygonCount +=
          static_cast<long>(pieces.size()) - 1;
    } else {
      rebuiltCore.push_back(coreCandidates[cell]);
      rebuiltCoreLabels.push_back(oldCellLabels[cell]);
    }
  }

  polygonList.swap(rebuiltCore);
  std::vector<long> cellLabels;
  cellLabels.swap(rebuiltCoreLabels);
  std::vector<long> cellBands(polygonList.size(), -1);
  const long rebuiltCoreCellCount =
      static_cast<long>(polygonList.size());
  long addedLayerCells = 0;
  const auto appendLayerCell =
      [&](const std::vector<long> &rawPolygon, long owner,
          long band, const char *failureMessage) {
        std::vector<long> polygon;
        for (std::size_t local = 0; local < rawPolygon.size(); ++local)
          appendUnique(polygon, rawPolygon[local]);
        if (polygon.size() > 1 && polygon.front() == polygon.back())
          polygon.pop_back();
        // A front loop can disappear between consecutive levels. Its former
        // strip cell then has no area and is intentionally omitted.
        if (polygon.size() < 3)
          return;
        double area = signedArea(polygon, pointList);
        if (simplePolygon(polygon, pointList, lengthTolerance,
                          areaTolerance) &&
            std::fabs(area) > areaTolerance &&
            area * originalArea[owner] < 0.0) {
          std::reverse(polygon.begin(), polygon.end());
          area = -area;
        }
        if (!simplePolygon(polygon, pointList, lengthTolerance,
                           areaTolerance) ||
            std::fabs(area) <= areaTolerance ||
            area * originalArea[owner] <= 0.0) {
          std::ostringstream message;
          message << failureMessage << " at total thickness "
                  << totalThickness
                  << " in band " << band << " near dual polygon "
                  << owner << " (simple="
                  << simplePolygon(polygon, pointList, lengthTolerance,
                                   areaTolerance)
                  << ", area=" << area
                  << ", reference=" << originalArea[owner]
                  << ", vertices";
          for (std::size_t local = 0; local < polygon.size(); ++local)
            message << " " << polygon[local] << "=("
                    << pointList[polygon[local]].x << ","
                    << pointList[polygon[local]].y << ")";
          message << ")"
                  << "; the requested layer exceeds the local available "
                  << "space or requires additional boundary refinement";
          ExecError(message.str());
        }
        polygonList.push_back(polygon);
        cellLabels.push_back(oldCellLabels[owner]);
        cellBands.push_back(band);
        ++addedLayerCells;
      };
  for (std::size_t edge = 0; edge < selectedEdges.size(); ++edge) {
    const SelectedEdge &selected = selectedEdges[edge];
    const std::vector<long> &firstLevels = endpointLayerNodes[
        SelectedEndpointKey(static_cast<long>(edge), selected.first)];
    const std::vector<long> &secondLevels = endpointLayerNodes[
        SelectedEndpointKey(static_cast<long>(edge), selected.second)];
    for (long layer = 0; layer < layerCount; ++layer) {
      std::vector<long> quadrilateral(4);
      quadrilateral[0] = firstLevels[layer];
      quadrilateral[1] = secondLevels[layer];
      quadrilateral[2] = secondLevels[layer + 1];
      quadrilateral[3] = firstLevels[layer + 1];
      appendLayerCell(
          quadrilateral, selected.owner, layer,
          "PdmtAddViscousLayers2D: layer-front trimming creates a folded layer cell");
    }
  }

  // Fill the rounded cap one angular sector at a time. The first band starts
  // at the physical point and is created as a temporary fan; the merge pass
  // below combines adjacent sectors into the largest convex polygons. Every
  // later band remains a structured row of quadrilaterals.
  for (std::set<long>::const_iterator corner = beveledCorner.begin();
       corner != beveledCorner.end(); ++corner) {
    const std::vector<SelectedEndpointKey> &rays =
        roundedCornerHandles[*corner];
    const long owner = selectedEdges[rays.front().first].owner;
    for (std::size_t sector = 0; sector + 1 < rays.size(); ++sector) {
      const std::vector<long> &firstLevels =
          endpointLayerNodes[rays[sector]];
      const std::vector<long> &secondLevels =
          endpointLayerNodes[rays[sector + 1]];
      std::vector<long> triangle(3);
      triangle[0] = *corner;
      triangle[1] = firstLevels[1];
      triangle[2] = secondLevels[1];
      appendLayerCell(
          triangle, owner, 0,
          "PdmtAddViscousLayers2D: cannot construct a valid rounded corner cap");
      for (long layer = 1; layer < layerCount; ++layer) {
        std::vector<long> quadrilateral(4);
        quadrilateral[0] = secondLevels[layer];
        quadrilateral[1] = firstLevels[layer];
        quadrilateral[2] = firstLevels[layer + 1];
        quadrilateral[3] = secondLevels[layer + 1];
        appendLayerCell(
            quadrilateral, owner, layer,
            "PdmtAddViscousLayers2D: rounded corner cap contains a folded cell");
      }
    }
  }

  // Triangle cells are useful internally when a front loop disappears or a
  // bevel starts at one wall vertex, but they need not remain in the output.
  // Merge each one across a same-material layer edge whenever the union is a
  // simple positive polygon. This retains the minimal topology change while
  // keeping the visible boundary-layer mesh polygonal and compact.
  long mergedLayerTriangleCount = 0;
  std::vector<char> activeLayerCell(addedLayerCells, 1);
  bool mergedTriangle = true;
  while (mergedTriangle) {
    mergedTriangle = false;
    std::map<EdgeKey, std::vector<long> > layerEdgeOwners;
    for (long localCell = 0; localCell < addedLayerCells; ++localCell) {
      if (!activeLayerCell[localCell])
        continue;
      const long cell = rebuiltCoreCellCount + localCell;
      for (std::size_t local = 0; local < polygonList[cell].size(); ++local)
        layerEdgeOwners[edgeKey(
            polygonList[cell][local],
            polygonList[cell][(local + 1) % polygonList[cell].size()])]
            .push_back(localCell);
    }
    for (long triangle = 0; triangle < addedLayerCells && !mergedTriangle;
         ++triangle) {
      if (!activeLayerCell[triangle] ||
          polygonList[rebuiltCoreCellCount + triangle].size() != 3)
        continue;
      const long triangleCell = rebuiltCoreCellCount + triangle;
      for (std::size_t local = 0;
           local < polygonList[triangleCell].size(); ++local) {
        const EdgeKey shared = edgeKey(
            polygonList[triangleCell][local],
            polygonList[triangleCell][(local + 1) %
                                      polygonList[triangleCell].size()]);
        const std::vector<long> &owners = layerEdgeOwners[shared];
        if (owners.size() != 2)
          continue;
        const long neighbour =
            owners[0] == triangle ? owners[1] : owners[0];
        if (neighbour == triangle || !activeLayerCell[neighbour])
          continue;
        const long neighbourCell = rebuiltCoreCellCount + neighbour;
        if (cellLabels[triangleCell] != cellLabels[neighbourCell] ||
            cellBands[triangleCell] != cellBands[neighbourCell])
          continue;
        const std::vector<long> merged = mergePolygonsAcrossEdge(
            polygonList[neighbourCell], polygonList[triangleCell], shared);
        const double mergedArea = signedArea(merged, pointList);
        const double neighbourArea =
            signedArea(polygonList[neighbourCell], pointList);
        if (!simplePolygon(merged, pointList, lengthTolerance,
                           areaTolerance) ||
            !convexPolygon(merged, pointList, areaTolerance) ||
            std::fabs(mergedArea) <= areaTolerance ||
            mergedArea * neighbourArea <= 0.0)
          continue;
        polygonList[neighbourCell] = merged;
        activeLayerCell[triangle] = 0;
        ++mergedLayerTriangleCount;
        mergedTriangle = true;
        break;
      }
    }
  }
  if (mergedLayerTriangleCount > 0) {
    std::vector<std::vector<long> > compactPolygons;
    std::vector<long> compactLabels;
    std::vector<long> compactBands;
    compactPolygons.reserve(polygonList.size() - mergedLayerTriangleCount);
    compactLabels.reserve(cellLabels.size() - mergedLayerTriangleCount);
    compactBands.reserve(cellBands.size() - mergedLayerTriangleCount);
    for (long cell = 0; cell < rebuiltCoreCellCount; ++cell) {
      compactPolygons.push_back(polygonList[cell]);
      compactLabels.push_back(cellLabels[cell]);
      compactBands.push_back(cellBands[cell]);
    }
    for (long localCell = 0; localCell < addedLayerCells; ++localCell) {
      if (!activeLayerCell[localCell])
        continue;
      const long cell = rebuiltCoreCellCount + localCell;
      compactPolygons.push_back(polygonList[cell]);
      compactLabels.push_back(cellLabels[cell]);
      compactBands.push_back(cellBands[cell]);
    }
    polygonList.swap(compactPolygons);
    cellLabels.swap(compactLabels);
    cellBands.swap(compactBands);
    addedLayerCells -= mergedLayerTriangleCount;
  }

  // Very thick fronts can change topology between two requested levels. The
  // individual cells remain topologically conforming, but their interior
  // spokes may cross geometrically. Remove such spokes by agglomerating the
  // smallest connected patch containing both crossed edges. The patch is
  // emitted as one simple polygon; no triangulation is introduced.
  long agglomeratedCrossingCellCount = 0;
  std::vector<char> activeCell(polygonList.size(), 1);
  std::vector<char> layerCell(polygonList.size(), 0);
  for (std::size_t cell = 0; cell < polygonList.size(); ++cell)
    layerCell[cell] = cellBands[cell] >= 0;
  for (std::size_t repair = 0; repair < polygonList.size(); ++repair) {
    std::map<EdgeKey, std::vector<long> > allEdgeOwners;
    for (std::size_t cell = 0; cell < polygonList.size(); ++cell) {
      if (!activeCell[cell])
        continue;
      for (std::size_t local = 0; local < polygonList[cell].size(); ++local)
        allEdgeOwners[edgeKey(
            polygonList[cell][local],
            polygonList[cell][(local + 1) % polygonList[cell].size()])]
            .push_back(static_cast<long>(cell));
    }

    std::vector<EdgeKey> uniqueEdges;
    uniqueEdges.reserve(allEdgeOwners.size());
    for (std::map<EdgeKey, std::vector<long> >::const_iterator edge =
             allEdgeOwners.begin(); edge != allEdgeOwners.end(); ++edge)
      uniqueEdges.push_back(edge->first);
    const double binSize = std::max(
        extent / 200.0, totalThickness / static_cast<double>(layerCount));
    std::map<std::pair<long, long>, std::vector<long> > edgeBins;
    for (std::size_t edge = 0; edge < uniqueEdges.size(); ++edge) {
      const Point &first = pointList[uniqueEdges[edge].first];
      const Point &second = pointList[uniqueEdges[edge].second];
      const long firstX =
          static_cast<long>(std::floor(std::min(first.x, second.x) / binSize));
      const long lastX =
          static_cast<long>(std::floor(std::max(first.x, second.x) / binSize));
      const long firstY =
          static_cast<long>(std::floor(std::min(first.y, second.y) / binSize));
      const long lastY =
          static_cast<long>(std::floor(std::max(first.y, second.y) / binSize));
      for (long x = firstX; x <= lastX; ++x)
        for (long y = firstY; y <= lastY; ++y)
          edgeBins[std::make_pair(x, y)].push_back(
              static_cast<long>(edge));
    }

    long crossedFirst = -1;
    long crossedSecond = -1;
    std::set<std::pair<long, long> > comparedEdges;
    for (std::map<std::pair<long, long>, std::vector<long> >::const_iterator
             bin = edgeBins.begin();
         bin != edgeBins.end() && crossedFirst < 0; ++bin)
      for (std::size_t first = 0;
           first < bin->second.size() && crossedFirst < 0; ++first)
        for (std::size_t second = first + 1;
             second < bin->second.size(); ++second) {
          long firstEdge = bin->second[first];
          long secondEdge = bin->second[second];
          if (firstEdge > secondEdge)
            std::swap(firstEdge, secondEdge);
          if (!comparedEdges.insert(
                  std::make_pair(firstEdge, secondEdge)).second)
            continue;
          const EdgeKey &a = uniqueEdges[firstEdge];
          const EdgeKey &b = uniqueEdges[secondEdge];
          if (a.first == b.first || a.first == b.second ||
              a.second == b.first || a.second == b.second)
            continue;
          const double abc = orientation(pointList[a.first],
                                         pointList[a.second],
                                         pointList[b.first]);
          const double abd = orientation(pointList[a.first],
                                         pointList[a.second],
                                         pointList[b.second]);
          const double cda = orientation(pointList[b.first],
                                         pointList[b.second],
                                         pointList[a.first]);
          const double cdb = orientation(pointList[b.first],
                                         pointList[b.second],
                                         pointList[a.second]);
          if (((abc > areaTolerance && abd < -areaTolerance) ||
               (abc < -areaTolerance && abd > areaTolerance)) &&
              ((cda > areaTolerance && cdb < -areaTolerance) ||
               (cda < -areaTolerance && cdb > areaTolerance))) {
            crossedFirst = firstEdge;
            crossedSecond = secondEdge;
            break;
          }
        }
    if (crossedFirst < 0)
      break;

    const std::vector<long> &allFirstOwners =
        allEdgeOwners[uniqueEdges[crossedFirst]];
    const std::vector<long> &allSecondOwners =
        allEdgeOwners[uniqueEdges[crossedSecond]];
    long transitionBand = -1;
    for (std::size_t first = 0;
         first < allFirstOwners.size() && transitionBand < 0; ++first)
      for (std::size_t second = 0;
           second < allSecondOwners.size(); ++second)
        if (cellBands[allFirstOwners[first]] >= 0 &&
            cellBands[allFirstOwners[first]] ==
                cellBands[allSecondOwners[second]]) {
          transitionBand = cellBands[allFirstOwners[first]];
          break;
        }
    if (transitionBand < 0)
      ExecError("PdmtAddViscousLayers2D: adjacent layer bands cross near a sharp end; increase --viscous-layer_count or refine that boundary");
    std::vector<long> firstOwners;
    std::vector<long> secondOwners;
    for (std::size_t owner = 0; owner < allFirstOwners.size(); ++owner)
      if (cellBands[allFirstOwners[owner]] == transitionBand)
        firstOwners.push_back(allFirstOwners[owner]);
    for (std::size_t owner = 0; owner < allSecondOwners.size(); ++owner)
      if (cellBands[allSecondOwners[owner]] == transitionBand)
        secondOwners.push_back(allSecondOwners[owner]);
    if (firstOwners.empty() || secondOwners.empty())
      ExecError("PdmtAddViscousLayers2D: cannot isolate a crossing layer band");
    const long materialLabel = cellLabels[firstOwners[0]];
    std::set<long> targetOwners(secondOwners.begin(), secondOwners.end());
    std::vector<long> predecessor(polygonList.size(), -2);
    std::queue<long> pending;
    for (std::size_t owner = 0; owner < firstOwners.size(); ++owner) {
      predecessor[firstOwners[owner]] = -1;
      pending.push(firstOwners[owner]);
    }
    long reached = -1;
    while (!pending.empty() && reached < 0) {
      const long cell = pending.front();
      pending.pop();
      if (targetOwners.find(cell) != targetOwners.end()) {
        reached = cell;
        break;
      }
      for (std::size_t local = 0; local < polygonList[cell].size(); ++local) {
        const EdgeKey edge = edgeKey(
            polygonList[cell][local],
            polygonList[cell][(local + 1) % polygonList[cell].size()]);
        const std::vector<long> &owners = allEdgeOwners[edge];
        for (std::size_t owner = 0; owner < owners.size(); ++owner) {
          const long neighbour = owners[owner];
          if (neighbour == cell || !activeCell[neighbour] ||
              predecessor[neighbour] != -2 ||
              cellLabels[neighbour] != materialLabel ||
              cellBands[neighbour] != transitionBand)
            continue;
          predecessor[neighbour] = cell;
          pending.push(neighbour);
        }
      }
    }
    if (reached < 0)
      ExecError("PdmtAddViscousLayers2D: opposing layer fronts collide across a material or boundary interface");

    std::set<long> cluster;
    for (std::size_t owner = 0; owner < firstOwners.size(); ++owner)
      cluster.insert(firstOwners[owner]);
    for (std::size_t owner = 0; owner < secondOwners.size(); ++owner)
      cluster.insert(secondOwners[owner]);
    for (long cell = reached; cell >= 0; cell = predecessor[cell])
      cluster.insert(cell);

    std::vector<long> patchBoundary;
    for (std::size_t expansion = 0;
         expansion < polygonList.size(); ++expansion) {
      std::map<EdgeKey, long> patchEdgeUse;
      for (std::set<long>::const_iterator cell = cluster.begin();
           cell != cluster.end(); ++cell)
        for (std::size_t local = 0; local < polygonList[*cell].size(); ++local)
          ++patchEdgeUse[edgeKey(
              polygonList[*cell][local],
              polygonList[*cell][(local + 1) % polygonList[*cell].size()])];
      std::map<long, std::vector<long> > boundaryNeighbours;
      long boundaryEdgeCount = 0;
      for (std::map<EdgeKey, long>::const_iterator edge = patchEdgeUse.begin();
           edge != patchEdgeUse.end(); ++edge)
        if (edge->second == 1) {
          boundaryNeighbours[edge->first.first].push_back(edge->first.second);
          boundaryNeighbours[edge->first.second].push_back(edge->first.first);
          ++boundaryEdgeCount;
        }
      bool validCycle = !boundaryNeighbours.empty();
      for (std::map<long, std::vector<long> >::const_iterator node =
               boundaryNeighbours.begin();
           node != boundaryNeighbours.end(); ++node)
        validCycle = validCycle && node->second.size() == 2;
      patchBoundary.clear();
      if (validCycle) {
        const long start = boundaryNeighbours.begin()->first;
        long previous = -1;
        long current = start;
        for (long step = 0; step <= boundaryEdgeCount; ++step) {
          if (step > 0 && current == start)
            break;
          patchBoundary.push_back(current);
          const std::vector<long> &neighbours = boundaryNeighbours[current];
          const long next = neighbours[0] == previous
                                ? neighbours[1] : neighbours[0];
          previous = current;
          current = next;
        }
        validCycle = static_cast<long>(patchBoundary.size()) ==
                         boundaryEdgeCount &&
                     simplePolygon(patchBoundary, pointList,
                                   lengthTolerance, areaTolerance);
      }
      if (validCycle)
        break;

      bool enlarged = false;
      if (!patchBoundary.empty()) {
        std::size_t firstCrossing = 0;
        std::size_t secondCrossing = 0;
        if (firstSelfIntersection(patchBoundary, pointList,
                                  lengthTolerance, areaTolerance,
                                  firstCrossing, secondCrossing)) {
          const std::size_t crossings[2] = {firstCrossing,
                                             secondCrossing};
          for (long crossing = 0; crossing < 2; ++crossing) {
            const std::size_t local = crossings[crossing];
            const EdgeKey edge = edgeKey(
                patchBoundary[local],
                patchBoundary[(local + 1) % patchBoundary.size()]);
            const std::vector<long> &owners = allEdgeOwners[edge];
            for (std::size_t owner = 0; owner < owners.size(); ++owner)
              if (cluster.find(owners[owner]) == cluster.end() &&
                  cellLabels[owners[owner]] == materialLabel &&
                  cellBands[owners[owner]] == transitionBand) {
                cluster.insert(owners[owner]);
                enlarged = true;
              }
          }
        }
      }
      if (!enlarged)
        ExecError("PdmtAddViscousLayers2D: cannot form a simple polygonal patch around crossing layer fronts");
    }
    if (patchBoundary.size() < 3 ||
        !simplePolygon(patchBoundary, pointList, lengthTolerance,
                       areaTolerance))
      ExecError("PdmtAddViscousLayers2D: crossing-front agglomeration did not converge");

    double referenceArea = 0.0;
    const long anchor = *cluster.begin();
    for (std::set<long>::const_iterator cell = cluster.begin();
         cell != cluster.end(); ++cell) {
      referenceArea += signedArea(polygonList[*cell], pointList);
      if (*cell != anchor)
        activeCell[*cell] = 0;
    }
    if (signedArea(patchBoundary, pointList) * referenceArea < 0.0)
      std::reverse(patchBoundary.begin(), patchBoundary.end());
    polygonList[anchor] = patchBoundary;
    layerCell[anchor] = 1;
    cellBands[anchor] = transitionBand;
    agglomeratedCrossingCellCount +=
        static_cast<long>(cluster.size()) - 1;
  }

  if (agglomeratedCrossingCellCount > 0) {
    std::vector<std::vector<long> > compactPolygons;
    std::vector<long> compactLabels;
    std::vector<long> compactBands;
    for (long layerPass = 0; layerPass < 2; ++layerPass)
      for (std::size_t cell = 0; cell < polygonList.size(); ++cell)
        if (activeCell[cell] && layerCell[cell] == (layerPass == 1)) {
          compactPolygons.push_back(polygonList[cell]);
          compactLabels.push_back(cellLabels[cell]);
          compactBands.push_back(cellBands[cell]);
        }
    polygonList.swap(compactPolygons);
    cellLabels.swap(compactLabels);
    cellBands.swap(compactBands);
    addedLayerCells = 0;
    for (std::size_t cell = 0; cell < activeCell.size(); ++cell)
      if (activeCell[cell] && layerCell[cell])
        ++addedLayerCells;
  }

  // Rebuild the explicit boundary. At a sliding endpoint, the lateral edge
  // of every layer is part of the adjacent side boundary, and the remaining
  // old side segment starts at the innermost copy.
  std::vector<std::vector<long> > edgeList;
  std::vector<long> edgeLabels;
  for (long edge = 0; edge < oldEdgeCount; ++edge) {
    std::vector<long> outputEdge(2);
    for (long endpoint = 0; endpoint < 2; ++endpoint) {
      const long node = edges[edge][endpoint];
      if (!selectedExplicitEdge[edge] &&
          slidingEndpointEdge.find(node) != slidingEndpointEdge.end())
        outputEdge[endpoint] = layerNodes[node][layerCount];
      else
        outputEdge[endpoint] = node;
    }
    edgeList.push_back(outputEdge);
    edgeLabels.push_back(oldEdgeLabels[edge]);
  }
  for (std::map<long, long>::const_iterator endpoint =
           slidingEndpointEdge.begin();
       endpoint != slidingEndpointEdge.end(); ++endpoint)
    for (long layer = 0; layer < layerCount; ++layer) {
      std::vector<long> lateral(2);
      lateral[0] = layerNodes[endpoint->first][layer];
      lateral[1] = layerNodes[endpoint->first][layer + 1];
      edgeList.push_back(lateral);
      edgeLabels.push_back(oldEdgeLabels[endpoint->second]);
    }

  // A conforming result has every polygon edge once (exterior) or twice
  // (interior), and its exterior exactly matches the rebuilt boundary.
  std::map<EdgeKey, long> outputEdgeUse;
  for (std::size_t cell = 0; cell < polygonList.size(); ++cell)
    for (std::size_t local = 0; local < polygonList[cell].size(); ++local)
      ++outputEdgeUse[edgeKey(
          polygonList[cell][local],
          polygonList[cell][(local + 1) % polygonList[cell].size()])];
  std::set<EdgeKey> explicitBoundary;
  for (std::size_t edge = 0; edge < edgeList.size(); ++edge)
    explicitBoundary.insert(edgeKey(edgeList[edge][0], edgeList[edge][1]));
  std::set<EdgeKey> polygonBoundary;
  for (std::map<EdgeKey, long>::const_iterator edge = outputEdgeUse.begin();
       edge != outputEdgeUse.end(); ++edge) {
    if (edge->second > 2)
      ExecError("PdmtAddViscousLayers2D: generated non-manifold layer topology");
    if (edge->second == 1)
      polygonBoundary.insert(edge->first);
  }
  if (polygonBoundary != explicitBoundary)
    ExecError("PdmtAddViscousLayers2D: deformed layers do not match the rebuilt exterior boundary");

  // Loop trimming intentionally leaves superseded candidate points behind.
  // Remove all unused nodes so coincident inactive candidates cannot be
  // mistaken for disconnected sharp-end vertices by downstream readers.
  std::vector<long> nodeRemap(pointList.size(), -1);
  std::vector<char> usedNode(pointList.size(), 0);
  for (std::size_t cell = 0; cell < polygonList.size(); ++cell)
    for (std::size_t local = 0; local < polygonList[cell].size(); ++local)
      usedNode[polygonList[cell][local]] = 1;
  for (std::size_t edge = 0; edge < edgeList.size(); ++edge)
    for (std::size_t endpoint = 0; endpoint < edgeList[edge].size();
         ++endpoint)
      usedNode[edgeList[edge][endpoint]] = 1;
  std::vector<Point> compactPoints;
  compactPoints.reserve(pointList.size());
  for (std::size_t oldNode = 0; oldNode < pointList.size(); ++oldNode)
    if (usedNode[oldNode]) {
      nodeRemap[oldNode] = static_cast<long>(compactPoints.size());
      compactPoints.push_back(pointList[oldNode]);
    }
  for (std::size_t cell = 0; cell < polygonList.size(); ++cell)
    for (std::size_t local = 0; local < polygonList[cell].size(); ++local)
      polygonList[cell][local] = nodeRemap[polygonList[cell][local]];
  for (std::size_t edge = 0; edge < edgeList.size(); ++edge)
    for (std::size_t endpoint = 0; endpoint < edgeList[edge].size();
         ++endpoint)
      edgeList[edge][endpoint] = nodeRemap[edgeList[edge][endpoint]];
  pointList.swap(compactPoints);

  nodes.resize(static_cast<long>(pointList.size()), 2);
  for (long node = 0; node < static_cast<long>(pointList.size()); ++node) {
    nodes(node, 0L) = pointList[node].x;
    nodes(node, 1L) = pointList[node].y;
  }
  cells.resize(static_cast<long>(polygonList.size()));
  for (long cell = 0; cell < static_cast<long>(polygonList.size()); ++cell) {
    cells[cell].resize(static_cast<long>(polygonList[cell].size()));
    for (long local = 0;
         local < static_cast<long>(polygonList[cell].size()); ++local)
      cells[cell][local] = polygonList[cell][local];
  }
  edges.resize(static_cast<long>(edgeList.size()));
  for (long edge = 0; edge < static_cast<long>(edgeList.size()); ++edge) {
    edges[edge].resize(2);
    edges[edge][0] = edgeList[edge][0];
    edges[edge][1] = edgeList[edge][1];
  }
  labels.resize(
      static_cast<long>(polygonList.size() + edgeList.size()));
  for (long cell = 0; cell < static_cast<long>(cellLabels.size()); ++cell)
    labels[cell] = cellLabels[cell];
  for (long edge = 0; edge < static_cast<long>(edgeLabels.size()); ++edge)
    labels[static_cast<long>(polygonList.size()) + edge] =
        edgeLabels[edge];

  if (verbosity) {
    std::cout << "PDMT 2D viscous layers: added "
              << addedLayerCells << " layer cells in " << layerCount
              << " equal layers across " << selectedEdges.size()
              << " selected dual boundary segments; total thickness "
              << totalThickness
              << " (remaining dual mesh deformed inward";
    if (!beveledCorner.empty())
      std::cout << "; replaced " << beveledCorner.size()
                << " sharp-corner miter(s) with rounded caps using "
                << roundedCornerRayCount << " intermediate ray(s)";
    if (trimmedFrontLoopCount > 0)
      std::cout << "; trimmed " << trimmedFrontLoopCount
                << " offset-front loop(s) into conforming fans";
    if (untwistedLayerConnectionCount > 0)
      std::cout << "; untwisted " << untwistedLayerConnectionCount
                << " sharp-end layer connection(s)";
    if (absorbedCorePolygonCount > 0)
      std::cout << "; absorbed " << absorbedCorePolygonCount
                << " crossed core polygon(s)";
    if (partitionedCorePolygonCount > 0)
      std::cout << "; repartitioned absorbed core patches into "
                << partitionedCorePolygonCount
                << " additional convex polygon(s)";
    if (mergedLayerTriangleCount > 0)
      std::cout << "; merged " << mergedLayerTriangleCount
                << " transition triangle(s) into adjacent layer polygons";
    if (agglomeratedCrossingCellCount > 0)
      std::cout << "; agglomerated " << agglomeratedCrossingCellCount
                << " crossing transition cell(s) into simple polygons";
    std::cout << ")" << std::endl;
  }
  return addedLayerCells;
}

} // namespace Pdmt2D

class pdmtAddViscousLayers2D_Op : public E_F0mps {
public:
  Expression mesh;

  static const int n_name_param = 9;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];

  pdmtAddViscousLayers2D_Op(const basicAC_F0 &args, Expression inputMesh)
      : mesh(inputMesh) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator()(Stack stack) const;
};

basicAC_F0::name_and_type pdmtAddViscousLayers2D_Op::name_param[] = {
    {"nodes", &typeid(KNM<double> *)},
    {"cells", &typeid(KN<KN<long> > *)},
    {"edges", &typeid(KN<KN<long> > *)},
    {"labels", &typeid(KN<long> *)},
    {"meshFile", &typeid(std::string *)},
    {"medMeshName", &typeid(std::string *)},
    {"viscousLayerGroups", &typeid(std::string *)},
    {"viscousLayerCount", &typeid(long)},
    {"viscousLayerThickness", &typeid(double)}};

class pdmtAddViscousLayers2D : public OneOperator {
public:
  pdmtAddViscousLayers2D()
      : OneOperator(atype<long>(), atype<pmesh>()) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new pdmtAddViscousLayers2D_Op(
        args, t[0]->CastTo(args[0]));
  }
};

AnyType pdmtAddViscousLayers2D_Op::operator()(Stack stack) const {
  const pmesh pTh = GetAny<pmesh>((*mesh)(stack));
  if (!pTh)
    ExecError("PdmtAddViscousLayers2D: null input mesh");
  for (int argument = 0; argument < n_name_param; ++argument)
    if (!nargs[argument])
      ExecError("PdmtAddViscousLayers2D: all named arguments are required");

  KNM<double> *nodes =
      GetAny<KNM<double> *>((*nargs[0])(stack));
  KN<KN<long> > *cells =
      GetAny<KN<KN<long> > *>((*nargs[1])(stack));
  KN<KN<long> > *edges =
      GetAny<KN<KN<long> > *>((*nargs[2])(stack));
  KN<long> *labels =
      GetAny<KN<long> *>((*nargs[3])(stack));
  const std::string meshFile =
      *GetAny<std::string *>((*nargs[4])(stack));
  const std::string medMeshName =
      *GetAny<std::string *>((*nargs[5])(stack));
  const std::string groups =
      *GetAny<std::string *>((*nargs[6])(stack));
  const long count = GetAny<long>((*nargs[7])(stack));
  const double thickness = GetAny<double>((*nargs[8])(stack));

  if (!nodes || !cells || !edges || !labels)
    ExecError("PdmtAddViscousLayers2D: null output array");
  return Pdmt2D::addViscousLayers(
      *pTh, meshFile, medMeshName, groups, count, thickness,
      *nodes, *cells, *edges, *labels);
}

#endif
