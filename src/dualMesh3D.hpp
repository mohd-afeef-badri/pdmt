/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     Build a conforming barycentric dual of a tetrahedral FreeFEM mesh.

*****************************************************************************/

#ifndef PDMT_DUAL_MESH_3D_HPP
#define PDMT_DUAL_MESH_3D_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <vector>

namespace Pdmt3D {

struct Point {
  double x, y, z;
};

typedef std::array<long, 2> EdgeKey;
typedef std::array<long, 3> FaceKey;

inline EdgeKey edgeKey(long a, long b) {
  EdgeKey key = {{a, b}};
  std::sort(key.begin(), key.end());
  return key;
}

inline FaceKey faceKey(long a, long b, long c) {
  FaceKey key = {{a, b, c}};
  std::sort(key.begin(), key.end());
  return key;
}

inline Point minus(const Point &a, const Point &b) {
  Point p = {a.x - b.x, a.y - b.y, a.z - b.z};
  return p;
}

inline Point cross(const Point &a, const Point &b) {
  Point p = {a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x};
  return p;
}

inline double dot(const Point &a, const Point &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

struct LocalFace {
  std::array<long, 3> oriented;
  int count;
  LocalFace() : oriented({{0, 0, 0}}), count(0) {}
};

inline bool sameOrientation(const std::array<long, 3> &a,
                            const std::array<long, 3> &b) {
  return (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]) ||
         (a[0] == b[1] && a[1] == b[2] && a[2] == b[0]) ||
         (a[0] == b[2] && a[1] == b[0] && a[2] == b[1]);
}

inline void addOrientedSubTetFaces(
    const std::array<long, 4> &subTet,
    const std::vector<Point> &points,
    std::map<FaceKey, LocalFace> &localFaces) {
  static const int oppositeFace[4][3] = {
      {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

  for (int opposite = 0; opposite < 4; ++opposite) {
    std::array<long, 3> tri = {{
        subTet[oppositeFace[opposite][0]],
        subTet[oppositeFace[opposite][1]],
        subTet[oppositeFace[opposite][2]]}};

    const Point ab = minus(points[tri[1]], points[tri[0]]);
    const Point ac = minus(points[tri[2]], points[tri[0]]);
    const Point ao = minus(points[subTet[opposite]], points[tri[0]]);
    if (dot(cross(ab, ac), ao) > 0.0)
      std::swap(tri[1], tri[2]);

    const FaceKey key = faceKey(tri[0], tri[1], tri[2]);
    LocalFace &entry = localFaces[key];
    if (entry.count == 0)
      entry.oriented = tri;
    ++entry.count;
  }
}

inline Point add(const Point &a, const Point &b) {
  Point p = {a.x + b.x, a.y + b.y, a.z + b.z};
  return p;
}

inline Point scale(const Point &a, double value) {
  Point p = {a.x * value, a.y * value, a.z * value};
  return p;
}

inline double norm(const Point &a) {
  return std::sqrt(dot(a, a));
}

inline Point polygonAreaVector(const std::vector<long> &polygon,
                               const std::vector<Point> &points) {
  Point area = {0.0, 0.0, 0.0};
  for (long i = 0; i < static_cast<long>(polygon.size()); ++i) {
    const Point &a = points[polygon[i]];
    const Point &b = points[polygon[(i + 1) % polygon.size()]];
    area = add(area, cross(a, b));
  }
  return area;
}

inline std::vector<std::vector<long> > connectedTriangleComponents(
    const std::vector<long> &faceIds,
    const std::vector<std::array<long, 3> > &faces,
    const std::set<EdgeKey> &splitEdges) {
  std::map<EdgeKey, std::vector<long> > edgeFaces;
  for (long local = 0; local < static_cast<long>(faceIds.size()); ++local) {
    const std::array<long, 3> &face = faces[faceIds[local]];
    for (int edge = 0; edge < 3; ++edge)
      edgeFaces[edgeKey(face[edge], face[(edge + 1) % 3])].push_back(local);
  }

  std::vector<std::vector<long> > adjacency(faceIds.size());
  for (std::map<EdgeKey, std::vector<long> >::const_iterator edge = edgeFaces.begin();
       edge != edgeFaces.end(); ++edge) {
    if (splitEdges.find(edge->first) != splitEdges.end())
      continue;
    const std::vector<long> &incident = edge->second;
    for (long i = 0; i < static_cast<long>(incident.size()); ++i)
      for (long j = i + 1; j < static_cast<long>(incident.size()); ++j) {
        adjacency[incident[i]].push_back(incident[j]);
        adjacency[incident[j]].push_back(incident[i]);
      }
  }

  std::vector<char> visited(faceIds.size(), 0);
  std::vector<std::vector<long> > components;
  for (long seed = 0; seed < static_cast<long>(faceIds.size()); ++seed) {
    if (visited[seed])
      continue;
    components.push_back(std::vector<long>());
    std::vector<long> work(1, seed);
    visited[seed] = 1;
    while (!work.empty()) {
      const long local = work.back();
      work.pop_back();
      components.back().push_back(faceIds[local]);
      for (std::vector<long>::const_iterator next = adjacency[local].begin();
           next != adjacency[local].end(); ++next)
        if (!visited[*next]) {
          visited[*next] = 1;
          work.push_back(*next);
        }
    }
  }
  return components;
}

inline std::vector<long> componentBoundaryLoop(
    const std::vector<long> &component,
    const std::vector<std::array<long, 3> > &faces) {
  std::map<EdgeKey, int> edgeCount;
  for (std::vector<long>::const_iterator id = component.begin(); id != component.end(); ++id) {
    const std::array<long, 3> &face = faces[*id];
    for (int edge = 0; edge < 3; ++edge)
      ++edgeCount[edgeKey(face[edge], face[(edge + 1) % 3])];
  }

  std::map<long, std::vector<long> > boundaryAdjacency;
  for (std::map<EdgeKey, int>::const_iterator edge = edgeCount.begin();
       edge != edgeCount.end(); ++edge) {
    if (edge->second == 1) {
      boundaryAdjacency[edge->first[0]].push_back(edge->first[1]);
      boundaryAdjacency[edge->first[1]].push_back(edge->first[0]);
    } else if (edge->second != 2) {
      ExecError("PdmtBuildDual3D: non-manifold triangle fan while merging a dual face");
    }
  }
  if (boundaryAdjacency.size() < 3)
    ExecError("PdmtBuildDual3D: a merged dual face has fewer than three boundary vertices");
  for (std::map<long, std::vector<long> >::const_iterator vertex = boundaryAdjacency.begin();
       vertex != boundaryAdjacency.end(); ++vertex)
    if (vertex->second.size() != 2)
      ExecError("PdmtBuildDual3D: a triangle fan does not have one simple boundary loop");

  std::vector<long> loop;
  const long start = boundaryAdjacency.begin()->first;
  long previous = -1;
  long current = start;
  do {
    loop.push_back(current);
    const std::vector<long> &neighbours = boundaryAdjacency[current];
    const long next = neighbours[0] == previous ? neighbours[1] : neighbours[0];
    previous = current;
    current = next;
    if (loop.size() > boundaryAdjacency.size())
      ExecError("PdmtBuildDual3D: failed to order a merged dual face");
  } while (current != start);

  if (loop.size() != boundaryAdjacency.size())
    ExecError("PdmtBuildDual3D: a merged dual face contains multiple boundary loops");
  return loop;
}

inline void mergeTriangleFans(
    const std::vector<Point> &points,
    const std::vector<std::array<long, 3> > &triangles,
    const std::vector<long> &triangleLabels,
    const std::vector<std::vector<long> > &triangleCells,
    const std::vector<char> &removablePoint,
    const std::set<EdgeKey> &boundarySplitEdges,
    std::vector<std::vector<long> > &polygons,
    std::vector<long> &polygonLabels,
    std::vector<std::vector<long> > &polygonCells) {
  typedef std::array<long, 3> GroupKey;
  std::vector<std::vector<std::pair<long, int> > > uses(triangles.size());
  for (long cell = 0; cell < static_cast<long>(triangleCells.size()); ++cell)
    for (std::vector<long>::const_iterator encoded = triangleCells[cell].begin();
         encoded != triangleCells[cell].end(); ++encoded) {
      const long face = std::labs(*encoded) - 1;
      uses[face].push_back(std::make_pair(cell, *encoded > 0 ? 1 : -1));
    }

  std::map<GroupKey, std::vector<long> > groups;
  for (long face = 0; face < static_cast<long>(triangles.size()); ++face) {
    if (uses[face].size() == 1) {
      const GroupKey key = {{1, uses[face][0].first, triangleLabels[face]}};
      groups[key].push_back(face);
    } else if (uses[face].size() == 2) {
      const long a = std::min(uses[face][0].first, uses[face][1].first);
      const long b = std::max(uses[face][0].first, uses[face][1].first);
      const GroupKey key = {{0, a, b}};
      groups[key].push_back(face);
    } else {
      ExecError("PdmtBuildDual3D: a triangular dual face has invalid cell incidence");
    }
  }

  polygonCells.assign(triangleCells.size(), std::vector<long>());
  for (std::map<GroupKey, std::vector<long> >::const_iterator group = groups.begin();
       group != groups.end(); ++group) {
    const std::set<EdgeKey> noSplitEdges;
    const std::vector<std::vector<long> > components = connectedTriangleComponents(
        group->second, triangles,
        group->first[0] ? boundarySplitEdges : noSplitEdges);
    for (std::vector<std::vector<long> >::const_iterator component = components.begin();
         component != components.end(); ++component) {
      std::vector<long> loop = componentBoundaryLoop(*component, triangles);
      std::vector<long> simplified;
      for (std::vector<long>::const_iterator vertex = loop.begin(); vertex != loop.end(); ++vertex)
        if (!removablePoint[*vertex])
          simplified.push_back(*vertex);
      if (simplified.size() >= 3)
        loop.swap(simplified);

      const long owner = group->first[1];
      Point desiredArea = {0.0, 0.0, 0.0};
      for (std::vector<long>::const_iterator oldFace = component->begin();
           oldFace != component->end(); ++oldFace) {
        int sign = 0;
        for (std::vector<std::pair<long, int> >::const_iterator use = uses[*oldFace].begin();
             use != uses[*oldFace].end(); ++use)
          if (use->first == owner)
            sign = use->second;
        if (!sign)
          ExecError("PdmtBuildDual3D: cannot orient a merged dual face");
        const std::array<long, 3> &tri = triangles[*oldFace];
        const Point ab = minus(points[tri[1]], points[tri[0]]);
        const Point ac = minus(points[tri[2]], points[tri[0]]);
        desiredArea = add(desiredArea, scale(cross(ab, ac), static_cast<double>(sign)));
      }
      if (dot(polygonAreaVector(loop, points), desiredArea) < 0.0)
        std::reverse(loop.begin(), loop.end());

      const long newFace = static_cast<long>(polygons.size());
      polygons.push_back(loop);
      polygonLabels.push_back(group->first[0] ? group->first[2] : 0);
      polygonCells[owner].push_back(newFace + 1);
      if (group->first[0] == 0)
        polygonCells[group->first[2]].push_back(-(newFace + 1));
    }
  }
}

} // namespace Pdmt3D

#include "gmshFeatureEdges.hpp"
#ifdef MEDCOUPLING
#include "medFeatureEdges.hpp"
#endif

class pdmtBuildDual3D_Op : public E_F0mps {
public:
  Expression mesh;

  static const int n_name_param = 10;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];

  pdmtBuildDual3D_Op(const basicAC_F0 &args, Expression inputMesh)
      : mesh(inputMesh) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator()(Stack stack) const;
};

basicAC_F0::name_and_type pdmtBuildDual3D_Op::name_param[] = {
    {"nodes", &typeid(KNM<double> *)},
    {"faces", &typeid(KN<KN<long> > *)},
    {"cells", &typeid(KN<KN<long> > *)},
    {"labels", &typeid(KN<long> *)},
    {"faceLabels", &typeid(KN<long> *)},
    {"featureAngle", &typeid(double)},
    {"meshFile", &typeid(std::string *)},
    {"conserveEdge", &typeid(std::string *)},
    {"mode", &typeid(std::string *)},
    {"medMeshName", &typeid(std::string *)}};

class pdmtBuildDual3D : public OneOperator {
public:
  pdmtBuildDual3D() : OneOperator(atype<long>(), atype<pmesh3>()) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new pdmtBuildDual3D_Op(args, t[0]->CastTo(args[0]));
  }
};

AnyType pdmtBuildDual3D_Op::operator()(Stack stack) const {
  using namespace Pdmt3D;

  const pmesh3 pTh = GetAny<pmesh3>((*mesh)(stack));
  if (!pTh)
    ExecError("PdmtBuildDual3D: null input mesh");
  if (!nargs[0] || !nargs[1] || !nargs[2])
    ExecError("PdmtBuildDual3D: nodes, faces and cells output arrays are required");

  KNM<double> *nodes = GetAny<KNM<double> *>((*nargs[0])(stack));
  KN<KN<long> > *faces = GetAny<KN<KN<long> > *>((*nargs[1])(stack));
  KN<KN<long> > *cells = GetAny<KN<KN<long> > *>((*nargs[2])(stack));
  KN<long> *labels = nargs[3] ? GetAny<KN<long> *>((*nargs[3])(stack)) : 0;
  KN<long> *faceLabels = nargs[4] ? GetAny<KN<long> *>((*nargs[4])(stack)) : 0;
  const double featureAngle = nargs[5] ? GetAny<double>((*nargs[5])(stack)) : 45.0;
  std::string meshFile;
  std::string conserveEdge;
  std::string mode = "smooth_dual";
  std::string medMeshName;
  if (nargs[6])
    meshFile = *GetAny<std::string *>((*nargs[6])(stack));
  if (nargs[7])
    conserveEdge = *GetAny<std::string *>((*nargs[7])(stack));
  if (nargs[8])
    mode = *GetAny<std::string *>((*nargs[8])(stack));
  if (nargs[9])
    medMeshName = *GetAny<std::string *>((*nargs[9])(stack));

  const Mesh3 &Th = *pTh;
  if (Th.nt == 0 || Th.nv == 0)
    ExecError("PdmtBuildDual3D: the tetrahedral mesh is empty");
  if (featureAngle < 0.0 || featureAngle > 180.0)
    ExecError("PdmtBuildDual3D: featureAngle must be between 0 and 180 degrees");
  if (mode != "subdivided_dual" && mode != "smooth_dual")
    ExecError("PdmtBuildDual3D: mode must be subdivided_dual or smooth_dual");
  const bool smoothDual = mode == "smooth_dual";

  std::set<EdgeKey> primalEdges;
  std::set<FaceKey> primalFaces;
  std::vector<std::vector<long> > incidentTets(Th.nv);
  std::vector<long> cellRegion(Th.nv, 0);
  std::vector<char> regionSet(Th.nv, 0);

  static const int tetEdges[6][2] = {
      {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  static const int tetFaces[4][3] = {
      {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};

  for (long t = 0; t < Th.nt; ++t) {
    long vertex[4];
    for (int i = 0; i < 4; ++i) {
      vertex[i] = Th(Th[t][i]);
      incidentTets[vertex[i]].push_back(t);
      if (!regionSet[vertex[i]]) {
        cellRegion[vertex[i]] = Th[t].lab;
        regionSet[vertex[i]] = 1;
      }
    }
    for (int i = 0; i < 6; ++i)
      primalEdges.insert(edgeKey(vertex[tetEdges[i][0]], vertex[tetEdges[i][1]]));
    for (int i = 0; i < 4; ++i)
      primalFaces.insert(faceKey(vertex[tetFaces[i][0]], vertex[tetFaces[i][1]],
                                 vertex[tetFaces[i][2]]));
  }

  std::set<EdgeKey> conservedPrimalEdges;
  if (!splitPhysicalNames(conserveEdge).empty()) {
    if (meshFile.empty())
      ExecError("PdmtBuildDual3D: meshFile is required when conserveEdge is set");
    if (meshFile.find(".med") != std::string::npos) {
#ifdef MEDCOUPLING
      conservedPrimalEdges = readMedFeatureEdges(
          meshFile, medMeshName, conserveEdge, Th, primalEdges, -2,
          "PdmtBuildDual3D");
#else
      ExecError("PdmtBuildDual3D: MED conserveEdge requires MEDCoupling support");
#endif
    } else {
      conservedPrimalEdges = readGmshFeatureEdges(meshFile, conserveEdge, Th, primalEdges);
    }
  }

  std::vector<Point> pointList;
  pointList.reserve(Th.nv + primalEdges.size() + primalFaces.size() + Th.nt);
  for (long v = 0; v < Th.nv; ++v) {
    Point p = {Th(v).x, Th(v).y, Th(v).z};
    pointList.push_back(p);
  }

  std::map<EdgeKey, long> edgeNodes;
  for (std::set<EdgeKey>::const_iterator it = primalEdges.begin();
       it != primalEdges.end(); ++it) {
    const Point a = pointList[(*it)[0]];
    const Point b = pointList[(*it)[1]];
    Point p = {(a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5};
    edgeNodes[*it] = static_cast<long>(pointList.size());
    pointList.push_back(p);
  }

  std::map<FaceKey, long> faceNodes;
  for (std::set<FaceKey>::const_iterator it = primalFaces.begin();
       it != primalFaces.end(); ++it) {
    const Point a = pointList[(*it)[0]];
    const Point b = pointList[(*it)[1]];
    const Point c = pointList[(*it)[2]];
    Point p = {(a.x + b.x + c.x) / 3.0,
               (a.y + b.y + c.y) / 3.0,
               (a.z + b.z + c.z) / 3.0};
    faceNodes[*it] = static_cast<long>(pointList.size());
    pointList.push_back(p);
  }

  std::vector<long> tetNodes(Th.nt);
  for (long t = 0; t < Th.nt; ++t) {
    Point p = {0.0, 0.0, 0.0};
    for (int i = 0; i < 4; ++i) {
      const long v = Th(Th[t][i]);
      p.x += pointList[v].x;
      p.y += pointList[v].y;
      p.z += pointList[v].z;
    }
    p.x *= 0.25;
    p.y *= 0.25;
    p.z *= 0.25;
    tetNodes[t] = static_cast<long>(pointList.size());
    pointList.push_back(p);
  }

  std::map<FaceKey, long> boundaryTriangleLabels;
  std::set<FaceKey> boundaryPrimalFaces;
  std::map<EdgeKey, std::vector<long> > boundaryEdgeFaces;
  for (long b = 0; b < Th.nbe; ++b) {
    long v[3] = {Th(Th.be(b)[0]), Th(Th.be(b)[1]), Th(Th.be(b)[2])};
    const FaceKey primalFace = faceKey(v[0], v[1], v[2]);
    boundaryPrimalFaces.insert(primalFace);
    const long fc = faceNodes[primalFace];
    for (int i = 0; i < 3; ++i) {
      const long current = v[i];
      const long next = v[(i + 1) % 3];
      const long previous = v[(i + 2) % 3];
      boundaryEdgeFaces[edgeKey(current, next)].push_back(b);
      boundaryTriangleLabels[faceKey(current, edgeNodes[edgeKey(current, next)], fc)] = Th.be(b).lab;
      boundaryTriangleLabels[faceKey(current, fc, edgeNodes[edgeKey(current, previous)])] = Th.be(b).lab;
    }
  }

  std::vector<char> removablePoint(pointList.size(), 0);
  if (smoothDual)
    for (std::map<FaceKey, long>::const_iterator face = faceNodes.begin();
         face != faceNodes.end(); ++face)
      if (boundaryPrimalFaces.find(face->first) == boundaryPrimalFaces.end())
        removablePoint[face->second] = 1;

  const double pi = 4.0 * std::atan(1.0);
  const double featureCos = std::cos(featureAngle * pi / 180.0);
  std::set<EdgeKey> boundarySplitEdges;
  for (std::map<EdgeKey, std::vector<long> >::const_iterator edge = boundaryEdgeFaces.begin();
       edge != boundaryEdgeFaces.end(); ++edge) {
    bool feature = conservedPrimalEdges.find(edge->first) != conservedPrimalEdges.end() ||
                   edge->second.size() != 2;
    if (!feature) {
      const long f0 = edge->second[0];
      const long f1 = edge->second[1];
      feature = Th.be(f0).lab != Th.be(f1).lab;
      if (!feature) {
        const long a0 = Th(Th.be(f0)[0]);
        const long b0 = Th(Th.be(f0)[1]);
        const long c0 = Th(Th.be(f0)[2]);
        const long a1 = Th(Th.be(f1)[0]);
        const long b1 = Th(Th.be(f1)[1]);
        const long c1 = Th(Th.be(f1)[2]);
        const Point n0 = cross(Pdmt3D::minus(pointList[b0], pointList[a0]),
                               Pdmt3D::minus(pointList[c0], pointList[a0]));
        const Point n1 = cross(Pdmt3D::minus(pointList[b1], pointList[a1]),
                               Pdmt3D::minus(pointList[c1], pointList[a1]));
        const double denominator = norm(n0) * norm(n1);
        if (denominator == 0.0)
          ExecError("PdmtBuildDual3D: degenerate boundary triangle");
        feature = dot(n0, n1) / denominator < featureCos;
      }
    }
    if (feature) {
      const long midpoint = edgeNodes[edge->first];
      boundarySplitEdges.insert(edgeKey(edge->first[0], midpoint));
      boundarySplitEdges.insert(edgeKey(edge->first[1], midpoint));
    } else if (smoothDual) {
      removablePoint[edgeNodes[edge->first]] = 1;
    }
  }

  for (std::set<EdgeKey>::const_iterator edge = conservedPrimalEdges.begin();
       edge != conservedPrimalEdges.end(); ++edge) {
    const long midpoint = edgeNodes[*edge];
    boundarySplitEdges.insert(edgeKey((*edge)[0], midpoint));
    boundarySplitEdges.insert(edgeKey((*edge)[1], midpoint));
    if (boundaryEdgeFaces.find(*edge) == boundaryEdgeFaces.end() && verbosity)
      std::cout << "PDMT 3D dual: conserved edge was not present as a "
                << "FreeFEM boundary edge; forcing the dual split from the "
                << "MED/Gmsh edge group" << std::endl;
  }

  std::vector<std::array<long, 3> > globalFaces;
  std::vector<long> globalFaceLabels;
  std::map<FaceKey, long> globalFaceIds;
  std::vector<std::vector<long> > cellFaceIds(Th.nv);

  for (long v = 0; v < Th.nv; ++v) {
    std::map<FaceKey, LocalFace> localFaces;

    for (std::vector<long>::const_iterator ti = incidentTets[v].begin();
         ti != incidentTets[v].end(); ++ti) {
      const long t = *ti;
      long tv[4];
      for (int i = 0; i < 4; ++i)
        tv[i] = Th(Th[t][i]);

      for (int ui = 0; ui < 4; ++ui) {
        const long u = tv[ui];
        if (u == v)
          continue;
        const long ec = edgeNodes[edgeKey(v, u)];
        for (int wi = 0; wi < 4; ++wi) {
          const long w = tv[wi];
          if (w == v || w == u)
            continue;
          const long fc = faceNodes[faceKey(v, u, w)];
          const std::array<long, 4> subTet = {{v, ec, fc, tetNodes[t]}};
          addOrientedSubTetFaces(subTet, pointList, localFaces);
        }
      }
    }

    for (std::map<FaceKey, LocalFace>::const_iterator it = localFaces.begin();
         it != localFaces.end(); ++it) {
      if (it->second.count == 2)
        continue;
      if (it->second.count != 1) {
        std::ostringstream message;
        message << "PdmtBuildDual3D: non-manifold barycentric face at primal vertex " << v;
        ExecError(message.str());
      }

      const std::array<long, 3> tri = it->second.oriented;
      std::map<FaceKey, long>::iterator existing = globalFaceIds.find(it->first);
      long faceId;
      long encodedFaceId;
      if (existing == globalFaceIds.end()) {
        faceId = static_cast<long>(globalFaces.size());
        globalFaceIds[it->first] = faceId;
        globalFaces.push_back(tri);
        std::map<FaceKey, long>::const_iterator boundary = boundaryTriangleLabels.find(it->first);
        globalFaceLabels.push_back(boundary == boundaryTriangleLabels.end() ? 0 : boundary->second);
        encodedFaceId = faceId + 1;
      } else {
        faceId = existing->second;
        encodedFaceId = sameOrientation(tri, globalFaces[faceId]) ? faceId + 1 : -(faceId + 1);
      }
      cellFaceIds[v].push_back(encodedFaceId);
    }
  }

  std::vector<std::vector<long> > polygonFaces;
  std::vector<long> polygonFaceLabels;
  std::vector<std::vector<long> > polygonCellFaces;
  mergeTriangleFans(pointList, globalFaces, globalFaceLabels, cellFaceIds,
                    removablePoint, boundarySplitEdges,
                    polygonFaces, polygonFaceLabels,
                    polygonCellFaces);
  globalFaceLabels.swap(polygonFaceLabels);
  cellFaceIds.swap(polygonCellFaces);

  // Interior primal vertices are used to assemble the barycentric pieces but
  // do not lie on the final dual-cell surfaces. Remove such orphan points and
  // remap the global face connectivity before returning the mesh.
  std::vector<long> oldToNew(pointList.size(), -1);
  for (std::vector<std::vector<long> >::const_iterator face = polygonFaces.begin();
       face != polygonFaces.end(); ++face)
    for (std::vector<long>::const_iterator vertex = face->begin(); vertex != face->end(); ++vertex)
      oldToNew[*vertex] = 0;

  std::vector<Point> compactPoints;
  compactPoints.reserve(pointList.size());
  for (long oldId = 0; oldId < static_cast<long>(pointList.size()); ++oldId) {
    if (oldToNew[oldId] < 0)
      continue;
    oldToNew[oldId] = static_cast<long>(compactPoints.size());
    compactPoints.push_back(pointList[oldId]);
  }
  for (std::vector<std::vector<long> >::iterator face = polygonFaces.begin();
       face != polygonFaces.end(); ++face)
    for (std::vector<long>::iterator vertex = face->begin(); vertex != face->end(); ++vertex)
      *vertex = oldToNew[*vertex];
  pointList.swap(compactPoints);

  nodes->resize(static_cast<long>(pointList.size()), 3);
  for (long i = 0; i < static_cast<long>(pointList.size()); ++i) {
    (*nodes)(i, 0L) = pointList[i].x;
    (*nodes)(i, 1L) = pointList[i].y;
    (*nodes)(i, 2L) = pointList[i].z;
  }

  faces->resize(static_cast<long>(polygonFaces.size()));
  for (long i = 0; i < static_cast<long>(polygonFaces.size()); ++i) {
    (*faces)[i].resize(static_cast<long>(polygonFaces[i].size()));
    for (long j = 0; j < static_cast<long>(polygonFaces[i].size()); ++j)
      (*faces)[i][j] = polygonFaces[i][j];
  }

  cells->resize(Th.nv);
  for (long i = 0; i < Th.nv; ++i) {
    (*cells)[i].resize(static_cast<long>(cellFaceIds[i].size()));
    for (long j = 0; j < static_cast<long>(cellFaceIds[i].size()); ++j)
      (*cells)[i][j] = cellFaceIds[i][j];
  }

  if (labels) {
    labels->resize(Th.nv);
    for (long i = 0; i < Th.nv; ++i)
      (*labels)[i] = cellRegion[i];
  }
  if (faceLabels) {
    faceLabels->resize(static_cast<long>(globalFaceLabels.size()));
    for (long i = 0; i < static_cast<long>(globalFaceLabels.size()); ++i)
      (*faceLabels)[i] = globalFaceLabels[i];
  }

  if (verbosity) {
    std::cout << "PDMT 3D " << mode << ": " << Th.nt << " tetrahedra -> " << Th.nv
              << " polyhedra, " << polygonFaces.size() << " polygonal faces and "
              << pointList.size() << " nodes" << std::endl;
    if (!conservedPrimalEdges.empty())
      std::cout << "PDMT 3D dual: conserved " << conservedPrimalEdges.size()
                << " primal boundary edges from physical group(s) "
                << conserveEdge << std::endl;
  }

  return static_cast<long>(Th.nv);
}

#endif
