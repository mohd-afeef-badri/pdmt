/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

     Build a conforming barycentric dual of a triangular surface mesh
     embedded in three dimensions.

*****************************************************************************/

#ifndef PDMT_DUAL_MESH_SURFACE_3D_HPP
#define PDMT_DUAL_MESH_SURFACE_3D_HPP

#include "dualMesh3D.hpp"

namespace Pdmt3S {

using Pdmt3D::Point;
using Pdmt3D::EdgeKey;

struct SurfaceEdge {
  std::vector<long> triangles;
  bool feature;
  long midpoint;
  SurfaceEdge() : feature(false), midpoint(-1) {}
};

inline std::vector<long> mergedWedgeBoundary(
    long primalVertex,
    const std::vector<long> &component,
    const std::vector<std::array<long, 3> > &triangles,
    const std::map<EdgeKey, SurfaceEdge> &edges,
    const std::vector<long> &triangleCentres) {
  std::map<EdgeKey, long> edgeCount;
  for (std::vector<long>::const_iterator triangle = component.begin();
       triangle != component.end(); ++triangle) {
    const std::array<long, 3> &vertices = triangles[*triangle];
    int local = -1;
    for (int i = 0; i < 3; ++i)
      if (vertices[i] == primalVertex)
        local = i;
    if (local < 0)
      ExecError("PdmtBuildDual3S: inconsistent vertex-to-triangle incidence");
    const long next = vertices[(local + 1) % 3];
    const long previous = vertices[(local + 2) % 3];
    const long wedge[4] = {
        primalVertex,
        edges.find(Pdmt3D::edgeKey(primalVertex, next))->second.midpoint,
        triangleCentres[*triangle],
        edges.find(Pdmt3D::edgeKey(primalVertex, previous))->second.midpoint};
    for (int i = 0; i < 4; ++i)
      ++edgeCount[Pdmt3D::edgeKey(wedge[i], wedge[(i + 1) % 4])];
  }

  std::map<long, std::vector<long> > adjacency;
  for (std::map<EdgeKey, long>::const_iterator edge = edgeCount.begin();
       edge != edgeCount.end(); ++edge) {
    if (edge->second == 1) {
      adjacency[edge->first[0]].push_back(edge->first[1]);
      adjacency[edge->first[1]].push_back(edge->first[0]);
    } else if (edge->second != 2) {
      ExecError("PdmtBuildDual3S: non-manifold wedge boundary");
    }
  }
  if (adjacency.size() < 3)
    ExecError("PdmtBuildDual3S: a dual polygon has fewer than three vertices");
  for (std::map<long, std::vector<long> >::const_iterator vertex = adjacency.begin();
       vertex != adjacency.end(); ++vertex)
    if (vertex->second.size() != 2)
      ExecError("PdmtBuildDual3S: a dual polygon does not have one simple boundary");

  std::vector<long> loop;
  const long start = adjacency.begin()->first;
  long previous = -1;
  long current = start;
  do {
    loop.push_back(current);
    const std::vector<long> &neighbours = adjacency[current];
    const long next = neighbours[0] == previous ? neighbours[1] : neighbours[0];
    previous = current;
    current = next;
    if (loop.size() > adjacency.size())
      ExecError("PdmtBuildDual3S: failed to order a dual polygon");
  } while (current != start);
  if (loop.size() != adjacency.size())
    ExecError("PdmtBuildDual3S: a dual polygon contains multiple boundary loops");
  return loop;
}

} // namespace Pdmt3S

class pdmtBuildDual3S_Op : public E_F0mps {
public:
  Expression mesh;

  static const int n_name_param = 8;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];

  pdmtBuildDual3S_Op(const basicAC_F0 &args, Expression inputMesh)
      : mesh(inputMesh) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator()(Stack stack) const;
};

basicAC_F0::name_and_type pdmtBuildDual3S_Op::name_param[] = {
    {"nodes", &typeid(KNM<double> *)},
    {"cells", &typeid(KN<KN<long> > *)},
    {"labels", &typeid(KN<long> *)},
    {"featureAngle", &typeid(double)},
    {"meshFile", &typeid(std::string *)},
    {"conserveEdge", &typeid(std::string *)},
    {"mode", &typeid(std::string *)},
    {"medMeshName", &typeid(std::string *)}};

class pdmtBuildDual3S : public OneOperator {
public:
  pdmtBuildDual3S() : OneOperator(atype<long>(), atype<pmeshS>()) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new pdmtBuildDual3S_Op(args, t[0]->CastTo(args[0]));
  }
};

AnyType pdmtBuildDual3S_Op::operator()(Stack stack) const {
  using namespace Pdmt3D;
  using namespace Pdmt3S;

  const pmeshS pTh = GetAny<pmeshS>((*mesh)(stack));
  if (!pTh)
    ExecError("PdmtBuildDual3S: null input mesh");
  if (!nargs[0] || !nargs[1])
    ExecError("PdmtBuildDual3S: nodes and cells output arrays are required");

  KNM<double> *nodes = GetAny<KNM<double> *>((*nargs[0])(stack));
  KN<KN<long> > *cells = GetAny<KN<KN<long> > *>((*nargs[1])(stack));
  KN<long> *labels = nargs[2] ? GetAny<KN<long> *>((*nargs[2])(stack)) : 0;
  const double featureAngle = nargs[3] ? GetAny<double>((*nargs[3])(stack)) : 45.0;
  std::string meshFile;
  std::string conserveEdge;
  std::string mode = "subdivided_dual";
  std::string medMeshName;
  if (nargs[4])
    meshFile = *GetAny<std::string *>((*nargs[4])(stack));
  if (nargs[5])
    conserveEdge = *GetAny<std::string *>((*nargs[5])(stack));
  if (nargs[6])
    mode = *GetAny<std::string *>((*nargs[6])(stack));
  if (nargs[7])
    medMeshName = *GetAny<std::string *>((*nargs[7])(stack));

  const MeshS &Th = *pTh;
  if (Th.nt == 0 || Th.nv == 0)
    ExecError("PdmtBuildDual3S: the triangular surface mesh is empty");
  if (featureAngle < 0.0 || featureAngle > 180.0)
    ExecError("PdmtBuildDual3S: featureAngle must be between 0 and 180 degrees");
  if (mode != "subdivided_dual" && mode != "smooth_dual")
    ExecError("PdmtBuildDual3S: mode must be subdivided_dual or smooth_dual");
  const bool smoothDual = mode == "smooth_dual";

  std::vector<Point> pointList;
  pointList.reserve(Th.nv + 3 * Th.nt + Th.nt);
  for (long vertex = 0; vertex < Th.nv; ++vertex) {
    const Point point = {Th(vertex).x, Th(vertex).y, Th(vertex).z};
    pointList.push_back(point);
  }

  std::vector<std::array<long, 3> > triangles(Th.nt);
  std::vector<Point> triangleNormals(Th.nt);
  std::vector<long> triangleCentres(Th.nt);
  std::vector<std::vector<long> > incidentTriangles(Th.nv);
  std::map<EdgeKey, SurfaceEdge> surfaceEdges;
  for (long triangle = 0; triangle < Th.nt; ++triangle) {
    for (int i = 0; i < 3; ++i) {
      triangles[triangle][i] = Th(Th[triangle][i]);
      incidentTriangles[triangles[triangle][i]].push_back(triangle);
    }
    const Point &a = pointList[triangles[triangle][0]];
    const Point &b = pointList[triangles[triangle][1]];
    const Point &c = pointList[triangles[triangle][2]];
    triangleNormals[triangle] = cross(Pdmt3D::minus(b, a), Pdmt3D::minus(c, a));
    if (norm(triangleNormals[triangle]) == 0.0)
      ExecError("PdmtBuildDual3S: degenerate input triangle");
    const Point centre = {(a.x + b.x + c.x) / 3.0,
                          (a.y + b.y + c.y) / 3.0,
                          (a.z + b.z + c.z) / 3.0};
    triangleCentres[triangle] = static_cast<long>(pointList.size());
    pointList.push_back(centre);
    for (int i = 0; i < 3; ++i)
      surfaceEdges[edgeKey(triangles[triangle][i],
                           triangles[triangle][(i + 1) % 3])].triangles.push_back(triangle);
  }

  std::set<EdgeKey> allEdges;
  for (std::map<EdgeKey, SurfaceEdge>::const_iterator edge = surfaceEdges.begin();
       edge != surfaceEdges.end(); ++edge)
    allEdges.insert(edge->first);
  std::set<EdgeKey> conservedEdges;
  if (!splitPhysicalNames(conserveEdge).empty()) {
    if (meshFile.empty())
      ExecError("PdmtBuildDual3S: meshFile is required when conserveEdge is set");
    if (meshFile.find(".med") != std::string::npos) {
#ifdef MEDCOUPLING
      conservedEdges = readMedFeatureEdges(
          meshFile, medMeshName, conserveEdge, Th, allEdges, -1,
          "PdmtBuildDual3S");
#else
      ExecError("PdmtBuildDual3S: MED conserveEdge requires MEDCoupling support");
#endif
    } else {
      conservedEdges = readGmshFeatureEdges(meshFile, conserveEdge, Th, allEdges);
    }
  }

  const double pi = 4.0 * std::atan(1.0);
  const double featureCos = std::cos(featureAngle * pi / 180.0);
  for (std::map<EdgeKey, SurfaceEdge>::iterator edge = surfaceEdges.begin();
       edge != surfaceEdges.end(); ++edge) {
    if (edge->second.triangles.empty() || edge->second.triangles.size() > 2)
      ExecError("PdmtBuildDual3S: the input surface is non-manifold at an edge");
    const Point &a = pointList[edge->first[0]];
    const Point &b = pointList[edge->first[1]];
    const Point midpoint = {(a.x + b.x) * 0.5,
                            (a.y + b.y) * 0.5,
                            (a.z + b.z) * 0.5};
    edge->second.midpoint = static_cast<long>(pointList.size());
    pointList.push_back(midpoint);

    edge->second.feature = edge->second.triangles.size() != 2 ||
        conservedEdges.find(edge->first) != conservedEdges.end();
    if (!edge->second.feature) {
      const long first = edge->second.triangles[0];
      const long second = edge->second.triangles[1];
      edge->second.feature = Th[first].lab != Th[second].lab;
      if (!edge->second.feature) {
        const double cosine = dot(triangleNormals[first], triangleNormals[second]) /
            (norm(triangleNormals[first]) * norm(triangleNormals[second]));
        edge->second.feature = cosine < featureCos;
      }
    }
  }

  std::vector<std::vector<long> > triangleAdjacency(Th.nt);
  std::set<long> smoothMidpointNodes;
  for (std::map<EdgeKey, SurfaceEdge>::const_iterator edge = surfaceEdges.begin();
       edge != surfaceEdges.end(); ++edge) {
    if (!edge->second.feature)
      smoothMidpointNodes.insert(edge->second.midpoint);
    if (!edge->second.feature) {
      const long first = edge->second.triangles[0];
      const long second = edge->second.triangles[1];
      triangleAdjacency[first].push_back(second);
      triangleAdjacency[second].push_back(first);
    }
  }

  std::vector<std::vector<long> > polygonList;
  std::vector<long> polygonLabels;
  for (long vertex = 0; vertex < Th.nv; ++vertex) {
    std::set<long> remaining(incidentTriangles[vertex].begin(),
                             incidentTriangles[vertex].end());
    while (!remaining.empty()) {
      const long seed = *remaining.begin();
      remaining.erase(seed);
      std::vector<long> component;
      std::vector<long> work(1, seed);
      while (!work.empty()) {
        const long triangle = work.back();
        work.pop_back();
        component.push_back(triangle);
        for (std::vector<long>::const_iterator adjacent = triangleAdjacency[triangle].begin();
             adjacent != triangleAdjacency[triangle].end(); ++adjacent) {
          std::set<long>::iterator found = remaining.find(*adjacent);
          if (found != remaining.end()) {
            remaining.erase(found);
            work.push_back(*adjacent);
          }
        }
      }

      std::vector<long> polygon = mergedWedgeBoundary(
          vertex, component, triangles, surfaceEdges, triangleCentres);
      if (smoothDual) {
        std::vector<long> simplified;
        for (std::vector<long>::const_iterator point = polygon.begin();
             point != polygon.end(); ++point)
          if (smoothMidpointNodes.find(*point) == smoothMidpointNodes.end())
            simplified.push_back(*point);
        if (simplified.size() < 3)
          ExecError("PdmtBuildDual3S: smooth_dual produced a polygon with fewer than three vertices");
        polygon.swap(simplified);
      }
      Point desiredNormal = {0.0, 0.0, 0.0};
      for (std::vector<long>::const_iterator triangle = component.begin();
           triangle != component.end(); ++triangle)
        desiredNormal = add(desiredNormal, triangleNormals[*triangle]);
      if (dot(polygonAreaVector(polygon, pointList), desiredNormal) < 0.0)
        std::reverse(polygon.begin(), polygon.end());
      polygonList.push_back(polygon);
      polygonLabels.push_back(Th[component[0]].lab);
    }
  }

  std::vector<long> oldToNew(pointList.size(), -1);
  for (std::vector<std::vector<long> >::const_iterator polygon = polygonList.begin();
       polygon != polygonList.end(); ++polygon)
    for (std::vector<long>::const_iterator vertex = polygon->begin();
         vertex != polygon->end(); ++vertex)
      oldToNew[*vertex] = 0;
  std::vector<Point> compactPoints;
  for (long oldId = 0; oldId < static_cast<long>(pointList.size()); ++oldId)
    if (oldToNew[oldId] >= 0) {
      oldToNew[oldId] = static_cast<long>(compactPoints.size());
      compactPoints.push_back(pointList[oldId]);
    }
  for (std::vector<std::vector<long> >::iterator polygon = polygonList.begin();
       polygon != polygonList.end(); ++polygon)
    for (std::vector<long>::iterator vertex = polygon->begin();
         vertex != polygon->end(); ++vertex)
      *vertex = oldToNew[*vertex];

  nodes->resize(static_cast<long>(compactPoints.size()), 3);
  for (long node = 0; node < static_cast<long>(compactPoints.size()); ++node) {
    (*nodes)(node, 0L) = compactPoints[node].x;
    (*nodes)(node, 1L) = compactPoints[node].y;
    (*nodes)(node, 2L) = compactPoints[node].z;
  }
  cells->resize(static_cast<long>(polygonList.size()));
  for (long cell = 0; cell < static_cast<long>(polygonList.size()); ++cell) {
    (*cells)[cell].resize(static_cast<long>(polygonList[cell].size()));
    for (long vertex = 0; vertex < static_cast<long>(polygonList[cell].size()); ++vertex)
      (*cells)[cell][vertex] = polygonList[cell][vertex];
  }
  if (labels) {
    labels->resize(static_cast<long>(polygonLabels.size()));
    for (long cell = 0; cell < static_cast<long>(polygonLabels.size()); ++cell)
      (*labels)[cell] = polygonLabels[cell];
  }

  if (verbosity) {
    std::cout << "PDMT 3S " << mode << ": " << Th.nt << " surface triangles -> "
              << polygonList.size() << " polygons and " << compactPoints.size()
              << " nodes" << std::endl;
    if (!conservedEdges.empty())
      std::cout << "PDMT 3S dual: conserved " << conservedEdges.size()
                << " primal surface edges from physical group(s) "
                << conserveEdge << std::endl;
  }
  return static_cast<long>(polygonList.size());
}

#endif
