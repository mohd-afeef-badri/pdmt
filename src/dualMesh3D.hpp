/*****************************************************************************

         This file is a part of PDMT (Parallel Dual Meshing Tool)

             Build a conforming dual of a tetrahedral mesh.

*****************************************************************************/

#ifndef PDMT_DUAL_MESH_3D_HPP
#define PDMT_DUAL_MESH_3D_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
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

inline bool finitePoint(const Point &point) {
  return std::isfinite(point.x) && std::isfinite(point.y) &&
         std::isfinite(point.z);
}

// The circumcentre of a triangle embedded in 3D.  It lies in the triangle's
// plane and is equidistant from all three vertices.
inline Point triangleCircumcenter(const Point &a, const Point &b,
                                  const Point &c) {
  const Point u = minus(b, a);
  const Point v = minus(c, a);
  const Point normal = cross(u, v);
  const double denominator = 2.0 * dot(normal, normal);
  if (denominator == 0.0)
    ExecError("PdmtBuildDual3D: cannot construct a circumcentric dual from a degenerate triangle");

  const Point offset = scale(
      add(scale(cross(v, normal), dot(u, u)),
          scale(cross(normal, u), dot(v, v))),
      1.0 / denominator);
  const Point center = add(a, offset);
  if (!finitePoint(center))
    ExecError("PdmtBuildDual3D: triangle circumcentre is not finite");
  return center;
}

// Solve the three perpendicular-bisector equations relative to a.  All
// circumcentres of tetrahedra incident to one primal edge consequently lie
// in that edge's single bisector plane, making the merged dual face planar.
inline Point tetraCircumcenter(const Point &a, const Point &b,
                               const Point &c, const Point &d) {
  const Point u = minus(b, a);
  const Point v = minus(c, a);
  const Point w = minus(d, a);
  const double denominator = 2.0 * dot(u, cross(v, w));
  if (denominator == 0.0)
    ExecError("PdmtBuildDual3D: cannot construct a circumcentric dual from a degenerate tetrahedron");

  const Point offset = scale(
      add(add(scale(cross(v, w), dot(u, u)),
              scale(cross(w, u), dot(v, v))),
          scale(cross(u, v), dot(w, w))),
      1.0 / denominator);
  const Point center = add(a, offset);
  if (!finitePoint(center))
    ExecError("PdmtBuildDual3D: tetrahedron circumcentre is not finite");
  return center;
}

inline double signedTetVolume6(const Point &a, const Point &b,
                               const Point &c, const Point &d) {
  return dot(cross(minus(b, a), minus(c, a)), minus(d, a));
}

inline Point weightedSimplexCenter(
    const long *vertices, int count,
    const std::vector<Point> &points,
    const std::vector<double> &cellPreference) {
  Point center = {0.0, 0.0, 0.0};
  double weightSum = 0.0;
  for (int local = 0; local < count; ++local) {
    // A cell asking for more volume must push the common dual node away from
    // its seed, hence the inverse preference in the affine combination.
    const double weight = 1.0 / cellPreference[vertices[local]];
    center = add(center, scale(points[vertices[local]], weight));
    weightSum += weight;
  }
  return scale(center, 1.0 / weightSum);
}

template<class MeshType>
inline std::vector<double> dualCellVolumes(
    const MeshType &mesh,
    const std::vector<Point> &points,
    const std::vector<double> &cellPreference) {
  std::vector<double> volumes(mesh.nv, 0.0);
  for (long tet = 0; tet < mesh.nt; ++tet) {
    long vertex[4];
    for (int local = 0; local < 4; ++local)
      vertex[local] = mesh(mesh[tet][local]);
    const Point tetCenter =
        weightedSimplexCenter(vertex, 4, points, cellPreference);
    for (int vi = 0; vi < 4; ++vi) {
      const long v = vertex[vi];
      for (int ui = 0; ui < 4; ++ui) {
        if (ui == vi)
          continue;
        const long edgeVertex[2] = {v, vertex[ui]};
        const Point edgeCenter =
            weightedSimplexCenter(edgeVertex, 2, points, cellPreference);
        for (int wi = 0; wi < 4; ++wi) {
          if (wi == vi || wi == ui)
            continue;
          const long faceVertex[3] = {v, vertex[ui], vertex[wi]};
          const Point faceCenter =
              weightedSimplexCenter(faceVertex, 3, points, cellPreference);
          volumes[v] += std::abs(signedTetVolume6(
              points[v], edgeCenter, faceCenter, tetCenter)) / 6.0;
        }
      }
    }
  }
  return volumes;
}

inline double coefficientOfVariation(const std::vector<double> &values) {
  if (values.empty())
    return 0.0;
  double mean = 0.0;
  for (std::vector<double>::const_iterator value = values.begin();
       value != values.end(); ++value)
    mean += *value;
  mean /= values.size();
  if (mean <= 0.0)
    return 0.0;
  double variance = 0.0;
  for (std::vector<double>::const_iterator value = values.begin();
       value != values.end(); ++value) {
    const double difference = *value - mean;
    variance += difference * difference;
  }
  return std::sqrt(variance / values.size()) / mean;
}

inline double boundaryInteriorVolumeRatio(
    const std::vector<double> &volumes,
    const std::vector<char> &boundaryVertex) {
  double boundarySum = 0.0;
  double interiorSum = 0.0;
  long boundaryCount = 0;
  long interiorCount = 0;
  for (long vertex = 0; vertex < static_cast<long>(volumes.size()); ++vertex) {
    if (boundaryVertex[vertex]) {
      boundarySum += volumes[vertex];
      ++boundaryCount;
    } else {
      interiorSum += volumes[vertex];
      ++interiorCount;
    }
  }
  if (!boundaryCount || !interiorCount || interiorSum <= 0.0)
    return 0.0;
  return (boundarySum / boundaryCount) / (interiorSum / interiorCount);
}

template<class MeshType>
inline long regularizeDualVolumes(
    const MeshType &mesh,
    const std::vector<Point> &points,
    long requestedIterations,
    double relaxation,
    std::vector<double> &cellPreference,
    double &initialCv,
    double &finalCv,
    double &initialBoundaryRatio,
    double &finalBoundaryRatio) {
  cellPreference.assign(mesh.nv, 1.0);
  initialCv = finalCv = 0.0;
  initialBoundaryRatio = finalBoundaryRatio = 0.0;
  if (requestedIterations <= 0)
    return 0;

  std::vector<std::set<long> > neighbourSets(mesh.nv);
  std::vector<char> boundaryVertex(mesh.nv, 0);
  for (long tet = 0; tet < mesh.nt; ++tet) {
    long vertex[4];
    for (int local = 0; local < 4; ++local)
      vertex[local] = mesh(mesh[tet][local]);
    for (int first = 0; first < 4; ++first)
      for (int second = first + 1; second < 4; ++second) {
        neighbourSets[vertex[first]].insert(vertex[second]);
        neighbourSets[vertex[second]].insert(vertex[first]);
      }
  }
  for (long face = 0; face < mesh.nbe; ++face)
    for (int local = 0; local < 3; ++local)
      boundaryVertex[mesh(mesh.be(face)[local])] = 1;

  std::vector<double> volumes =
      dualCellVolumes(mesh, points, cellPreference);
  initialCv = coefficientOfVariation(volumes);
  initialBoundaryRatio =
      boundaryInteriorVolumeRatio(volumes, boundaryVertex);
  long completedIterations = 0;
  for (long iteration = 0; iteration < requestedIterations; ++iteration) {
    std::vector<double> updatedPreference(cellPreference);
    double maximumLogUpdate = 0.0;
    for (long vertex = 0; vertex < mesh.nv; ++vertex) {
      if (volumes[vertex] <= 0.0 || neighbourSets[vertex].empty())
        continue;
      double localTarget = 0.0;
      long targetCount = 0;
      // Boundary cells have a truncated vertex star. Compare them directly
      // with adjacent interior cells instead of averaging mostly with other
      // small boundary cells. Interior targets are smoothed within the
      // interior population so that the two updates do not cancel each other.
      for (std::set<long>::const_iterator neighbour =
               neighbourSets[vertex].begin();
           neighbour != neighbourSets[vertex].end(); ++neighbour) {
        if (!boundaryVertex[*neighbour]) {
          localTarget += volumes[*neighbour];
          ++targetCount;
        }
      }
      if (!boundaryVertex[vertex]) {
        localTarget += volumes[vertex];
        ++targetCount;
      }
      if (!targetCount) {
        localTarget = volumes[vertex];
        targetCount = 1;
        for (std::set<long>::const_iterator neighbour =
                 neighbourSets[vertex].begin();
             neighbour != neighbourSets[vertex].end(); ++neighbour) {
          localTarget += volumes[*neighbour];
          ++targetCount;
        }
      }
      localTarget /= targetCount;
      const double logUpdate =
          relaxation * std::log(localTarget / volumes[vertex]);
      maximumLogUpdate = std::max(maximumLogUpdate, std::abs(logUpdate));
      updatedPreference[vertex] =
          cellPreference[vertex] * std::exp(logUpdate);
    }

    // Only relative preferences matter. Normalization also prevents drift
    // and keeps all weighted centers comfortably inside their simplices.
    double meanLogPreference = 0.0;
    for (long vertex = 0; vertex < mesh.nv; ++vertex)
      meanLogPreference += std::log(updatedPreference[vertex]);
    meanLogPreference /= mesh.nv;
    for (long vertex = 0; vertex < mesh.nv; ++vertex)
      cellPreference[vertex] = std::max(
          0.05, std::min(20.0,
              std::exp(std::log(updatedPreference[vertex]) -
                       meanLogPreference)));

    volumes = dualCellVolumes(mesh, points, cellPreference);
    ++completedIterations;
    if (maximumLogUpdate <= 1.e-10)
      break;
  }
  finalCv = coefficientOfVariation(volumes);
  finalBoundaryRatio =
      boundaryInteriorVolumeRatio(volumes, boundaryVertex);
  return completedIterations;
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

inline double squaredDistance(const Point &a, const Point &b) {
  const Point difference = minus(a, b);
  return dot(difference, difference);
}

inline double pointCoordinate(const Point &point, int axis) {
  if (axis == 0)
    return point.x;
  if (axis == 1)
    return point.y;
  return point.z;
}

class PointKdTree {
  struct Node {
    long point;
    long left;
    long right;
    int axis;
  };

  const std::vector<Point> &points_;
  std::vector<long> order_;
  std::vector<Node> nodes_;
  long root_;

  struct CoordinateLess {
    const std::vector<Point> &points;
    int axis;
    CoordinateLess(const std::vector<Point> &inputPoints, int inputAxis)
        : points(inputPoints), axis(inputAxis) {}
    bool operator()(long first, long second) const {
      return pointCoordinate(points[first], axis) <
             pointCoordinate(points[second], axis);
    }
  };

  long build(long begin, long end, int depth) {
    if (begin >= end)
      return -1;
    const int axis = depth % 3;
    const long middle = begin + (end - begin) / 2;
    std::nth_element(order_.begin() + begin, order_.begin() + middle,
                     order_.begin() + end,
                     CoordinateLess(points_, axis));
    const long node = static_cast<long>(nodes_.size());
    Node entry = {order_[middle], -1, -1, axis};
    nodes_.push_back(entry);
    const long left = build(begin, middle, depth + 1);
    const long right = build(middle + 1, end, depth + 1);
    nodes_[node].left = left;
    nodes_[node].right = right;
    return node;
  }

  void nearest(long node, const Point &query, long &bestPoint,
               double &bestSquaredDistance) const {
    if (node < 0)
      return;
    const Node &entry = nodes_[node];
    const double candidateSquaredDistance =
        squaredDistance(query, points_[entry.point]);
    if (candidateSquaredDistance < bestSquaredDistance) {
      bestSquaredDistance = candidateSquaredDistance;
      bestPoint = entry.point;
    }
    const double difference =
        pointCoordinate(query, entry.axis) -
        pointCoordinate(points_[entry.point], entry.axis);
    const long nearChild = difference <= 0.0 ? entry.left : entry.right;
    const long farChild = difference <= 0.0 ? entry.right : entry.left;
    nearest(nearChild, query, bestPoint, bestSquaredDistance);
    if (difference * difference <= bestSquaredDistance)
      nearest(farChild, query, bestPoint, bestSquaredDistance);
  }

public:
  explicit PointKdTree(const std::vector<Point> &points)
      : points_(points), order_(points.size()), root_(-1) {
    for (long point = 0; point < static_cast<long>(points.size()); ++point)
      order_[point] = point;
    nodes_.reserve(points.size());
    root_ = build(0, static_cast<long>(order_.size()), 0);
  }

  long nearest(const Point &query, double &bestSquaredDistance) const {
    long bestPoint = -1;
    bestSquaredDistance = std::numeric_limits<double>::infinity();
    nearest(root_, query, bestPoint, bestSquaredDistance);
    return bestPoint;
  }
};

inline double polygonLengthScale(const std::vector<long> &polygon,
                                 const std::vector<Point> &points) {
  if (polygon.empty())
    return 0.0;
  Point minimum = points[polygon[0]];
  Point maximum = minimum;
  for (std::vector<long>::const_iterator vertex = polygon.begin() + 1;
       vertex != polygon.end(); ++vertex) {
    const Point &point = points[*vertex];
    minimum.x = std::min(minimum.x, point.x);
    minimum.y = std::min(minimum.y, point.y);
    minimum.z = std::min(minimum.z, point.z);
    maximum.x = std::max(maximum.x, point.x);
    maximum.y = std::max(maximum.y, point.y);
    maximum.z = std::max(maximum.z, point.z);
  }
  return norm(minus(maximum, minimum));
}

// A right primal triangle has its circumcentre on an edge midpoint.  The two
// points have different topological ids, but retaining both creates a
// zero-length polygon edge that MEDCoupling reports as overlapping.  Collapse
// only consecutive coincident points; this does not change the face boundary.
inline void removeConsecutiveCoincidentPoints(
    std::vector<long> &polygon, const std::vector<Point> &points) {
  if (polygon.size() <= 1)
    return;
  const double scale = polygonLengthScale(polygon, points);
  const double tolerance = std::max(
      scale * 1.e-12, 64.0 * std::numeric_limits<double>::epsilon());
  const double toleranceSquared = tolerance * tolerance;
  bool changed = true;
  while (changed && polygon.size() > 1) {
    changed = false;
    for (long vertex = 0; vertex < static_cast<long>(polygon.size()); ++vertex) {
      const long next = (vertex + 1) % polygon.size();
      if (vertex != next &&
          squaredDistance(points[polygon[vertex]], points[polygon[next]]) <=
          toleranceSquared) {
        polygon.erase(polygon.begin() + next);
        changed = true;
        break;
      }
    }
  }
}

struct Point2D {
  double x, y;
};

inline Point2D projectPoint(const Point &point, int droppedCoordinate) {
  Point2D projected;
  if (droppedCoordinate == 0) {
    projected.x = point.y;
    projected.y = point.z;
  } else if (droppedCoordinate == 1) {
    projected.x = point.x;
    projected.y = point.z;
  } else {
    projected.x = point.x;
    projected.y = point.y;
  }
  return projected;
}

inline double orientation2D(const Point2D &a, const Point2D &b,
                            const Point2D &c) {
  return (b.x - a.x) * (c.y - a.y) -
         (b.y - a.y) * (c.x - a.x);
}

inline bool pointOnSegment2D(const Point2D &point, const Point2D &a,
                             const Point2D &b, double lengthTolerance,
                             double areaTolerance) {
  if (std::abs(orientation2D(a, b, point)) > areaTolerance)
    return false;
  return point.x >= std::min(a.x, b.x) - lengthTolerance &&
         point.x <= std::max(a.x, b.x) + lengthTolerance &&
         point.y >= std::min(a.y, b.y) - lengthTolerance &&
         point.y <= std::max(a.y, b.y) + lengthTolerance;
}

inline int orientationSign(double value, double tolerance) {
  if (value > tolerance)
    return 1;
  if (value < -tolerance)
    return -1;
  return 0;
}

// Circumcentres may lie outside non-well-centred tetrahedra.  The topological
// boundary of a merged triangle fan can then cross itself even though all its
// vertices are coplanar.  Detect those rings before emitting a POLYGON cell.
inline bool polygonHasSelfIntersection(
    const std::vector<long> &polygon, const std::vector<Point> &points,
    const Point &planeNormal) {
  if (polygon.size() <= 3)
    return false;

  const double scale = polygonLengthScale(polygon, points);
  if (scale == 0.0)
    return true;
  const double lengthTolerance = scale * 1.e-12;
  const double areaTolerance = scale * scale * 1.e-12;

  const double normalComponent[3] = {
      std::abs(planeNormal.x), std::abs(planeNormal.y),
      std::abs(planeNormal.z)};
  int droppedCoordinate = 0;
  if (normalComponent[1] > normalComponent[droppedCoordinate])
    droppedCoordinate = 1;
  if (normalComponent[2] > normalComponent[droppedCoordinate])
    droppedCoordinate = 2;
  if (normalComponent[droppedCoordinate] <= areaTolerance)
    return true;

  std::vector<Point2D> projected;
  projected.reserve(polygon.size());
  for (std::vector<long>::const_iterator vertex = polygon.begin();
       vertex != polygon.end(); ++vertex)
    projected.push_back(projectPoint(points[*vertex], droppedCoordinate));

  const long edgeCount = static_cast<long>(projected.size());
  for (long first = 0; first < edgeCount; ++first) {
    const long firstNext = (first + 1) % edgeCount;
    const Point2D &a = projected[first];
    const Point2D &b = projected[firstNext];
    const double edgeDx = b.x - a.x;
    const double edgeDy = b.y - a.y;
    if (edgeDx * edgeDx + edgeDy * edgeDy <=
        lengthTolerance * lengthTolerance)
      return true;

    // Adjacent collinear edges which reverse direction overlap in their
    // interiors.  Ordinary straight-through collinear vertices are valid.
    const long previous = (first + edgeCount - 1) % edgeCount;
    const Point2D &p = projected[previous];
    if (std::abs(orientation2D(p, a, b)) <= areaTolerance) {
      const double incomingX = a.x - p.x;
      const double incomingY = a.y - p.y;
      if (incomingX * edgeDx + incomingY * edgeDy <
          -lengthTolerance * lengthTolerance)
        return true;
    }

    for (long second = first + 1; second < edgeCount; ++second) {
      const long secondNext = (second + 1) % edgeCount;
      if (second == firstNext || secondNext == first)
        continue;
      const Point2D &c = projected[second];
      const Point2D &d = projected[secondNext];
      const double abc = orientation2D(a, b, c);
      const double abd = orientation2D(a, b, d);
      const double cda = orientation2D(c, d, a);
      const double cdb = orientation2D(c, d, b);
      const int abcSign = orientationSign(abc, areaTolerance);
      const int abdSign = orientationSign(abd, areaTolerance);
      const int cdaSign = orientationSign(cda, areaTolerance);
      const int cdbSign = orientationSign(cdb, areaTolerance);
      if ((abcSign * abdSign < 0 && cdaSign * cdbSign < 0) ||
          (abcSign == 0 && pointOnSegment2D(c, a, b, lengthTolerance,
                                            areaTolerance)) ||
          (abdSign == 0 && pointOnSegment2D(d, a, b, lengthTolerance,
                                            areaTolerance)) ||
          (cdaSign == 0 && pointOnSegment2D(a, c, d, lengthTolerance,
                                            areaTolerance)) ||
          (cdbSign == 0 && pointOnSegment2D(b, c, d, lengthTolerance,
                                            areaTolerance)))
        return true;
    }
  }
  return false;
}

inline bool polygonIsPlanar(const std::vector<long> &polygon,
                            const std::vector<Point> &points,
                            const Point &planeNormal) {
  if (polygon.size() <= 3)
    return true;
  const double normalLength = norm(planeNormal);
  const double scale = polygonLengthScale(polygon, points);
  if (normalLength == 0.0 || scale == 0.0)
    return false;
  const Point &origin = points[polygon[0]];
  // Boundary-patch vertices can be assembled through different tetrahedra.
  // Their absolute plane residual stays near roundoff, but normalizing by a
  // very small clipped face can amplify it.  This tolerance still limits the
  // normalized warp to 1e-8; visibly warped or crossed polygons remain
  // rejected, and internal bisector faces are normally many orders tighter.
  const double tolerance = scale * 1.e-8;
  for (std::vector<long>::const_iterator vertex = polygon.begin() + 1;
       vertex != polygon.end(); ++vertex)
    if (std::abs(dot(minus(points[*vertex], origin), planeNormal)) /
            normalLength >
        tolerance)
      return false;
  return true;
}

inline Point polygonPlaneNormal(const std::vector<long> &polygon,
                                const std::vector<Point> &points) {
  Point best = {0.0, 0.0, 0.0};
  double bestSquaredLength = 0.0;
  if (polygon.size() < 3)
    return best;
  const Point &origin = points[polygon[0]];
  for (long first = 1; first + 1 < static_cast<long>(polygon.size()); ++first)
    for (long second = first + 1;
         second < static_cast<long>(polygon.size()); ++second) {
      const Point candidate = cross(
          minus(points[polygon[first]], origin),
          minus(points[polygon[second]], origin));
      const double squaredLength = dot(candidate, candidate);
      if (squaredLength > bestSquaredLength) {
        bestSquaredLength = squaredLength;
        best = candidate;
      }
    }
  return best;
}

struct RestrictedFace {
  std::vector<long> vertices;
  // kind == 0: a face inherited from the primal tetrahedron.
  // kind == 1: a cap on the bisector with neighbourSeed.
  // kind == 2: an artificial bounding-box face.
  int kind;
  FaceKey primalFace;
  long neighbourSeed;
};

struct RestrictedPolyhedron {
  std::vector<Point> points;
  std::vector<RestrictedFace> faces;
};

inline RestrictedPolyhedron makeRestrictedBoundingBox(
    const Point &minimum, const Point &maximum) {
  RestrictedPolyhedron result;
  const Point corners[8] = {
      {minimum.x, minimum.y, minimum.z},
      {maximum.x, minimum.y, minimum.z},
      {maximum.x, maximum.y, minimum.z},
      {minimum.x, maximum.y, minimum.z},
      {minimum.x, minimum.y, maximum.z},
      {maximum.x, minimum.y, maximum.z},
      {maximum.x, maximum.y, maximum.z},
      {minimum.x, maximum.y, maximum.z}};
  for (int corner = 0; corner < 8; ++corner)
    result.points.push_back(corners[corner]);
  static const int boxFaces[6][4] = {
      {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
      {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  for (int side = 0; side < 6; ++side) {
    RestrictedFace face;
    face.kind = 2;
    face.primalFace = FaceKey{{-1, -1, -1}};
    face.neighbourSeed = -1;
    for (int vertex = 0; vertex < 4; ++vertex)
      face.vertices.push_back(boxFaces[side][vertex]);
    result.faces.push_back(face);
  }
  return result;
}

inline void removeRepeatedPolygonVertices(std::vector<long> &polygon) {
  if (polygon.empty())
    return;
  std::vector<long> cleaned;
  cleaned.reserve(polygon.size());
  for (std::vector<long>::const_iterator vertex = polygon.begin();
       vertex != polygon.end(); ++vertex)
    if (cleaned.empty() || cleaned.back() != *vertex)
      cleaned.push_back(*vertex);
  if (cleaned.size() > 1 && cleaned.front() == cleaned.back())
    cleaned.pop_back();
  polygon.swap(cleaned);
}

inline void orderCapPolygon(std::vector<long> &polygon,
                            const std::vector<Point> &points,
                            const Point &outwardNormal) {
  Point centroid = {0.0, 0.0, 0.0};
  for (std::vector<long>::const_iterator vertex = polygon.begin();
       vertex != polygon.end(); ++vertex)
    centroid = add(centroid, points[*vertex]);
  centroid = scale(centroid, 1.0 / polygon.size());

  Point firstAxis = minus(points[polygon[0]], centroid);
  double firstLength = norm(firstAxis);
  for (long vertex = 1;
       firstLength == 0.0 && vertex < static_cast<long>(polygon.size());
       ++vertex) {
    firstAxis = minus(points[polygon[vertex]], centroid);
    firstLength = norm(firstAxis);
  }
  if (firstLength == 0.0)
    return;
  firstAxis = scale(firstAxis, 1.0 / firstLength);
  Point secondAxis = cross(outwardNormal, firstAxis);
  const double secondLength = norm(secondAxis);
  if (secondLength == 0.0)
    return;
  secondAxis = scale(secondAxis, 1.0 / secondLength);

  std::vector<std::pair<double, long> > angular;
  angular.reserve(polygon.size());
  for (std::vector<long>::const_iterator vertex = polygon.begin();
       vertex != polygon.end(); ++vertex) {
    const Point radius = minus(points[*vertex], centroid);
    angular.push_back(std::make_pair(
        std::atan2(dot(radius, secondAxis), dot(radius, firstAxis)),
        *vertex));
  }
  std::sort(angular.begin(), angular.end());
  for (long vertex = 0; vertex < static_cast<long>(polygon.size()); ++vertex)
    polygon[vertex] = angular[vertex].second;
  if (dot(polygonAreaVector(polygon, points), outwardNormal) < 0.0)
    std::reverse(polygon.begin(), polygon.end());
}

inline void clipRestrictedPolyhedron(RestrictedPolyhedron &polyhedron,
                                     const Point &normal,
                                     double offset,
                                     long neighbourSeed,
                                     double tolerance,
                                     const FaceKey *primalFace = 0) {
  std::vector<double> signedDistance(polyhedron.points.size());
  for (long vertex = 0;
       vertex < static_cast<long>(polyhedron.points.size()); ++vertex)
    signedDistance[vertex] =
        dot(polyhedron.points[vertex], normal) - offset;

  std::map<EdgeKey, long> intersectionVertex;
  std::vector<long> capVertices;
  std::vector<RestrictedFace> clippedFaces;
  for (std::vector<RestrictedFace>::const_iterator oldFace =
           polyhedron.faces.begin();
       oldFace != polyhedron.faces.end(); ++oldFace) {
    RestrictedFace clipped = *oldFace;
    clipped.vertices.clear();
    const long count = static_cast<long>(oldFace->vertices.size());
    for (long edge = 0; edge < count; ++edge) {
      const long a = oldFace->vertices[edge];
      const long b = oldFace->vertices[(edge + 1) % count];
      const bool aInside = signedDistance[a] <= tolerance;
      const bool bInside = signedDistance[b] <= tolerance;
      if (aInside)
        clipped.vertices.push_back(a);
      if (aInside == bInside)
        continue;

      const EdgeKey key = edgeKey(a, b);
      std::map<EdgeKey, long>::const_iterator found =
          intersectionVertex.find(key);
      long intersection;
      if (found != intersectionVertex.end()) {
        intersection = found->second;
      } else {
        const double denominator = signedDistance[a] - signedDistance[b];
        if (denominator == 0.0)
          continue;
        const double parameter = signedDistance[a] / denominator;
        const Point point = add(
            polyhedron.points[a],
            scale(minus(polyhedron.points[b], polyhedron.points[a]),
                  parameter));
        intersection = static_cast<long>(polyhedron.points.size());
        polyhedron.points.push_back(point);
        signedDistance.push_back(0.0);
        intersectionVertex[key] = intersection;
      }
      clipped.vertices.push_back(intersection);
      capVertices.push_back(intersection);
    }
    removeRepeatedPolygonVertices(clipped.vertices);
    if (clipped.vertices.size() >= 3)
      clippedFaces.push_back(clipped);
  }

  std::sort(capVertices.begin(), capVertices.end());
  capVertices.erase(std::unique(capVertices.begin(), capVertices.end()),
                    capVertices.end());
  if (capVertices.size() >= 3) {
    orderCapPolygon(capVertices, polyhedron.points, normal);
    RestrictedFace cap;
    cap.vertices.swap(capVertices);
    if (primalFace) {
      cap.kind = 0;
      cap.primalFace = *primalFace;
      cap.neighbourSeed = -1;
    } else {
      cap.kind = 1;
      cap.primalFace = FaceKey{{-1, -1, -1}};
      cap.neighbourSeed = neighbourSeed;
    }
    clippedFaces.push_back(cap);
  }
  polyhedron.faces.swap(clippedFaces);
}

inline bool restrictedPolyhedronHasVolume(
    const RestrictedPolyhedron &polyhedron, double distanceTolerance) {
  std::set<long> usedVertices;
  for (std::vector<RestrictedFace>::const_iterator face =
           polyhedron.faces.begin();
       face != polyhedron.faces.end(); ++face)
    usedVertices.insert(face->vertices.begin(), face->vertices.end());
  if (usedVertices.size() < 4)
    return false;

  for (std::vector<RestrictedFace>::const_iterator face =
           polyhedron.faces.begin();
       face != polyhedron.faces.end(); ++face) {
    if (face->vertices.size() < 3)
      continue;
    const Point normal = cross(
        minus(polyhedron.points[face->vertices[1]],
              polyhedron.points[face->vertices[0]]),
        minus(polyhedron.points[face->vertices[2]],
              polyhedron.points[face->vertices[0]]));
    const double normalLength = norm(normal);
    if (normalLength == 0.0)
      continue;
    const Point &origin = polyhedron.points[face->vertices[0]];
    for (std::set<long>::const_iterator vertex = usedVertices.begin();
         vertex != usedVertices.end(); ++vertex)
      if (std::abs(dot(minus(polyhedron.points[*vertex], origin), normal)) >
          distanceTolerance * normalLength)
        return true;
  }
  return false;
}

typedef std::array<long long, 3> QuantizedPointKey;

inline long findOrAddRestrictedPoint(
    const Point &point, double tolerance,
    std::vector<Point> &points,
    std::map<QuantizedPointKey, std::vector<long> > &buckets) {
  QuantizedPointKey centre = {{
      static_cast<long long>(std::floor(point.x / tolerance)),
      static_cast<long long>(std::floor(point.y / tolerance)),
      static_cast<long long>(std::floor(point.z / tolerance))}};
  const double toleranceSquared = tolerance * tolerance;
  for (int dx = -1; dx <= 1; ++dx)
    for (int dy = -1; dy <= 1; ++dy)
      for (int dz = -1; dz <= 1; ++dz) {
        QuantizedPointKey key = {{
            centre[0] + dx, centre[1] + dy, centre[2] + dz}};
        std::map<QuantizedPointKey, std::vector<long> >::const_iterator bucket =
            buckets.find(key);
        if (bucket == buckets.end())
          continue;
        for (std::vector<long>::const_iterator candidate =
                 bucket->second.begin();
             candidate != bucket->second.end(); ++candidate)
          if (squaredDistance(points[*candidate], point) <= toleranceSquared)
            return *candidate;
      }
  const long id = static_cast<long>(points.size());
  points.push_back(point);
  buckets[centre].push_back(id);
  return id;
}

inline bool pointLiesOnSegment(const Point &point, const Point &a,
                               const Point &b, double tolerance) {
  const Point direction = minus(b, a);
  const double lengthSquared = dot(direction, direction);
  if (lengthSquared == 0.0)
    return false;
  const double parameter = dot(minus(point, a), direction) / lengthSquared;
  if (parameter < -tolerance || parameter > 1.0 + tolerance)
    return false;
  const Point projection = add(a, scale(direction, parameter));
  return squaredDistance(point, projection) <=
         tolerance * tolerance * lengthSquared;
}

inline void removeCollinearPolygonVertices(
    std::vector<long> &polygon, const std::vector<Point> &points,
    double tolerance) {
  bool changed = true;
  while (changed && polygon.size() > 3) {
    changed = false;
    std::vector<long> simplified;
    simplified.reserve(polygon.size());
    const long count = static_cast<long>(polygon.size());
    for (long vertex = 0; vertex < count; ++vertex) {
      const long previous = polygon[(vertex + count - 1) % count];
      const long current = polygon[vertex];
      const long next = polygon[(vertex + 1) % count];
      if (pointLiesOnSegment(
              points[current], points[previous], points[next], tolerance)) {
        changed = true;
        continue;
      }
      simplified.push_back(current);
    }
    if (simplified.size() < 3)
      break;
    polygon.swap(simplified);
  }
}

inline bool trianglesAreCoplanar(
    const std::array<long, 3> &first,
    const std::array<long, 3> &second,
    const std::vector<Point> &points) {
  const Point firstNormal = cross(
      minus(points[first[1]], points[first[0]]),
      minus(points[first[2]], points[first[0]]));
  const Point secondNormal = cross(
      minus(points[second[1]], points[second[0]]),
      minus(points[second[2]], points[second[0]]));
  const double denominator = norm(firstNormal) * norm(secondNormal);
  if (denominator == 0.0)
    return false;
  const double cosine = std::abs(dot(firstNormal, secondNormal)) / denominator;
  if (cosine < 1.0 - 1.e-10)
    return false;

  double scale = 0.0;
  for (int firstVertex = 0; firstVertex < 3; ++firstVertex)
    for (int secondVertex = firstVertex + 1; secondVertex < 3;
         ++secondVertex) {
      scale = std::max(
          scale, norm(minus(points[first[firstVertex]],
                            points[first[secondVertex]])));
      scale = std::max(
          scale, norm(minus(points[second[firstVertex]],
                            points[second[secondVertex]])));
    }
  if (scale == 0.0)
    return false;
  const double planeDistance = std::abs(dot(
      firstNormal, minus(points[second[0]], points[first[0]]))) /
      norm(firstNormal);
  return planeDistance <= scale * 1.e-9;
}

inline std::vector<std::vector<long> > connectedTriangleComponents(
    const std::vector<long> &faceIds,
    const std::vector<std::array<long, 3> > &faces,
    const std::vector<Point> &points,
    const std::set<EdgeKey> &splitEdges,
    const std::set<EdgeKey> &forcedMergeEdges,
    bool splitNonCoplanar,
    bool onlyForcedMerges) {
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
        if (splitNonCoplanar &&
            forcedMergeEdges.find(edge->first) == forcedMergeEdges.end()) {
          if (onlyForcedMerges ||
              !trianglesAreCoplanar(faces[faceIds[incident[i]]],
                                    faces[faceIds[incident[j]]], points))
            continue;
        }
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
    const std::vector<std::array<long, 3> > &faces,
    long groupType, long groupFirst, long groupSecond) {
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
      std::ostringstream message;
      message << "PdmtBuildDual3D: non-manifold triangle fan while merging "
              << (groupType ? "boundary cell/patch " : "cell interface ")
              << groupFirst << "/" << groupSecond << "; edge "
              << edge->first[0] << "-" << edge->first[1] << " has "
              << edge->second << " incident triangles";
      ExecError(message.str());
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

inline long disjointSetRoot(std::vector<long> &parent, long item) {
  long root = item;
  while (parent[root] != root)
    root = parent[root];
  while (parent[item] != item) {
    const long next = parent[item];
    parent[item] = root;
    item = next;
  }
  return root;
}

inline void conformPolygonCellEdges(
    const std::vector<Point> &points,
    std::vector<std::vector<long> > &polygons,
    const std::vector<std::vector<long> > &polygonCells) {
  for (int iteration = 0; iteration < 8; ++iteration) {
    std::map<EdgeKey, std::set<long> > edgeSplits;
    for (std::vector<std::vector<long> >::const_iterator cell =
             polygonCells.begin();
         cell != polygonCells.end(); ++cell) {
      std::set<long> cellVertices;
      for (std::vector<long>::const_iterator encodedFace = cell->begin();
           encodedFace != cell->end(); ++encodedFace) {
        const std::vector<long> &face =
            polygons[std::labs(*encodedFace) - 1];
        cellVertices.insert(face.begin(), face.end());
      }
      for (std::vector<long>::const_iterator encodedFace = cell->begin();
           encodedFace != cell->end(); ++encodedFace) {
        const std::vector<long> &face =
            polygons[std::labs(*encodedFace) - 1];
        for (long edge = 0; edge < static_cast<long>(face.size()); ++edge) {
          const long a = face[edge];
          const long b = face[(edge + 1) % face.size()];
          for (std::set<long>::const_iterator candidate =
                   cellVertices.begin();
               candidate != cellVertices.end(); ++candidate) {
            if (*candidate == a || *candidate == b)
              continue;
            if (pointLiesOnSegment(
                    points[*candidate], points[a], points[b], 1.e-10))
              edgeSplits[edgeKey(a, b)].insert(*candidate);
          }
        }
      }
    }
    if (edgeSplits.empty())
      return;

    bool inserted = false;
    for (std::vector<std::vector<long> >::iterator face = polygons.begin();
         face != polygons.end(); ++face) {
      std::vector<long> conformed;
      for (long edge = 0; edge < static_cast<long>(face->size()); ++edge) {
        const long a = (*face)[edge];
        const long b = (*face)[(edge + 1) % face->size()];
        conformed.push_back(a);
        const std::map<EdgeKey, std::set<long> >::const_iterator splits =
            edgeSplits.find(edgeKey(a, b));
        if (splits == edgeSplits.end())
          continue;
        const Point direction = minus(points[b], points[a]);
        const double lengthSquared = dot(direction, direction);
        std::vector<std::pair<double, long> > ordered;
        for (std::set<long>::const_iterator candidate =
                 splits->second.begin();
             candidate != splits->second.end(); ++candidate) {
          const double parameter =
              dot(minus(points[*candidate], points[a]), direction) /
              lengthSquared;
          if (parameter > 1.e-10 && parameter < 1.0 - 1.e-10)
            ordered.push_back(std::make_pair(parameter, *candidate));
        }
        if (ordered.empty())
          continue;
        std::sort(ordered.begin(), ordered.end());
        for (std::vector<std::pair<double, long> >::const_iterator split =
                 ordered.begin();
             split != ordered.end(); ++split)
          conformed.push_back(split->second);
        inserted = true;
      }
      face->swap(conformed);
    }
    if (!inserted)
      return;
  }
  ExecError("PdmtBuildDual3D: failed to conform polygon edges within a Voronoi cell");
}

inline void mergeTriangleFans(
    const std::vector<Point> &points,
    const std::vector<std::array<long, 3> > &triangles,
    const std::vector<long> &triangleLabels,
    const std::vector<long> &trianglePatches,
    const std::vector<std::vector<long> > &triangleCells,
    const std::vector<char> &removablePoint,
    const std::set<EdgeKey> &boundarySplitEdges,
    const std::set<EdgeKey> &forcedMergeEdges,
    bool splitNonCoplanar,
    bool onlyForcedMerges,
    bool validateCircumcentricFans,
    std::vector<std::vector<long> > &polygons,
    std::vector<long> &polygonLabels,
    std::vector<std::vector<long> > &polygonCells) {
  typedef std::array<long, 3> GroupKey;
  if (triangleLabels.size() != triangles.size() ||
      trianglePatches.size() != triangles.size())
    ExecError("PdmtBuildDual3D: invalid triangle face metadata");
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
      const GroupKey key = {{1, uses[face][0].first, trianglePatches[face]}};
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
        group->second, triangles, points,
        group->first[0] ? boundarySplitEdges : noSplitEdges,
        forcedMergeEdges,
        splitNonCoplanar, onlyForcedMerges);
    for (std::vector<std::vector<long> >::const_iterator component = components.begin();
         component != components.end(); ++component) {
      std::vector<long> loop = componentBoundaryLoop(
          *component, triangles, group->first[0],
          group->first[1], group->first[2]);
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

      if (validateCircumcentricFans)
        removeConsecutiveCoincidentPoints(loop, points);
      if (validateCircumcentricFans)
        removeCollinearPolygonVertices(loop, points, 1.e-10);
      if (loop.size() < 3)
        continue;

      const Point geometricNormal = polygonPlaneNormal(loop, points);
      if (validateCircumcentricFans && norm(geometricNormal) == 0.0)
        continue;
      const bool invalidPlanarity = validateCircumcentricFans &&
          !polygonIsPlanar(loop, points, geometricNormal);
      const bool invalidSimplicity = validateCircumcentricFans &&
          polygonHasSelfIntersection(loop, points, geometricNormal);
      if (invalidPlanarity || invalidSimplicity) {
        std::ostringstream message;
        message << "PdmtBuildDual3D: the circumcentric "
                << (group->first[0] ? "domain-boundary" : "cell-interface")
                << " face for cell " << group->first[1];
        if (!group->first[0])
          message << " and cell " << group->first[2];
        message << " is not a "
                << (invalidPlanarity ? "planar" : "simple")
                << " Voronoi polygon (" << loop.size() << " vertices)";
        ExecError(message.str());
      }

      const long newFace = static_cast<long>(polygons.size());
      polygons.push_back(loop);
      polygonLabels.push_back(
          group->first[0] ? triangleLabels[component->front()] : 0);
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

  static const int n_name_param = 12;
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
    {"medMeshName", &typeid(std::string *)},
    {"smoothIterations", &typeid(long)},
    {"smoothRelaxation", &typeid(double)}};

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
  const long smoothIterations =
      nargs[10] ? GetAny<long>((*nargs[10])(stack)) : 0;
  const double smoothRelaxation =
      nargs[11] ? GetAny<double>((*nargs[11])(stack)) : 0.3;
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
  if (mode != "subdivided_dual" && mode != "smooth_dual" &&
      mode != "circumcentric_dual")
    ExecError("PdmtBuildDual3D: mode must be subdivided_dual, smooth_dual, or circumcentric_dual");
  if (smoothIterations < 0)
    ExecError("PdmtBuildDual3D: smoothIterations must be non-negative");
  if (smoothRelaxation <= 0.0 || smoothRelaxation > 1.0)
    ExecError("PdmtBuildDual3D: smoothRelaxation must be in (0,1]");
  const bool smoothDual = mode == "smooth_dual";
  const bool circumcentricDual = mode == "circumcentric_dual";
  if (circumcentricDual && smoothIterations > 0)
    ExecError("PdmtBuildDual3D: smoothIterations is incompatible with circumcentric_dual because volume regularization destroys face planarity");
  const bool simplifyDual = smoothDual || circumcentricDual;

  std::set<EdgeKey> primalEdges;
  std::set<FaceKey> primalFaces;
  std::map<FaceKey, std::vector<std::pair<long, long> > > primalFaceUses;
  std::vector<std::vector<long> > incidentTets(Th.nv);
  std::vector<std::set<long> > primalNeighbours(Th.nv);
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
    for (int i = 0; i < 6; ++i) {
      const long a = vertex[tetEdges[i][0]];
      const long b = vertex[tetEdges[i][1]];
      primalEdges.insert(edgeKey(a, b));
      primalNeighbours[a].insert(b);
      primalNeighbours[b].insert(a);
    }
    for (int i = 0; i < 4; ++i) {
      const FaceKey face = faceKey(
          vertex[tetFaces[i][0]], vertex[tetFaces[i][1]],
          vertex[tetFaces[i][2]]);
      primalFaces.insert(face);
      primalFaceUses[face].push_back(std::make_pair(t, vertex[i]));
    }
  }

  std::vector<std::vector<long> > adjacentTets(Th.nt);
  for (std::map<FaceKey, std::vector<std::pair<long, long> > >::const_iterator
           face = primalFaceUses.begin();
       face != primalFaceUses.end(); ++face) {
    if (face->second.size() == 2) {
      const long first = face->second[0].first;
      const long second = face->second[1].first;
      adjacentTets[first].push_back(second);
      adjacentTets[second].push_back(first);
    }
  }

  std::set<FaceKey> boundaryPrimalFaces;
  for (long boundary = 0; boundary < Th.nbe; ++boundary) {
    const long a = Th(Th.be(boundary)[0]);
    const long b = Th(Th.be(boundary)[1]);
    const long c = Th(Th.be(boundary)[2]);
    boundaryPrimalFaces.insert(faceKey(a, b, c));
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

  if (circumcentricDual) {
    std::map<FaceKey, long> boundaryLabel;
    std::vector<FaceKey> boundaryFacesById;
    std::map<FaceKey, long> boundaryFaceId;
    std::map<EdgeKey, std::vector<long> > boundaryFacesByEdge;
    for (long boundary = 0; boundary < Th.nbe; ++boundary) {
      const FaceKey face = faceKey(
          Th(Th.be(boundary)[0]), Th(Th.be(boundary)[1]),
          Th(Th.be(boundary)[2]));
      boundaryLabel[face] = Th.be(boundary).lab;
      if (boundaryFaceId.find(face) == boundaryFaceId.end()) {
        const long id = static_cast<long>(boundaryFacesById.size());
        boundaryFaceId[face] = id;
        boundaryFacesById.push_back(face);
        for (int edge = 0; edge < 3; ++edge)
          boundaryFacesByEdge[edgeKey(
              face[edge], face[(edge + 1) % 3])].push_back(id);
      }
    }

    std::vector<long> boundaryParent(boundaryFacesById.size());
    for (long face = 0;
         face < static_cast<long>(boundaryParent.size()); ++face)
      boundaryParent[face] = face;
    for (std::map<EdgeKey, std::vector<long> >::const_iterator edge =
             boundaryFacesByEdge.begin();
         edge != boundaryFacesByEdge.end(); ++edge) {
      if (conservedPrimalEdges.find(edge->first) !=
          conservedPrimalEdges.end())
        continue;
      const std::vector<long> &incident = edge->second;
      for (long first = 0; first < static_cast<long>(incident.size());
           ++first)
        for (long second = first + 1;
             second < static_cast<long>(incident.size()); ++second) {
          const FaceKey &firstFace = boundaryFacesById[incident[first]];
          const FaceKey &secondFace = boundaryFacesById[incident[second]];
          if (boundaryLabel[firstFace] != boundaryLabel[secondFace] ||
              !trianglesAreCoplanar(
                  firstFace, secondFace, pointList))
            continue;
          const long firstRoot =
              disjointSetRoot(boundaryParent, incident[first]);
          const long secondRoot =
              disjointSetRoot(boundaryParent, incident[second]);
          if (firstRoot != secondRoot)
            boundaryParent[secondRoot] = firstRoot;
        }
    }
    std::map<long, long> rootPatch;
    std::map<FaceKey, long> boundaryPatch;
    for (long face = 0;
         face < static_cast<long>(boundaryFacesById.size()); ++face) {
      const long root = disjointSetRoot(boundaryParent, face);
      if (rootPatch.find(root) == rootPatch.end())
        rootPatch[root] = static_cast<long>(rootPatch.size());
      boundaryPatch[boundaryFacesById[face]] = rootPatch[root];
    }

    Point minimum = pointList[0];
    Point maximum = pointList[0];
    for (std::vector<Point>::const_iterator point = pointList.begin() + 1;
         point != pointList.end(); ++point) {
      minimum.x = std::min(minimum.x, point->x);
      minimum.y = std::min(minimum.y, point->y);
      minimum.z = std::min(minimum.z, point->z);
      maximum.x = std::max(maximum.x, point->x);
      maximum.y = std::max(maximum.y, point->y);
      maximum.z = std::max(maximum.z, point->z);
    }
    const double domainScale = norm(Pdmt3D::minus(maximum, minimum));
    if (domainScale == 0.0)
      ExecError("PdmtBuildDual3D: the tetrahedral mesh has zero extent");
    const double pointTolerance = std::max(
        domainScale * 1.e-12,
        64.0 * std::numeric_limits<double>::epsilon());
    const double distanceTolerance = domainScale * 1.e-12;

    std::vector<Point> restrictedPoints;
    std::map<QuantizedPointKey, std::vector<long> > pointBuckets;
    std::vector<std::array<long, 3> > restrictedTriangles;
    std::vector<long> restrictedTriangleLabels;
    std::vector<long> restrictedTrianglePatches;
    std::vector<std::vector<long> > restrictedCellFaces(Th.nv);
    std::set<std::array<long, 5> > restrictedInterfaceTriangles;

    Point boxMinimum = {
        minimum.x - domainScale, minimum.y - domainScale,
        minimum.z - domainScale};
    Point boxMaximum = {
        maximum.x + domainScale, maximum.y + domainScale,
        maximum.z + domainScale};
    const double constraintTolerance =
        domainScale * domainScale * 1.e-11;
    const PointKdTree siteTree(pointList);

    // A Voronoi cell of an obtuse Delaunay mesh is not necessarily contained
    // in the tetrahedral star of its seed.  Construct the complete convex cell
    // first, then walk through every connected tetrahedron it intersects.
    // This avoids leaving an interior primal face as a triangular hole when a
    // circumcentre lies outside its tetrahedron.
    std::vector<long> visitedTet(Th.nt, -1);
    for (long seed = 0; seed < Th.nv; ++seed) {
      RestrictedPolyhedron globalCell =
          makeRestrictedBoundingBox(boxMinimum, boxMaximum);
      std::set<long> appliedConstraints;
      for (std::set<long>::const_iterator otherIt =
               primalNeighbours[seed].begin();
           otherIt != primalNeighbours[seed].end(); ++otherIt) {
        const long other = *otherIt;
        appliedConstraints.insert(other);
        const Point normal =
            Pdmt3D::minus(pointList[other], pointList[seed]);
        const double offset = 0.5 *
            (dot(pointList[other], pointList[other]) -
             dot(pointList[seed], pointList[seed]));
        clipRestrictedPolyhedron(
            globalCell, normal, offset, other,
            norm(normal) * distanceTolerance);
      }

      // A non-Delaunay tetrahedral topology can omit sites whose bisectors
      // actually support this geometric Voronoi cell.  A convex polyhedron
      // satisfies every nearest-site half-space exactly when all its vertices
      // do.  Query those vertices in a kd-tree and add the most violated
      // missing constraint until none remains.  This keeps the construction
      // local without an O(number-of-sites^2) all-pairs clipping pass.
      while (true) {
        std::set<long> usedVertices;
        for (std::vector<RestrictedFace>::const_iterator face =
                 globalCell.faces.begin();
             face != globalCell.faces.end(); ++face)
          usedVertices.insert(face->vertices.begin(), face->vertices.end());
        long missingSite = -1;
        double maximumViolation = constraintTolerance;
        for (std::set<long>::const_iterator vertex = usedVertices.begin();
             vertex != usedVertices.end(); ++vertex) {
          const Point &candidate = globalCell.points[*vertex];
          double nearestSquaredDistance = 0.0;
          const long nearestSite =
              siteTree.nearest(candidate, nearestSquaredDistance);
          if (nearestSite == seed ||
              appliedConstraints.find(nearestSite) !=
                  appliedConstraints.end())
            continue;
          const double violation =
              squaredDistance(candidate, pointList[seed]) -
              nearestSquaredDistance;
          if (violation > maximumViolation) {
            maximumViolation = violation;
            missingSite = nearestSite;
          }
        }
        if (missingSite < 0)
          break;
        appliedConstraints.insert(missingSite);
        const Point normal =
            Pdmt3D::minus(pointList[missingSite], pointList[seed]);
        const double offset = 0.5 *
            (dot(pointList[missingSite], pointList[missingSite]) -
             dot(pointList[seed], pointList[seed]));
        clipRestrictedPolyhedron(
            globalCell, normal, offset, missingSite,
            norm(normal) * distanceTolerance);
      }
      if (!restrictedPolyhedronHasVolume(globalCell, distanceTolerance))
        ExecError("PdmtBuildDual3D: a primal vertex has an empty Voronoi cell");

      std::vector<long> work = incidentTets[seed];
      while (!work.empty()) {
        const long tet = work.back();
        work.pop_back();
        if (visitedTet[tet] == seed)
          continue;
        visitedTet[tet] = seed;

        long vertex[4];
        for (int local = 0; local < 4; ++local)
          vertex[local] = Th(Th[tet][local]);
        RestrictedPolyhedron clipped = globalCell;
        for (int opposite = 0; opposite < 4; ++opposite) {
          const FaceKey primalFace = faceKey(
              vertex[tetFaces[opposite][0]],
              vertex[tetFaces[opposite][1]],
              vertex[tetFaces[opposite][2]]);
          const Point &a = pointList[vertex[tetFaces[opposite][0]]];
          const Point &b = pointList[vertex[tetFaces[opposite][1]]];
          const Point &c = pointList[vertex[tetFaces[opposite][2]]];
          Point normal = cross(Pdmt3D::minus(b, a),
                               Pdmt3D::minus(c, a));
          if (dot(normal,
                  Pdmt3D::minus(pointList[vertex[opposite]], a)) > 0.0)
            normal = scale(normal, -1.0);
          clipRestrictedPolyhedron(
              clipped, normal, dot(a, normal), -1,
              norm(normal) * distanceTolerance, &primalFace);
          if (clipped.faces.empty())
            break;
        }
        if (!restrictedPolyhedronHasVolume(clipped, distanceTolerance))
          continue;

        for (std::vector<long>::const_iterator adjacent =
                 adjacentTets[tet].begin();
             adjacent != adjacentTets[tet].end(); ++adjacent)
          if (visitedTet[*adjacent] != seed)
            work.push_back(*adjacent);

        for (std::vector<RestrictedFace>::const_iterator face =
                 clipped.faces.begin();
             face != clipped.faces.end(); ++face) {
          bool boundaryFace = false;
          long label = 0;
          long neighbour = -1;
          if (face->kind == 0) {
            std::map<FaceKey, long>::const_iterator boundary =
                boundaryLabel.find(face->primalFace);
            if (boundary == boundaryLabel.end())
              continue;
            boundaryFace = true;
            label = boundary->second;
          } else if (face->kind == 1) {
            neighbour = face->neighbourSeed;
            if (seed > neighbour)
              continue;
          } else {
            ExecError("PdmtBuildDual3D: the artificial Voronoi bounding box intersects the domain");
          }

          std::vector<long> globalVertices;
          globalVertices.reserve(face->vertices.size());
          Point faceCentroid = {0.0, 0.0, 0.0};
          for (std::vector<long>::const_iterator localVertex =
                   face->vertices.begin();
               localVertex != face->vertices.end(); ++localVertex) {
            const Point &point = clipped.points[*localVertex];
            globalVertices.push_back(findOrAddRestrictedPoint(
                point, pointTolerance, restrictedPoints, pointBuckets));
            faceCentroid = add(faceCentroid, point);
          }
          removeRepeatedPolygonVertices(globalVertices);
          if (globalVertices.size() < 3)
            continue;
          faceCentroid = scale(
              faceCentroid, 1.0 / face->vertices.size());
          const long centroid = findOrAddRestrictedPoint(
              faceCentroid, pointTolerance, restrictedPoints, pointBuckets);

          for (long edge = 0;
               edge < static_cast<long>(globalVertices.size()); ++edge) {
            const long a = globalVertices[edge];
            const long b = globalVertices[
                (edge + 1) % globalVertices.size()];
            if (a == b || a == centroid || b == centroid)
              continue;
            const Point triangleNormal = cross(
                Pdmt3D::minus(
                    restrictedPoints[a], restrictedPoints[centroid]),
                Pdmt3D::minus(
                    restrictedPoints[b], restrictedPoints[centroid]));
            if (norm(triangleNormal) <=
                domainScale * domainScale * 1.e-24)
              continue;
            if (!boundaryFace) {
              std::array<long, 3> triangleVertices = {{centroid, a, b}};
              std::sort(triangleVertices.begin(), triangleVertices.end());
              const std::array<long, 5> key = {{
                  std::min(seed, neighbour), std::max(seed, neighbour),
                  triangleVertices[0], triangleVertices[1],
                  triangleVertices[2]}};
              if (!restrictedInterfaceTriangles.insert(key).second)
                continue;
            }
            const long triangle =
                static_cast<long>(restrictedTriangles.size());
            restrictedTriangles.push_back(
                std::array<long, 3>{{centroid, a, b}});
            restrictedTriangleLabels.push_back(label);
            restrictedTrianglePatches.push_back(
                boundaryFace ? boundaryPatch[face->primalFace] : 0);
            restrictedCellFaces[seed].push_back(triangle + 1);
            if (!boundaryFace)
              restrictedCellFaces[neighbour].push_back(-(triangle + 1));
          }
        }
      }
    }

    std::vector<char> noRemovablePoint(restrictedPoints.size(), 0);
    const std::set<EdgeKey> noBoundarySplitEdges;
    const std::set<EdgeKey> noForcedMergeEdges;
    std::vector<std::vector<long> > polygonFaces;
    std::vector<long> polygonFaceLabels;
    std::vector<std::vector<long> > polygonCellFaces;
    mergeTriangleFans(
        restrictedPoints, restrictedTriangles, restrictedTriangleLabels,
        restrictedTrianglePatches,
        restrictedCellFaces, noRemovablePoint,
        noBoundarySplitEdges, noForcedMergeEdges,
        true, false, true,
        polygonFaces, polygonFaceLabels, polygonCellFaces);
    conformPolygonCellEdges(
        restrictedPoints, polygonFaces, polygonCellFaces);

    std::vector<long> oldToNew(restrictedPoints.size(), -1);
    for (std::vector<std::vector<long> >::const_iterator face =
             polygonFaces.begin();
         face != polygonFaces.end(); ++face)
      for (std::vector<long>::const_iterator vertexId = face->begin();
           vertexId != face->end(); ++vertexId)
        oldToNew[*vertexId] = 0;
    std::vector<Point> compactPoints;
    compactPoints.reserve(restrictedPoints.size());
    for (long old = 0; old < static_cast<long>(restrictedPoints.size()); ++old)
      if (oldToNew[old] >= 0) {
        oldToNew[old] = static_cast<long>(compactPoints.size());
        compactPoints.push_back(restrictedPoints[old]);
      }
    for (std::vector<std::vector<long> >::iterator face =
             polygonFaces.begin();
         face != polygonFaces.end(); ++face)
      for (std::vector<long>::iterator vertexId = face->begin();
           vertexId != face->end(); ++vertexId)
        *vertexId = oldToNew[*vertexId];

    nodes->resize(static_cast<long>(compactPoints.size()), 3);
    for (long point = 0; point < static_cast<long>(compactPoints.size());
         ++point) {
      (*nodes)(point, 0L) = compactPoints[point].x;
      (*nodes)(point, 1L) = compactPoints[point].y;
      (*nodes)(point, 2L) = compactPoints[point].z;
    }
    faces->resize(static_cast<long>(polygonFaces.size()));
    for (long face = 0; face < static_cast<long>(polygonFaces.size()); ++face) {
      (*faces)[face].resize(static_cast<long>(polygonFaces[face].size()));
      for (long vertexId = 0;
           vertexId < static_cast<long>(polygonFaces[face].size()); ++vertexId)
        (*faces)[face][vertexId] = polygonFaces[face][vertexId];
    }
    cells->resize(Th.nv);
    for (long cell = 0; cell < Th.nv; ++cell) {
      (*cells)[cell].resize(static_cast<long>(polygonCellFaces[cell].size()));
      for (long face = 0;
           face < static_cast<long>(polygonCellFaces[cell].size()); ++face)
        (*cells)[cell][face] = polygonCellFaces[cell][face];
    }
    if (labels) {
      labels->resize(Th.nv);
      for (long cell = 0; cell < Th.nv; ++cell)
        (*labels)[cell] = cellRegion[cell];
    }
    if (faceLabels) {
      faceLabels->resize(static_cast<long>(polygonFaceLabels.size()));
      for (long face = 0;
           face < static_cast<long>(polygonFaceLabels.size()); ++face)
        (*faceLabels)[face] = polygonFaceLabels[face];
    }
    if (verbosity)
      std::cout << "PDMT 3D circumcentric_dual: restricted Voronoi clipping "
                << Th.nt << " tetrahedra -> " << Th.nv
                << " polyhedra, " << polygonFaces.size()
                << " polygonal faces and " << compactPoints.size()
                << " nodes" << std::endl;
    return static_cast<long>(Th.nv);
  }

  std::vector<double> cellPreference;
  double initialVolumeCv = 0.0;
  double finalVolumeCv = 0.0;
  double initialBoundaryRatio = 0.0;
  double finalBoundaryRatio = 0.0;
  const long completedSmoothIterations = regularizeDualVolumes(
      Th, pointList, smoothIterations, smoothRelaxation, cellPreference,
      initialVolumeCv, finalVolumeCv, initialBoundaryRatio,
      finalBoundaryRatio);
  if (cellPreference.empty())
    cellPreference.assign(Th.nv, 1.0);

  std::map<EdgeKey, long> edgeNodes;
  for (std::set<EdgeKey>::const_iterator it = primalEdges.begin();
       it != primalEdges.end(); ++it) {
    const long vertex[2] = {(*it)[0], (*it)[1]};
    const Point p = circumcentricDual
        ? scale(add(pointList[vertex[0]], pointList[vertex[1]]), 0.5)
        : weightedSimplexCenter(vertex, 2, pointList, cellPreference);
    edgeNodes[*it] = static_cast<long>(pointList.size());
    pointList.push_back(p);
  }

  std::map<FaceKey, long> faceNodes;
  for (std::set<FaceKey>::const_iterator it = primalFaces.begin();
       it != primalFaces.end(); ++it) {
    const long vertex[3] = {(*it)[0], (*it)[1], (*it)[2]};
    const Point p = circumcentricDual
        ? triangleCircumcenter(pointList[vertex[0]], pointList[vertex[1]],
                               pointList[vertex[2]])
        : weightedSimplexCenter(vertex, 3, pointList, cellPreference);
    faceNodes[*it] = static_cast<long>(pointList.size());
    pointList.push_back(p);
  }

  std::vector<long> tetNodes(Th.nt);
  for (long t = 0; t < Th.nt; ++t) {
    long vertex[4];
    for (int i = 0; i < 4; ++i)
      vertex[i] = Th(Th[t][i]);
    const Point p = circumcentricDual
        ? tetraCircumcenter(pointList[vertex[0]], pointList[vertex[1]],
                            pointList[vertex[2]], pointList[vertex[3]])
        : weightedSimplexCenter(vertex, 4, pointList, cellPreference);
    tetNodes[t] = static_cast<long>(pointList.size());
    pointList.push_back(p);
  }

  std::map<FaceKey, long> boundaryTriangleLabels;
  std::map<EdgeKey, std::vector<long> > boundaryEdgeFaces;
  for (long b = 0; b < Th.nbe; ++b) {
    long v[3] = {Th(Th.be(b)[0]), Th(Th.be(b)[1]), Th(Th.be(b)[2])};
    const FaceKey primalFace = faceKey(v[0], v[1], v[2]);
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
  if (simplifyDual)
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
        const double normalCosine = dot(n0, n1) / denominator;
        // A circumcentric interior face is planar by construction.  On the
        // domain boundary, do not undo that guarantee by merging triangles
        // from distinct planes, regardless of the requested feature angle.
        feature = circumcentricDual
            ? normalCosine < 1.0 - 64.0 * std::numeric_limits<double>::epsilon()
            : normalCosine < featureCos;
      }
    }
    if (feature) {
      const long midpoint = edgeNodes[edge->first];
      boundarySplitEdges.insert(edgeKey(edge->first[0], midpoint));
      boundarySplitEdges.insert(edgeKey(edge->first[1], midpoint));
    } else if (simplifyDual) {
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

  const std::set<EdgeKey> forcedMergeEdges;

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
  mergeTriangleFans(pointList, globalFaces, globalFaceLabels,
                    globalFaceLabels, cellFaceIds,
                    removablePoint, boundarySplitEdges, forcedMergeEdges,
                    false, false,
                    circumcentricDual,
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
    if (smoothIterations > 0)
      std::cout << "PDMT 3D regularization: completed "
                << completedSmoothIterations << "/" << smoothIterations
                << " dual-volume iterations at relaxation "
                << smoothRelaxation
                << "; cell-volume CV " << initialVolumeCv
                << " -> " << finalVolumeCv
                << "; mean boundary/interior volume ratio "
                << initialBoundaryRatio << " -> " << finalBoundaryRatio
                << " (primal boundary and conserved edges unchanged)"
                << std::endl;
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
