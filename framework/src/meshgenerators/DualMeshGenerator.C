//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DualMeshGenerator.h"
#include "Conversion.h"
#include "CastUniquePointer.h"
#include "MooseUtils.h"
#include "libmesh/cell_c0polyhedron.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/face_c0polygon.h"
#include "libmesh/libmesh_exceptions.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/node_elem.h"
#include "libmesh/poly2tri_triangulator.h"
#include "libmesh/type_vector.h"
#include "libmesh/elem.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

registerMooseObject("MooseApp", DualMeshGenerator);

InputParameters
DualMeshGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  MooseEnum dual_mesh_type("voronoi barycentric", "barycentric");
  params.addParam<MooseEnum>("dual_mesh_type",
                             dual_mesh_type,
                             "Whether to output a barycentric dual or a Voronoi dual of the "
                             "Delaunay triangulation of the primal mesh.");
  params.addRangeCheckedParam<Real>(
      "boundary_node_angular_tol",
      1e-8,
      "boundary_node_angular_tol>=0",
      "Angular tolerance, in radians, used to decide whether a primal boundary node is collinear "
      "with its two adjacent boundary edges. Nodes whose boundary angle differs from pi by more "
      "than this tolerance are treated as primal boundary vertices.");
  params.addRangeCheckedParam<Real>(
      "geometry_relative_tol",
      1e-12,
      "geometry_relative_tol>=0",
      "Relative tolerance used for geometric point comparison, intersection, and area checks. The "
      "generator scales this value by the input mesh bounding-box size.");
  params.addClassDescription("Takes a 2D or 3D mesh as input and returns a dual mesh, i.e., "
                             "changes each input node into an element and each input element "
                             "into a node located at its circumcenter or centroid.");
  return params;
}

DualMeshGenerator::DualMeshGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _dual_mesh_type(getParam<MooseEnum>("dual_mesh_type")),
    _boundary_node_angular_tol(getParam<Real>("boundary_node_angular_tol")),
    _geometry_relative_tol(getParam<Real>("geometry_relative_tol"))
{
}

// True only when two elements share a full edge.
static bool
elementsShareTwoNodes(const Elem * a, const Elem * b)
{
  unsigned int shared_nodes = 0;

  for (const auto i : a->node_index_range())
    for (const auto j : b->node_index_range())
      if (a->node_id(i) == b->node_id(j))
        ++shared_nodes;

  mooseAssert(
      (shared_nodes < 3),
      "Detected elements sharing > 3 nodes."); // This is typically a sign of intertwined elements

  return shared_nodes == 2;
}

// Cross product, used for colinearity, segment tests, and concavity tests
static Real
cross2D(const Point & a, const Point & b, const Point & c)
{
  return (b(0) - a(0)) * (c(1) - a(1)) - (b(1) - a(1)) * (c(0) - a(0));
}

static void
addUniquePoint(std::vector<Point> & points, const Point & point, const Real length_tol = 1e-12)
{
  for (const auto & existing_point : points)
    if (MooseUtils::absoluteFuzzyEqual(existing_point, point, length_tol))
      return;

  points.push_back(point);
}

static bool
pointOnSegment2D(const Point & point,
                 const Point & a,
                 const Point & b,
                 const Real length_tol,
                 const Real area_tol)
{
  if (std::abs(cross2D(a, b, point)) > area_tol)
    return false;

  return (point(0) >= std::min(a(0), b(0)) - length_tol &&
          point(0) <= std::max(a(0), b(0)) + length_tol &&
          point(1) >= std::min(a(1), b(1)) - length_tol &&
          point(1) <= std::max(a(1), b(1)) + length_tol);
}

// We use intersections when we clip a mesh back to the primal boundary. Those intersection points
// are necessary to keep the elements being reasonable within the boundary
static void
addSegmentIntersections2D(std::vector<Point> & points,
                          const Point & p0,
                          const Point & p1,
                          const Point & q0,
                          const Point & q1,
                          const Real length_tol,
                          const Real area_tol,
                          const Real parameter_tol)
{
  const Point r = p1 - p0;
  const Point s = q1 - q0;
  const Real denom = r(0) * s(1) - r(1) * s(0);

  if (std::abs(denom) < area_tol)
  {
    if (std::abs(cross2D(p0, p1, q0)) > area_tol)
      return;

    if (pointOnSegment2D(q0, p0, p1, length_tol, area_tol))
      addUniquePoint(points, q0, length_tol);
    if (pointOnSegment2D(q1, p0, p1, length_tol, area_tol))
      addUniquePoint(points, q1, length_tol);
    if (pointOnSegment2D(p0, q0, q1, length_tol, area_tol))
      addUniquePoint(points, p0, length_tol);
    if (pointOnSegment2D(p1, q0, q1, length_tol, area_tol))
      addUniquePoint(points, p1, length_tol);

    return;
  }

  const Point qp = q0 - p0;
  const Real t = (qp(0) * s(1) - qp(1) * s(0)) / denom;
  const Real u = (qp(0) * r(1) - qp(1) * r(0)) / denom;

  if (t >= -parameter_tol && t <= 1.0 + parameter_tol && u >= -parameter_tol &&
      u <= 1.0 + parameter_tol)
    addUniquePoint(points, p0 + t * r, length_tol);
}

// Ray casting test to see if points are inside a polygon
static bool
pointInPolygon2D(const Point & point,
                 const std::vector<Point> & polygon,
                 const Real length_tol,
                 const Real area_tol)
{
  if (polygon.size() < 3)
    return false;

  bool inside = false;

  for (const auto i : index_range(polygon))
  {
    const auto j = (i + polygon.size() - 1) % polygon.size();
    const Point & pi = polygon[i];
    const Point & pj = polygon[j];

    if (pointOnSegment2D(point, pj, pi, length_tol, area_tol))
      return true;

    if ((pi(1) > point(1)) != (pj(1) > point(1)) &&
        point(0) < (pj(0) - pi(0)) * (point(1) - pi(1)) / (pj(1) - pi(1)) + pi(0))
      inside = !inside;
  }

  return inside;
}

static Real
polygonSignedArea2D(const std::vector<Point> & polygon)
{
  if (polygon.size() < 3)
    return 0.0;

  Real area = 0.0;

  for (const auto i : index_range(polygon))
  {
    const auto j = (i + polygon.size() - 1) % polygon.size();
    area += polygon[j](0) * polygon[i](1) - polygon[i](0) * polygon[j](1);
  }

  return 0.5 * area;
}

struct BoundarySegment
{
  dof_id_type node0;
  dof_id_type node1;
  Point p0;
  Point p1;
};

// Helps maps treat an edge the same regardless of direction
static std::pair<dof_id_type, dof_id_type>
edgeKey(const dof_id_type node0, const dof_id_type node1)
{
  return {std::min(node0, node1), std::max(node0, node1)};
}

// Used to reject degenerate polygon faces
static bool
hasNonzeroArea3D(const std::vector<Point> & points, const Real tol = 1e-12)
{
  for (const auto i : index_range(points))
    for (const auto j : make_range(i + 1, points.size()))
      for (const auto k : make_range(j + 1, points.size()))
        if (((points[j] - points[i]).cross(points[k] - points[i])).norm() > tol)
          return true;

  return false;
}

// For adding polygon faces to polyhedrons
static void
addSidePoints3D(std::vector<std::vector<Point>> & sides,
                const std::vector<Point> & side_points,
                const Real tol = 1e-12)
{
  std::vector<Point> unique_side_points;

  for (const auto & point : side_points)
    addUniquePoint(unique_side_points, point, tol);

  if (unique_side_points.size() < 3 || !hasNonzeroArea3D(unique_side_points, tol))
    return;

  sides.push_back(unique_side_points);
}

// For finding a normal of the first non-collinear triple in a face's points
static Point
faceNormal3D(const std::vector<Point> & points, const Real tol = 1e-12)
{
  for (const auto i : index_range(points))
    for (const auto j : make_range(i + 1, points.size()))
      for (const auto k : make_range(j + 1, points.size()))
      {
        const Point normal = (points[j] - points[i]).cross(points[k] - points[i]);

        if (normal.norm() > tol)
          return normal;
      }

  return Point();
}

static void
addUniqueDirection3D(std::vector<Point> & directions,
                     const Point & direction,
                     const Real tol = 1e-8)
{
  if (direction.norm() <= tol)
    return;

  const Point unit_direction = direction / direction.norm();

  for (const auto & existing_direction : directions)
    if ((unit_direction - existing_direction).norm() <= tol)
      return;

  directions.push_back(unit_direction);
}

// Sorting, but 3D
static std::vector<Point>
sortPointsAroundAxis3D(const std::vector<Point> & unsorted_points,
                       const Point & axis,
                       const Real tol = 1e-12)
{
  std::vector<Point> points;

  for (const auto & point : unsorted_points)
    addUniquePoint(points, point, tol);

  if (points.size() < 3 || axis.norm() <= tol)
    return points;

  const Point axis_unit = axis / axis.norm();
  Point reference = std::abs(axis_unit(0)) < 0.9 ? Point(1.0, 0.0, 0.0) : Point(0.0, 1.0, 0.0);
  Point e1 = reference - (reference * axis_unit) * axis_unit;

  if (e1.norm() <= tol)
  {
    reference = Point(0.0, 0.0, 1.0);
    e1 = reference - (reference * axis_unit) * axis_unit;
  }

  if (e1.norm() <= tol)
    return points;

  e1 /= e1.norm();

  Point e2 = axis_unit.cross(e1);

  if (e2.norm() <= tol)
    return points;

  e2 /= e2.norm();

  Point center;

  for (const auto & point : points)
    center += point;

  center /= points.size();

  std::sort(points.begin(),
            points.end(),
            [&center, &e1, &e2](const Point & a, const Point & b)
            {
              const Point da = a - center;
              const Point db = b - center;

              return std::atan2(da * e2, da * e1) < std::atan2(db * e2, db * e1);
            });

  if (axis * faceNormal3D(points) < 0.0)
    std::reverse(points.begin(), points.end());

  return points;
}

static Real
polyhedronScale3D(const std::vector<std::vector<Point>> & side_points)
{
  Real scale = 1.0;

  for (const auto & side : side_points)
    for (const auto i : index_range(side))
      scale = std::max(scale, (side[(i + 1) % side.size()] - side[i]).norm());

  return scale;
}

static bool
pointInsideTriangulatedSurface3D(const Point & point,
                                 const std::vector<std::vector<Point>> & surface_triangles,
                                 const Real tol = 1e-12)
{
  Real solid_angle = 0.0;

  for (const auto & triangle : surface_triangles)
  {
    if (triangle.size() != 3)
      continue;

    const Point a = triangle[0] - point;
    const Point b = triangle[1] - point;
    const Point c = triangle[2] - point;
    const Real a_norm = a.norm();
    const Real b_norm = b.norm();
    const Real c_norm = c.norm();

    if (a_norm <= tol || b_norm <= tol || c_norm <= tol)
      return true;

    const Real numerator = a * b.cross(c);
    const Real denominator =
        a_norm * b_norm * c_norm + (a * b) * c_norm + (b * c) * a_norm + (c * a) * b_norm;

    if (std::abs(numerator) <= tol && std::abs(denominator) <= tol)
      continue;

    solid_angle += 2.0 * std::atan2(numerator, denominator);
  }

  return std::abs(solid_angle) > libMesh::pi;
}

static Real
tetVolume6(const Point & point0, const Point & point1, const Point & point2, const Point & point3)
{
  return (point1 - point0).cross(point2 - point0) * (point3 - point0);
}

static bool
surfaceTriangles3D(const std::vector<std::vector<Point>> & side_points,
                   std::vector<std::vector<Point>> & surface_triangles,
                   const Real tol = 1e-12)
{
  surface_triangles.clear();

  for (const auto & side : side_points)
  {
    if (side.size() < 3)
      continue;

    for (const auto i : make_range(std::size_t(1), side.size() - 1))
    {
      const std::vector<Point> triangle = {side[0], side[i], side[i + 1]};

      if (hasNonzeroArea3D(triangle, tol))
        surface_triangles.push_back(triangle);
    }
  }

  return !surface_triangles.empty();
}

std::unique_ptr<MeshBase>
DualMeshGenerator::generate()
{
  const auto input_mesh = std::move(_input);

  if (!input_mesh->is_prepared())
    input_mesh->find_neighbors();

  const bool use_voronoi = _dual_mesh_type == "voronoi";
  const unsigned int mesh_dimension = input_mesh->mesh_dimension();

  if (mesh_dimension != 2 && mesh_dimension != 3)
    mooseError("DualMeshGenerator currently only supports 2D and 3D Meshes");
  if (mesh_dimension == 3 && use_voronoi)
    mooseError("DualMeshGenerator does not support Voronoi duals for 3D meshes");
  if (mesh_dimension == 3)
  {
    auto dualMesh = buildReplicatedMesh(3);

    std::unordered_map<dof_id_type, std::vector<Point>> boundary_node_normals;
    std::map<std::pair<dof_id_type, dof_id_type>, std::vector<Point>> boundary_edge_normals;

    for (const auto & elem : input_mesh->element_ptr_range())
    {
      const Point elem_centroid = elem->true_centroid();

      for (const auto side : elem->side_index_range())
      {
        if (elem->neighbor_ptr(side) != nullptr)
          continue;

        auto side_elem = elem->build_side_ptr(side);
        std::vector<Point> side_points;
        side_points.reserve(side_elem->n_vertices());

        for (const auto n : make_range(side_elem->n_vertices()))
          side_points.push_back(side_elem->point(n));

        Point normal = faceNormal3D(side_points);

        if (normal.norm() <= 1e-12)
          continue;

        const Point face_centroid = side_elem->true_centroid();

        if (normal * (elem_centroid - face_centroid) > 0.0)
          normal = -1.0 * normal;

        for (const auto n : make_range(side_elem->n_vertices()))
          addUniqueDirection3D(boundary_node_normals[side_elem->node_id(n)], normal);

        for (const auto n : make_range(side_elem->n_vertices()))
        {
          const dof_id_type node0 = side_elem->node_id(n);
          const dof_id_type node1 = side_elem->node_id((n + 1) % side_elem->n_vertices());

          addUniqueDirection3D(boundary_edge_normals[edgeKey(node0, node1)], normal);
        }
      }
    }

    std::unordered_set<dof_id_type> boundary_vertex_nodes;

    for (const auto & node_normals : boundary_node_normals)
      if (node_normals.second.size() > 1)
        boundary_vertex_nodes.insert(node_normals.first);

    std::map<std::pair<dof_id_type, dof_id_type>, Point> boundary_edge_midpoints;

    for (const auto & edge_normals : boundary_edge_normals)
    {
      const auto & edge = edge_normals.first;

      if (edge_normals.second.size() > 1 && boundary_vertex_nodes.count(edge.first) &&
          boundary_vertex_nodes.count(edge.second))
        boundary_edge_midpoints[edge] =
            0.5 * (*input_mesh->node_ptr(edge.first) + *input_mesh->node_ptr(edge.second));
    }

    std::unordered_map<dof_id_type, std::vector<const Elem *>> source_node_to_elems;

    for (const auto & elem : input_mesh->element_ptr_range())
      for (const auto n : make_range(elem->n_vertices()))
        source_node_to_elems[elem->node_id(n)].push_back(elem);

    const auto tryAddPolyhedron = [&](MeshBase & mesh,
                                      const std::vector<std::vector<Point>> & side_points) -> bool
    {
      if (side_points.size() < 4)
        return false;

      std::vector<Node *> local_nodes;
      std::vector<std::shared_ptr<libMesh::Polygon>> sides;

      const auto getLocalNode = [&](const Point & point)
      {
        for (auto * const node : local_nodes)
          if (MooseUtils::absoluteFuzzyEqual(*node, point))
            return node;

        Node * const node = mesh.add_point(point);
        local_nodes.push_back(node);

        return node;
      };

      sides.reserve(side_points.size());

      for (const auto & side : side_points)
      {
        auto polygon = std::make_shared<libMesh::C0Polygon>(side.size());

        for (const auto i : index_range(side))
          polygon->set_node(i, getLocalNode(side[i]));

        sides.push_back(polygon);
      }

      libmesh_try
      {
        std::unique_ptr<libMesh::Node> mid_elem_node;
        auto dual_elem = std::make_unique<libMesh::C0Polyhedron>(sides, mid_elem_node);

        if (mid_elem_node)
          mesh.add_node(std::move(mid_elem_node));

        mesh.add_elem(std::move(dual_elem));
      }
      libmesh_catch(const libMesh::NotImplemented &) { return false; }
      libmesh_catch(const libMesh::LogicError &) { return false; }

      return true;
    };

    const auto addPolyhedron = [&](const std::vector<std::vector<Point>> & side_points)
    {
      auto trial_mesh = buildReplicatedMesh(3);

      if (!tryAddPolyhedron(*trial_mesh, side_points))
        return false;

      return tryAddPolyhedron(*dualMesh, side_points);
    };

    const auto addTetrahedralizedPolyhedron =
        [&](const std::vector<std::vector<Point>> & side_points) -> bool
    {
      if (side_points.size() < 4)
        return false;

      const Real tol = std::max(_geometry_relative_tol, Real(1e-12));
      const Real length_tol = tol * polyhedronScale3D(side_points);
      const Real volume_tol = length_tol * length_tol * length_tol;
      std::vector<std::vector<Point>> surface_triangles;

      if (!surfaceTriangles3D(side_points, surface_triangles, length_tol))
        return false;

      std::vector<Point> unique_points;

      for (const auto & side : side_points)
        for (const auto & point : side)
          addUniquePoint(unique_points, point, length_tol);

      if (unique_points.size() < 4)
        return false;

      Point vertex_center;

      for (const auto & point : unique_points)
        vertex_center += point;

      vertex_center /= unique_points.size();

      Point surface_center;
      Real surface_weight = 0.0;

      for (const auto & triangle : surface_triangles)
      {
        const Real triangle_weight =
            (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]).norm();

        surface_center += triangle_weight * (triangle[0] + triangle[1] + triangle[2]) / 3.0;
        surface_weight += triangle_weight;
      }

      if (surface_weight > length_tol * length_tol)
        surface_center /= surface_weight;
      else
        surface_center = vertex_center;

      std::vector<Point> candidate_points;

      const auto addCandidatePoint = [&](const Point & point)
      {
        if (!pointInsideTriangulatedSurface3D(point, surface_triangles, length_tol))
          return;

        addUniquePoint(candidate_points, point, length_tol);
      };

      addCandidatePoint(vertex_center);
      addCandidatePoint(surface_center);
      addCandidatePoint(0.5 * (vertex_center + surface_center));

      for (const auto & point : unique_points)
        for (const auto fraction : {0.25, 0.5, 0.75})
        {
          addCandidatePoint(vertex_center + fraction * (point - vertex_center));
          addCandidatePoint(surface_center + fraction * (point - surface_center));
        }

      for (const auto & triangle : surface_triangles)
      {
        const Point triangle_center = (triangle[0] + triangle[1] + triangle[2]) / 3.0;
        Point triangle_normal = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);

        for (const auto fraction : {0.25, 0.5, 0.75})
        {
          addCandidatePoint(vertex_center + fraction * (triangle_center - vertex_center));
          addCandidatePoint(surface_center + fraction * (triangle_center - surface_center));
        }

        if (triangle_normal.norm() > length_tol)
        {
          triangle_normal /= triangle_normal.norm();

          Real triangle_scale = 1.0;

          for (const auto i : index_range(triangle))
            triangle_scale = std::max(triangle_scale,
                                      (triangle[(i + 1) % triangle.size()] - triangle[i]).norm());

          for (const auto sign : {-1.0, 1.0})
            for (const auto fraction : {0.05, 0.1, 0.2, 0.35})
              addCandidatePoint(triangle_center +
                                sign * fraction * triangle_scale * triangle_normal);
        }
      }

      if (candidate_points.empty())
        return false;

      const auto candidateCoversTriangle =
          [&](const Point & candidate_point, const std::vector<Point> & triangle)
      {
        if (std::abs(tetVolume6(candidate_point, triangle[0], triangle[1], triangle[2])) <=
            volume_tol)
          return false;

        const Point triangle_center = (triangle[0] + triangle[1] + triangle[2]) / 3.0;

        for (const auto fraction : {0.25, 0.5})
          if (!pointInsideTriangulatedSurface3D(candidate_point +
                                                    fraction * (triangle_center - candidate_point),
                                                surface_triangles,
                                                length_tol))
            return false;

        const Point tet_center = (candidate_point + triangle[0] + triangle[1] + triangle[2]) / 4.0;

        return pointInsideTriangulatedSurface3D(tet_center, surface_triangles, length_tol);
      };

      const auto candidateWeaklyCoversTriangle =
          [&](const Point & candidate_point, const std::vector<Point> & triangle)
      {
        if (std::abs(tetVolume6(candidate_point, triangle[0], triangle[1], triangle[2])) <=
            volume_tol)
          return false;

        const Point tet_center = (candidate_point + triangle[0] + triangle[1] + triangle[2]) / 4.0;

        return pointInsideTriangulatedSurface3D(tet_center, surface_triangles, length_tol);
      };

      std::vector<bool> covered_triangles(surface_triangles.size(), false);
      std::vector<std::pair<Point, std::vector<std::size_t>>> selected_candidate_triangles;

      while (std::find(covered_triangles.begin(), covered_triangles.end(), false) !=
             covered_triangles.end())
      {
        std::vector<std::size_t> best_triangle_indices;
        Point best_candidate_point;

        for (const auto & candidate_point : candidate_points)
        {
          std::vector<std::size_t> candidate_triangle_indices;

          for (const auto triangle_i : index_range(surface_triangles))
            if (!covered_triangles[triangle_i] &&
                candidateCoversTriangle(candidate_point, surface_triangles[triangle_i]))
              candidate_triangle_indices.push_back(triangle_i);

          if (candidate_triangle_indices.size() > best_triangle_indices.size())
          {
            best_candidate_point = candidate_point;
            best_triangle_indices = candidate_triangle_indices;
          }
        }

        if (best_triangle_indices.empty())
          for (const auto & candidate_point : candidate_points)
          {
            std::vector<std::size_t> candidate_triangle_indices;

            for (const auto triangle_i : index_range(surface_triangles))
              if (!covered_triangles[triangle_i] &&
                  candidateWeaklyCoversTriangle(candidate_point, surface_triangles[triangle_i]))
                candidate_triangle_indices.push_back(triangle_i);

            if (candidate_triangle_indices.size() > best_triangle_indices.size())
            {
              best_candidate_point = candidate_point;
              best_triangle_indices = candidate_triangle_indices;
            }
          }

        if (best_triangle_indices.empty())
          return false;

        for (const auto triangle_i : best_triangle_indices)
          covered_triangles[triangle_i] = true;

        selected_candidate_triangles.push_back({best_candidate_point, best_triangle_indices});
      }

      std::vector<Node *> local_nodes;

      const auto getLocalNode = [&](const Point & point)
      {
        for (auto * const node : local_nodes)
          if (MooseUtils::absoluteFuzzyEqual(*node, point, length_tol))
            return node;

        Node * const node = dualMesh->add_point(point);
        local_nodes.push_back(node);

        return node;
      };

      for (const auto & candidate_triangles : selected_candidate_triangles)
      {
        Node * const interior_node = getLocalNode(candidate_triangles.first);

        for (const auto triangle_i : candidate_triangles.second)
        {
          const auto & triangle = surface_triangles[triangle_i];

          if (std::abs(tetVolume6(
                  candidate_triangles.first, triangle[0], triangle[1], triangle[2])) <= volume_tol)
            continue;

          auto tet = std::make_unique<Tet4>();
          tet->set_node(0) = interior_node;
          tet->set_node(1) = getLocalNode(triangle[0]);

          if (tetVolume6(candidate_triangles.first, triangle[0], triangle[1], triangle[2]) > 0.0)
          {
            tet->set_node(2) = getLocalNode(triangle[1]);
            tet->set_node(3) = getLocalNode(triangle[2]);
          }
          else
          {
            tet->set_node(2) = getLocalNode(triangle[2]);
            tet->set_node(3) = getLocalNode(triangle[1]);
          }

          dualMesh->add_elem(std::move(tet));
        }
      }

      return true;
    };

    const auto addPolyhedronOrTetrahedralize =
        [&](const std::vector<std::vector<Point>> & side_points)
    {
      if (addPolyhedron(side_points))
        return;

      if (!addTetrahedralizedPolyhedron(side_points))
        mooseError("Could not tetrahedralize rejected non-convex 3D dual polyhedron.");
    };

    // Build one dual cell around each primal node. Concave cells are broken into TET4s using
    // checked interior points.
    for (const auto & node_elems : source_node_to_elems)
    {
      const dof_id_type source_node_id = node_elems.first;
      const Point & source_point = *input_mesh->node_ptr(source_node_id);
      std::map<std::pair<dof_id_type, dof_id_type>, std::vector<Point>> edge_to_points;
      std::map<std::pair<dof_id_type, dof_id_type>, std::vector<Point>>
          midpoint_boundary_face_centroids;
      std::vector<Point> boundary_face_centroids;
      Point boundary_normal;

      for (const auto & elem : node_elems.second)
      {
        const Point elem_centroid = elem->true_centroid();

        for (const auto side : elem->side_index_range())
        {
          auto side_elem = elem->build_side_ptr(side);
          const Point face_centroid = side_elem->true_centroid();

          std::vector<dof_id_type> current_side_node_ids;
          std::vector<Point> current_side_points;
          current_side_node_ids.reserve(side_elem->n_vertices());
          current_side_points.reserve(side_elem->n_vertices());

          // For interior dual elements
          for (const auto n : make_range(side_elem->n_vertices()))
          {
            current_side_node_ids.push_back(side_elem->node_id(n));
            current_side_points.push_back(side_elem->point(n));
          }

          const auto source_node_it =
              std::find(current_side_node_ids.begin(), current_side_node_ids.end(), source_node_id);

          if (source_node_it == current_side_node_ids.end())
            continue;

          const auto source_side_index =
              cast_int<unsigned int>(source_node_it - current_side_node_ids.begin());
          const dof_id_type previous_node_id =
              current_side_node_ids[(source_side_index + current_side_node_ids.size() - 1) %
                                    current_side_node_ids.size()];
          const dof_id_type next_node_id =
              current_side_node_ids[(source_side_index + 1) % current_side_node_ids.size()];

          auto & previous_edge_points = edge_to_points[edgeKey(source_node_id, previous_node_id)];
          auto & next_edge_points = edge_to_points[edgeKey(source_node_id, next_node_id)];

          addUniquePoint(previous_edge_points, elem_centroid);
          addUniquePoint(next_edge_points, elem_centroid);

          if (elem->neighbor_ptr(side) == nullptr)
          {
            addUniquePoint(previous_edge_points, face_centroid);
            addUniquePoint(next_edge_points, face_centroid);
            addUniquePoint(boundary_face_centroids, face_centroid);

            const auto previous_midpoint_it =
                boundary_edge_midpoints.find(edgeKey(source_node_id, previous_node_id));
            const auto next_midpoint_it =
                boundary_edge_midpoints.find(edgeKey(source_node_id, next_node_id));

            if (boundary_vertex_nodes.count(source_node_id))
            {
              if (previous_midpoint_it != boundary_edge_midpoints.end())
              {
                addUniquePoint(previous_edge_points, previous_midpoint_it->second);
                addUniquePoint(midpoint_boundary_face_centroids[previous_midpoint_it->first],
                               face_centroid);
              }

              if (next_midpoint_it != boundary_edge_midpoints.end())
              {
                addUniquePoint(next_edge_points, next_midpoint_it->second);
                addUniquePoint(midpoint_boundary_face_centroids[next_midpoint_it->first],
                               face_centroid);
              }
            }

            Point normal = faceNormal3D(current_side_points);

            if (normal.norm() > 1e-12)
            {
              if (normal * (elem_centroid - face_centroid) > 0.0)
                normal = -1.0 * normal;

              boundary_normal += normal / normal.norm();
            }
          }
        }
      }

      std::vector<std::vector<Point>> polyhedron_side_points;
      std::vector<std::pair<Point, Point>> midpoint_split_boundary_faces;

      for (const auto & edge_points : edge_to_points)
      {
        if (edge_points.second.size() < 3)
          continue;

        const dof_id_type other_node_id = edge_points.first.first == source_node_id
                                              ? edge_points.first.second
                                              : edge_points.first.first;
        const Point edge_axis = *input_mesh->node_ptr(other_node_id) - source_point;
        const auto sorted_edge_points = sortPointsAroundAxis3D(edge_points.second, edge_axis);

        addSidePoints3D(polyhedron_side_points, sorted_edge_points);
      }

      for (const auto & midpoint_face_centroids : midpoint_boundary_face_centroids)
      {
        const auto midpoint_it = boundary_edge_midpoints.find(midpoint_face_centroids.first);

        if (midpoint_it == boundary_edge_midpoints.end())
          continue;

        for (const auto & face_centroid : midpoint_face_centroids.second)
          addSidePoints3D(polyhedron_side_points,
                          {source_point, midpoint_it->second, face_centroid});

        for (const auto i : index_range(midpoint_face_centroids.second))
          for (const auto j : make_range(i + 1, midpoint_face_centroids.second.size()))
            midpoint_split_boundary_faces.push_back(
                {midpoint_face_centroids.second[i], midpoint_face_centroids.second[j]});
      }

      if (boundary_face_centroids.size() >= 2)
      {
        const Point boundary_axis = boundary_normal.norm() > 1e-12 ? boundary_normal : source_point;
        const auto sorted_boundary_points =
            sortPointsAroundAxis3D(boundary_face_centroids, boundary_axis);

        const auto boundaryFaceWasSplit = [&](const Point & point0, const Point & point1)
        {
          for (const auto & split_boundary_face : midpoint_split_boundary_faces)
            if ((MooseUtils::absoluteFuzzyEqual(point0, split_boundary_face.first) &&
                 MooseUtils::absoluteFuzzyEqual(point1, split_boundary_face.second)) ||
                (MooseUtils::absoluteFuzzyEqual(point0, split_boundary_face.second) &&
                 MooseUtils::absoluteFuzzyEqual(point1, split_boundary_face.first)))
              return true;

          return false;
        };

        if (boundary_vertex_nodes.count(source_node_id))
        {
          if (sorted_boundary_points.size() == 2)
          {
            if (!boundaryFaceWasSplit(sorted_boundary_points[0], sorted_boundary_points[1]))
              addSidePoints3D(polyhedron_side_points,
                              {source_point, sorted_boundary_points[0], sorted_boundary_points[1]});
          }
          else
            for (const auto i : index_range(sorted_boundary_points))
            {
              const Point & point0 = sorted_boundary_points[i];
              const Point & point1 =
                  sorted_boundary_points[(i + 1) % sorted_boundary_points.size()];

              if (!boundaryFaceWasSplit(point0, point1))
                addSidePoints3D(polyhedron_side_points, {source_point, point0, point1});
            }
        }
        else if (sorted_boundary_points.size() >= 3)
          addSidePoints3D(polyhedron_side_points, sorted_boundary_points);
      }

      if (polyhedron_side_points.size() < 4)
        continue;

      addPolyhedronOrTetrahedralize(polyhedron_side_points);
    }

    dualMesh->unset_is_prepared();
    return dynamic_pointer_cast<MeshBase>(dualMesh);
  }

  // Rest of file is 2D mesh treatment
  // Calculating tolerance, so they are relative to the input mesh bounding box
  const auto input_bounding_box = MeshTools::create_bounding_box(*input_mesh);
  const Point mesh_extent = input_bounding_box.max() - input_bounding_box.min();
  const Real mesh_scale = std::max(std::max(std::abs(mesh_extent(0)), std::abs(mesh_extent(1))),
                                   std::numeric_limits<Real>::min());
  const Real length_tol = _geometry_relative_tol * mesh_scale;
  const Real area_tol = _geometry_relative_tol * mesh_scale * mesh_scale;
  const Real parameter_tol = _geometry_relative_tol;

  std::unordered_map<dof_id_type, Point> boundary_node_points;
  std::unordered_map<dof_id_type, std::vector<Point>> boundary_node_midpoints;
  std::vector<BoundarySegment> physical_boundary_segments;

  // Looping through primal elements
  for (const auto & elem : input_mesh->element_ptr_range())
  {
    for (const auto side : elem->side_index_range())
    {
      if (elem->neighbor_ptr(side) == nullptr)
      {
        // If element has nullptr neighbor we want to record vertices, midpoints, and other info
        auto side_elem = elem->build_side_ptr(side);

        if (side_elem->n_nodes() == 2)
        {
          const dof_id_type node0 = side_elem->node_id(0);
          const dof_id_type node1 = side_elem->node_id(1);

          physical_boundary_segments.push_back(
              {node0, node1, side_elem->point(0), side_elem->point(1)});

          boundary_node_points[node0] = side_elem->point(0);
          boundary_node_points[node1] = side_elem->point(1);

          const Point midpoint = 0.5 * (side_elem->point(0) + side_elem->point(1));
          boundary_node_midpoints[node0].push_back(midpoint);
          boundary_node_midpoints[node1].push_back(midpoint);
        }
      }
    }
  }

  std::unordered_map<dof_id_type, std::vector<std::size_t>> boundary_node_to_segments;
  // Building map for boundary segments
  for (const auto i : index_range(physical_boundary_segments))
  {
    boundary_node_to_segments[physical_boundary_segments[i].node0].push_back(i);
    boundary_node_to_segments[physical_boundary_segments[i].node1].push_back(i);
  }
  // Helper so given a node we can get another node than shares a line segment with it
  const auto otherBoundaryNode = [&](const BoundarySegment & segment,
                                     const dof_id_type node_id) -> dof_id_type
  { return segment.node0 == node_id ? segment.node1 : segment.node0; };
  // Set of vertices we want to add or preserve
  std::unordered_set<dof_id_type> boundary_vertex_nodes;

  for (const auto & node_segments : boundary_node_to_segments)
  {
    const dof_id_type node_id = node_segments.first;
    const auto & segment_ids = node_segments.second;

    if (segment_ids.size() != 2)
    {
      boundary_vertex_nodes.insert(node_id);
      continue;
    }

    const Point & p = boundary_node_points[node_id];
    const Point v0 = boundary_node_points[otherBoundaryNode(
                         physical_boundary_segments[segment_ids[0]], node_id)] -
                     p;
    const Point v1 = boundary_node_points[otherBoundaryNode(
                         physical_boundary_segments[segment_ids[1]], node_id)] -
                     p;

    const Real norm_product = v0.norm() * v1.norm();

    if (norm_product < length_tol * length_tol)
    {
      boundary_vertex_nodes.insert(node_id);
      continue;
    }

    Real cos_angle = (v0 * v1) / norm_product;
    cos_angle = std::max(Real(-1.0), std::min(Real(1.0), cos_angle));

    if (std::abs(libMesh::pi - std::acos(cos_angle)) > _boundary_node_angular_tol)
      boundary_vertex_nodes.insert(node_id); // Records vertices based on tolerance criteria
  }

  std::vector<std::pair<Point, Point>> boundary_clip_segments;
  std::vector<bool> used_boundary_segments(physical_boundary_segments.size(), false);
  // Helper to find the boundary line we want to clip the mesh along
  const auto traceBoundaryClipSegment =
      [&](const dof_id_type start_node, const std::size_t start_segment_id)
  {
    dof_id_type current_node = start_node;
    std::size_t current_segment_id = start_segment_id;

    while (current_segment_id < physical_boundary_segments.size() &&
           !used_boundary_segments[current_segment_id])
    {
      used_boundary_segments[current_segment_id] = true;

      const auto & segment = physical_boundary_segments[current_segment_id];
      const dof_id_type next_node = otherBoundaryNode(segment, current_node);

      if (next_node == start_node || boundary_vertex_nodes.count(next_node))
      {
        boundary_clip_segments.push_back(
            {boundary_node_points[start_node], boundary_node_points[next_node]});
        return;
      }

      const auto node_segment_it = boundary_node_to_segments.find(next_node);

      if (node_segment_it == boundary_node_to_segments.end())
        return;

      std::size_t next_segment_id = physical_boundary_segments.size();

      for (const auto candidate_segment_id : node_segment_it->second)
        if (!used_boundary_segments[candidate_segment_id])
        {
          next_segment_id = candidate_segment_id;
          break;
        }

      current_node = next_node;
      current_segment_id = next_segment_id;
    }
  };
  // Now we actually trace these lines
  for (const auto boundary_vertex_node : boundary_vertex_nodes)
  {
    const auto node_segment_it = boundary_node_to_segments.find(boundary_vertex_node);

    if (node_segment_it == boundary_node_to_segments.end())
      continue;

    for (const auto segment_id : node_segment_it->second)
      if (!used_boundary_segments[segment_id])
        traceBoundaryClipSegment(boundary_vertex_node, segment_id);
  }
  // Now we clip along the traced lines
  for (const auto segment_id : index_range(physical_boundary_segments))
    if (!used_boundary_segments[segment_id])
      traceBoundaryClipSegment(physical_boundary_segments[segment_id].node0, segment_id);

  // Now for actually finding the dual
  std::vector<Point> dual_centers;
  std::unordered_map<dof_id_type, dof_id_type> source_elem_to_center_id;
  std::unordered_map<dof_id_type, std::vector<const Elem *>> source_node_to_elems;
  std::unique_ptr<ReplicatedMesh> tri_mesh;

  if (use_voronoi)
  {
    // tri_mesh input is the pure primal mesh
    tri_mesh = buildReplicatedMesh(2);

    for (const auto & node : input_mesh->node_ptr_range())
    {
      Node * new_node = tri_mesh->add_point(*node);

      auto node_elem = std::make_unique<NodeElem>();
      node_elem->set_node(0) = new_node;
      tri_mesh->add_elem(std::move(node_elem));
    }

    // We're circumscribing this inside a large box that is scaled to the bounding box of the primal
    // mesh, then we'll take the Delaunay triangulation
    const Real outer_padding = 10.0 * mesh_scale;
    const Point outer_min(input_bounding_box.min()(0) - outer_padding,
                          input_bounding_box.min()(1) - outer_padding,
                          0.0);
    const Point outer_max(input_bounding_box.max()(0) + outer_padding,
                          input_bounding_box.max()(1) + outer_padding,
                          0.0);

    Node * p0 = tri_mesh->add_point(Point(outer_min(0), outer_min(1), 0.0));
    Node * p1 = tri_mesh->add_point(Point(outer_max(0), outer_min(1), 0.0));
    Node * p2 = tri_mesh->add_point(Point(outer_max(0), outer_max(1), 0.0));
    Node * p3 = tri_mesh->add_point(Point(outer_min(0), outer_max(1), 0.0));

    auto big_square = std::make_unique<Quad4>();

    big_square->set_node(0) = p0;
    big_square->set_node(1) = p1;
    big_square->set_node(2) = p2;
    big_square->set_node(3) = p3;

    tri_mesh->add_elem(std::move(big_square));

    Poly2TriTriangulator triangulator(dynamic_cast<UnstructuredMesh &>(*tri_mesh));
    triangulator.triangulation_type() = libMesh::TriangulatorInterface::PSLG;
    triangulator.minimum_angle() = 0;
    triangulator.desired_area() = 0;

    triangulator.insert_extra_points() = false;     // Don't add new nodes!
    triangulator.smooth_after_generating() = false; // Don't mess up geometry

    triangulator.triangulate();

    // Now for Voronoi, we loop over the triangles and get the circumcenters, etc
    for (const auto & tri_elem : tri_mesh->element_ptr_range())
    {
      if (tri_elem->n_vertices() != 3)
        continue;

      const dof_id_type center_id = dual_centers.size();

      dual_centers.push_back(
          libMesh::circumcenter(tri_elem->point(0), tri_elem->point(1), tri_elem->point(2)));
      source_elem_to_center_id[tri_elem->id()] = center_id;

      for (const auto n : make_range(tri_elem->n_nodes()))
        source_node_to_elems[tri_elem->node_id(n)].push_back(tri_elem);
    }
  }
  else // For barycentric, loop over elements and get centroids
  {
    for (const auto & elem : input_mesh->element_ptr_range())
    {
      const dof_id_type center_id = dual_centers.size();

      dual_centers.push_back(elem->true_centroid());
      source_elem_to_center_id[elem->id()] = center_id;

      for (const auto n : make_range(elem->n_nodes()))
        source_node_to_elems[elem->node_id(n)].push_back(elem);
    }
  }
  // Now we've gathered our centers and boundary info, we can now build the dual
  auto dualMesh = buildReplicatedMesh(2);
  // Helper to determine if points are within the primal boundary
  const auto pointInsideBoundary = [&](const Point & point)
  {
    bool inside = false;

    for (const auto & segment : physical_boundary_segments)
    {
      if (pointOnSegment2D(point, segment.p0, segment.p1, length_tol, area_tol))
        return true;

      if ((segment.p0(1) > point(1)) != (segment.p1(1) > point(1)))
      {
        const Real intersection_x = (segment.p1(0) - segment.p0(0)) * (point(1) - segment.p0(1)) /
                                        (segment.p1(1) - segment.p0(1)) +
                                    segment.p0(0);

        if (point(0) < intersection_x)
          inside = !inside;
      }
    }

    return inside;
  };
  // Actual clipping algorithm
  const auto clipDualPolygonToBoundary = [&](const std::vector<Point> & dual_points)
  {
    std::vector<Point> clipped_points;

    if (dual_points.size() < 3)
      return clipped_points;

    for (const auto & point : dual_points)
      if (pointInsideBoundary(point))
        addUniquePoint(clipped_points, point, length_tol);

    for (const auto i : index_range(dual_points))
    {
      const Point & p0 = dual_points[i];
      const Point & p1 = dual_points[(i + 1) % dual_points.size()];

      for (const auto & boundary_segment : boundary_clip_segments)
        addSegmentIntersections2D(clipped_points,
                                  p0,
                                  p1,
                                  boundary_segment.first,
                                  boundary_segment.second,
                                  length_tol,
                                  area_tol,
                                  parameter_tol);
    }

    for (const auto & boundary_segment : boundary_clip_segments)
    {
      if (pointInPolygon2D(boundary_segment.first, dual_points, length_tol, area_tol))
        addUniquePoint(clipped_points, boundary_segment.first, length_tol);

      if (pointInPolygon2D(boundary_segment.second, dual_points, length_tol, area_tol))
        addUniquePoint(clipped_points, boundary_segment.second, length_tol);
    }

    if (clipped_points.size() < 3)
      return clipped_points;

    Point center;

    for (const auto & point : clipped_points)
      center += point;

    center /= clipped_points.size();

    std::sort(clipped_points.begin(),
              clipped_points.end(),
              [&center](const Point & a, const Point & b)
              {
                return std::atan2(a(1) - center(1), a(0) - center(0)) <
                       std::atan2(b(1) - center(1), b(0) - center(0));
              });

    std::vector<Point> unique_clipped_points;

    for (const auto & point : clipped_points)
      addUniquePoint(unique_clipped_points, point, length_tol);

    if (unique_clipped_points.size() > 1 &&
        MooseUtils::absoluteFuzzyEqual(
            unique_clipped_points.front(), unique_clipped_points.back(), length_tol))
      unique_clipped_points.pop_back();

    return unique_clipped_points;
  };

  // Helper for finding if a node is a boundary vertex
  const auto isBoundaryVertexPoint = [&](const Point & point)
  {
    for (const auto boundary_vertex_node : boundary_vertex_nodes)
    {
      const auto point_it = boundary_node_points.find(boundary_vertex_node);

      if (point_it != boundary_node_points.end() &&
          MooseUtils::absoluteFuzzyEqual(point, point_it->second, length_tol))
        return true;
    }

    return false;
  };

  // Helper for finding if a node lies on a boundary segment, for our clipping routine
  const auto isBoundarySegmentPoint = [&](const Point & point)
  {
    for (const auto & boundary_segment : boundary_clip_segments)
      if (pointOnSegment2D(
              point, boundary_segment.first, boundary_segment.second, length_tol, area_tol))
        return true;

    return false;
  };
  // Concave element indexing helper, used in our later decomposition of concave elems.
  const auto concaveBoundaryVertexIndex = [&](const std::vector<Point> & points) -> std::size_t
  {
    if (points.size() < 4)
      return points.size();

    const Real signed_area = polygonSignedArea2D(points);

    if (std::abs(signed_area) < area_tol)
      return points.size();

    const Real orientation = signed_area > 0.0 ? 1.0 : -1.0;

    for (const auto i : index_range(points))
    {
      const Point & previous = points[(i + points.size() - 1) % points.size()];
      const Point & current = points[i];
      const Point & next = points[(i + 1) % points.size()];

      if (isBoundaryVertexPoint(current) &&
          orientation * cross2D(previous, current, next) < -area_tol)
        return i;
    }

    return points.size();
  };

  const auto addDualElement = [&](const std::vector<Point> & points)
  {
    auto dual_elem = std::make_unique<libMesh::C0Polygon>(points.size());

    for (const auto i : index_range(points))
      dual_elem->set_node(i, dualMesh->add_point(points[i]));

    // Fixing flipped elems
    if (dual_elem->is_flipped())
    {
      auto reversed_elem = std::make_unique<libMesh::C0Polygon>(points.size());

      for (const auto i : index_range(points))
        reversed_elem->set_node(i, dualMesh->add_point(points[points.size() - 1 - i]));

      dual_elem = std::move(reversed_elem);
    }

    if (!dual_elem->is_flipped())
      dualMesh->add_elem(std::move(dual_elem));
  };

  // Build one dual element around each source node.
  for (const auto & node_elems : source_node_to_elems)
  {
    const dof_id_type source_node_id = node_elems.first;
    const auto & incident_elems = node_elems.second;

    std::vector<std::vector<unsigned int>> adjacency(incident_elems.size());

    for (const auto i : index_range(incident_elems))
      for (const auto j : make_range(i + 1, incident_elems.size()))
        if (elementsShareTwoNodes(incident_elems[i], incident_elems[j]))
        {
          adjacency[i].push_back(j);
          adjacency[j].push_back(i);
        }

    // Start at an endpoint for boundary chains, otherwise start anywhere.
    unsigned int start = 0;

    for (const auto i : index_range(adjacency))
      if (adjacency[i].size() == 1)
      {
        start = i;
        break;
      }

    std::vector<bool> used(incident_elems.size(), false);
    std::vector<dof_id_type> ordered_center_ids;

    unsigned int current = start;
    unsigned int previous = libMesh::invalid_uint;

    // Walk element adjacency to collect dual centers in connected order.
    while (true)
    {
      const Elem * elem = incident_elems[current];

      auto it = source_elem_to_center_id.find(elem->id());
      if (it != source_elem_to_center_id.end())
      {
        const dof_id_type center_id = it->second;

        if (std::find(ordered_center_ids.begin(), ordered_center_ids.end(), center_id) ==
            ordered_center_ids.end())
          ordered_center_ids.push_back(center_id);
      }

      used[current] = true;

      unsigned int next = libMesh::invalid_uint;

      for (const auto candidate : adjacency[current])
        if (candidate != previous && !used[candidate])
        {
          next = candidate;
          break;
        }

      if (next == libMesh::invalid_uint)
        break;

      previous = current;
      current = next;
    }

    std::vector<Point> dual_points;
    std::vector<Point> source_center_points;

    for (const auto center_id : ordered_center_ids)
    {
      addUniquePoint(dual_points, dual_centers[center_id], length_tol);
      source_center_points.push_back(dual_centers[center_id]);
    }

    if (!use_voronoi)
    {
      const auto boundary_midpoint_it = boundary_node_midpoints.find(source_node_id);

      if (boundary_midpoint_it != boundary_node_midpoints.end())
      {
        const auto boundary_point_it = boundary_node_points.find(source_node_id);

        if (boundary_vertex_nodes.count(source_node_id) &&
            boundary_point_it != boundary_node_points.end())
          addUniquePoint(dual_points, boundary_point_it->second, length_tol); // Add primal vertices

        for (const auto & midpoint : boundary_midpoint_it->second)
          addUniquePoint(dual_points, midpoint, length_tol);

        if (boundary_point_it != boundary_node_points.end() && !source_center_points.empty())
        {
          Point center_average; // Our source center is not quite the centroid, it's the midpoint
                                // between the dual nodes' centroid and the primal boundary point.

          for (const auto & center_point : source_center_points)
            center_average += center_point;

          center_average /= source_center_points.size();

          const Point sort_center = 0.5 * (boundary_point_it->second + center_average);

          // Sorting is necessary, as the centroid ordering is jumbled after adding vertices and
          // midpoints.
          std::sort(dual_points.begin(),
                    dual_points.end(),
                    [&sort_center](const Point & a, const Point & b)
                    {
                      return std::atan2(a(1) - sort_center(1), a(0) - sort_center(0)) <
                             std::atan2(b(1) - sort_center(1), b(0) - sort_center(0));
                    });
        }
      }
    }

    if (dual_points.size() >= 3)
      // Finally, clipping back to boundary.
      dual_points = clipDualPolygonToBoundary(dual_points);

    if (dual_points.size() < 3)
      continue;

    const std::size_t concave_vertex_index = concaveBoundaryVertexIndex(dual_points);

    if (concave_vertex_index < dual_points.size())
    {
      const Point corner_point = dual_points[concave_vertex_index];
      std::vector<std::pair<Point, Real>> sorted_points;

      // We use phi sorting only for concave dual elements, since the id's are always too jumbled to
      // sort by adjacency
      for (const auto i : index_range(dual_points))
      {
        if (i == concave_vertex_index)
          continue;

        const Real phi =
            std::atan2(dual_points[i](1) - corner_point(1), dual_points[i](0) - corner_point(0));
        sorted_points.push_back({dual_points[i], phi});
      }

      std::sort(sorted_points.begin(),
                sorted_points.end(),
                [](const auto & a, const auto & b) { return a.second < b.second; });

      std::vector<Point> fan_points = {corner_point};

      for (const auto i : index_range(sorted_points))
      {
        if (!isBoundarySegmentPoint(sorted_points[i].first))
          continue;

        const std::size_t next_i = (i + 1) % sorted_points.size();
        const std::size_t prev_i = (i + sorted_points.size() - 1) % sorted_points.size();

        // We pick the direction of fan-triangulating concave polygons such that we never create a
        // triangle that bridges across a boundary
        if (!isBoundarySegmentPoint(sorted_points[next_i].first))
        {
          // Now we can fan out triangles to break up otherwise concave dual elements
          for (const auto k : index_range(sorted_points))
            fan_points.push_back(sorted_points[(i + k) % sorted_points.size()].first);

          break;
        }

        if (!isBoundarySegmentPoint(sorted_points[prev_i].first))
        {
          for (const auto k : index_range(sorted_points))
            fan_points.push_back(
                sorted_points[(i + sorted_points.size() - k) % sorted_points.size()].first);

          break;
        }
      }

      if (fan_points.size() >= 3)
      {
        for (const auto i : make_range(std::size_t(1), fan_points.size() - 1))
        {
          const std::vector<Point> triangle_points = {
              fan_points[0], fan_points[i], fan_points[i + 1]};

          if (std::abs(cross2D(triangle_points[0], triangle_points[1], triangle_points[2])) >
              area_tol)
            addDualElement(triangle_points);
        }
      }

      continue;
    }

    addDualElement(dual_points);
  }

  dualMesh->unset_is_prepared();

  return dynamic_pointer_cast<MeshBase>(dualMesh);
}
