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
#include "MooseMeshUtils.h"
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
#include <array>
#include <cmath>
#include <limits>
#include <map>

registerMooseObject("MooseApp", DualMeshGenerator);

InputParameters
DualMeshGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  MooseEnum dual_mesh_type("voronoi barycentric", "barycentric");
  params.addParam<MooseEnum>(
      "dual_mesh_type",
      dual_mesh_type,
      "Whether to place dual nodes at Delaunay circumcenters or primal element centroids.");
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
    _boundary_node_angular_tol(getParam<Real>("boundary_node_angular_tol")),
    _geometry_relative_tol(getParam<Real>("geometry_relative_tol")),
    _dual_mesh_type(getParam<MooseEnum>("dual_mesh_type"))
{
}

// True only when two elements share a full edge.
static bool
elementsShareTwoNodes(const Elem * a, const Elem * b)
{
  unsigned int shared_nodes = 0;

  for (unsigned int i = 0; i < a->n_nodes(); ++i)
    for (unsigned int j = 0; j < b->n_nodes(); ++j)
      if (a->node_id(i) == b->node_id(j))
        ++shared_nodes;

  return shared_nodes == 2;
}

static Real
cross2D(const Point & a, const Point & b, const Point & c)
{
  return (b(0) - a(0)) * (c(1) - a(1)) - (b(1) - a(1)) * (c(0) - a(0));
}

static bool
samePoint(const Point & a, const Point & b, const Real length_tol = 1e-12)
{
  return (a - b).norm() <= length_tol;
}

static void
addUniquePoint(std::vector<Point> & points, const Point & point, const Real length_tol = 1e-12)
{
  for (const auto & existing_point : points)
    if (samePoint(existing_point, point, length_tol))
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

static bool
pointInPolygon2D(const Point & point,
                 const std::vector<Point> & polygon,
                 const Real length_tol,
                 const Real area_tol)
{
  if (polygon.size() < 3)
    return false;

  bool inside = false;

  for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++)
  {
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

  for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++)
    area += polygon[j](0) * polygon[i](1) - polygon[i](0) * polygon[j](1);

  return 0.5 * area;
}

struct BoundarySegment
{
  dof_id_type node0;
  dof_id_type node1;
  Point p0;
  Point p1;
};

static std::pair<dof_id_type, dof_id_type>
edgeKey(const dof_id_type node0, const dof_id_type node1)
{
  return {std::min(node0, node1), std::max(node0, node1)};
}

static bool
hasNonzeroArea3D(const std::vector<Point> & points, const Real tol = 1e-12)
{
  for (std::size_t i = 0; i < points.size(); ++i)
    for (std::size_t j = i + 1; j < points.size(); ++j)
      for (std::size_t k = j + 1; k < points.size(); ++k)
        if (((points[j] - points[i]).cross(points[k] - points[i])).norm() > tol)
          return true;

  return false;
}

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

static Point
faceNormal3D(const std::vector<Point> & points, const Real tol = 1e-12)
{
  for (std::size_t i = 0; i < points.size(); ++i)
    for (std::size_t j = i + 1; j < points.size(); ++j)
      for (std::size_t k = j + 1; k < points.size(); ++k)
      {
        const Point normal = (points[j] - points[i]).cross(points[k] - points[i]);

        if (normal.norm() > tol)
          return normal;
      }

  return Point();
}

static Real
tetVolume6(const Point & a, const Point & b, const Point & c, const Point & d)
{
  return (b - a).cross(c - a) * (d - a);
}

static bool
findConcaveEdge3D(const std::vector<std::vector<Point>> & side_points,
                  const std::vector<Point> & unique_points,
                  std::pair<unsigned int, unsigned int> & concave_edge,
                  const Real tol = 1e-12)
{
  if (side_points.size() < 4 || unique_points.size() < 4)
    return false;

  Point polyhedron_center;

  for (const auto & point : unique_points)
    polyhedron_center += point;

  polyhedron_center /= unique_points.size();

  const auto pointIndex = [&](const Point & point, unsigned int & index) -> bool
  {
    for (const auto i : make_range(unique_points.size()))
      if (samePoint(unique_points[i], point))
      {
        index = cast_int<unsigned int>(i);
        return true;
      }

    return false;
  };

  std::vector<Point> face_centers(side_points.size());

  for (const auto side_i : make_range(side_points.size()))
  {
    for (const auto & point : side_points[side_i])
      face_centers[side_i] += point;

    face_centers[side_i] /= side_points[side_i].size();

    if (faceNormal3D(side_points[side_i], tol).norm() <= tol)
      return false;
  }

  std::map<std::pair<unsigned int, unsigned int>, std::vector<unsigned int>> edge_to_sides;

  for (const auto side_i : make_range(side_points.size()))
    for (const auto point_i : make_range(side_points[side_i].size()))
    {
      unsigned int p0 = 0;
      unsigned int p1 = 0;

      if (!pointIndex(side_points[side_i][point_i], p0) ||
          !pointIndex(side_points[side_i][(point_i + 1) % side_points[side_i].size()], p1))
        return false;

      edge_to_sides[{std::min(p0, p1), std::max(p0, p1)}].push_back(
          cast_int<unsigned int>(side_i));
    }

  bool found_concave_edge = false;
  Real best_concave_score = 0.0;

  for (const auto & edge_sides : edge_to_sides)
  {
    const auto & adjacent_sides = edge_sides.second;

    if (adjacent_sides.size() != 2)
      continue;

    const unsigned int side0 = adjacent_sides[0];
    const unsigned int side1 = adjacent_sides[1];
    const Point edge_point0 = unique_points[edge_sides.first.first];
    const Point edge_point1 = unique_points[edge_sides.first.second];
    const Point edge_vector = edge_point1 - edge_point0;

    if (edge_vector.norm() <= tol)
      continue;

    const Point edge_axis = edge_vector / edge_vector.norm();
    const Point edge_midpoint = 0.5 * (edge_point0 + edge_point1);

    const auto edgeRadialDirection = [&](const Point & point) -> Point
    {
      const Point offset = point - edge_midpoint;
      return offset - (offset * edge_axis) * edge_axis;
    };

    Point side_direction0 = edgeRadialDirection(face_centers[side0]);
    Point side_direction1 = edgeRadialDirection(face_centers[side1]);

    if (side_direction0.norm() <= tol || side_direction1.norm() <= tol)
      continue;

    side_direction0 /= side_direction0.norm();
    side_direction1 /= side_direction1.norm();

    const auto signedAngle = [&](const Point & from, const Point & to) -> Real
    { return std::atan2(edge_axis * from.cross(to), from * to); };

    const Real side_angle = signedAngle(side_direction0, side_direction1);

    if (std::abs(side_angle) <= 1e-10 ||
        std::abs(std::abs(side_angle) - libMesh::pi) <= 1e-10)
      continue;

    const auto inSmallerWedge = [&](Point direction) -> bool
    {
      if (direction.norm() <= tol)
        return true;

      direction /= direction.norm();

      const Real angle = signedAngle(side_direction0, direction);

      if (side_angle > 0.0)
        return angle >= -1e-10 && angle <= side_angle + 1e-10;
      else
        return angle <= 1e-10 && angle >= side_angle - 1e-10;
    };

    unsigned int inside_count = 0;
    unsigned int outside_count = 0;

    for (const auto & point : unique_points)
    {
      if (samePoint(point, edge_point0) || samePoint(point, edge_point1))
        continue;

      const Point direction = edgeRadialDirection(point);

      if (direction.norm() <= tol)
        continue;

      if (inSmallerWedge(direction))
        ++inside_count;
      else
        ++outside_count;
    }

    const Point center_direction = edgeRadialDirection(polyhedron_center);
    const bool center_outside =
        center_direction.norm() > tol && !inSmallerWedge(center_direction);

    if (!center_outside && outside_count <= inside_count)
      continue;

    const Real concave_score =
        (center_outside ? 1000.0 : 0.0) + static_cast<Real>(outside_count) -
        static_cast<Real>(inside_count);

    if (concave_score > best_concave_score)
    {
      found_concave_edge = true;
      best_concave_score = concave_score;
      concave_edge = edge_sides.first;
    }
  }

  return found_concave_edge;
}

static bool
isConvexPolyhedron3D(const std::vector<std::vector<Point>> & sides, const Real tol = 1e-10)
{
  if (sides.size() < 4)
    return false;

  std::vector<Point> points;

  for (const auto & side : sides)
    for (const auto & point : side)
      addUniquePoint(points, point, tol);

  if (points.size() < 4)
    return false;

  Point center;

  for (const auto & point : points)
    center += point;

  center /= points.size();

  for (const auto & side : sides)
  {
    if (side.size() < 3)
      return false;

    const Point normal = faceNormal3D(side, tol);

    if (normal.norm() <= tol)
      return false;

    Point outward_normal = normal;

    if (normal * (center - side[0]) > tol)
      outward_normal = -1.0 * normal;

    for (const auto & point : points)
    {
      const Point offset = point - side[0];

      if (outward_normal * offset > tol * std::max(Real(1.0), offset.norm()))
        return false;
    }
  }

  return true;
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

static void
addTriangulatedSidePoints3D(std::vector<std::vector<Point>> & sides,
                            const std::vector<Point> & side_points,
                            const Real tol = 1e-12)
{
  std::vector<Point> unique_side_points;

  for (const auto & point : side_points)
    addUniquePoint(unique_side_points, point, tol);

  if (unique_side_points.size() < 3 || !hasNonzeroArea3D(unique_side_points, tol))
    return;

  const Point normal = faceNormal3D(unique_side_points, tol);
  std::vector<Point> sorted_side_points =
      normal.norm() > tol ? sortPointsAroundAxis3D(unique_side_points, normal, tol)
                          : unique_side_points;

  if (sorted_side_points.size() == 3)
  {
    sides.push_back(sorted_side_points);
    return;
  }

  for (std::size_t i = 1; i + 1 < sorted_side_points.size(); ++i)
  {
    const std::vector<Point> triangle_points = {
        sorted_side_points[0], sorted_side_points[i], sorted_side_points[i + 1]};

    if (hasNonzeroArea3D(triangle_points, tol))
      sides.push_back(triangle_points);
  }
}

std::unique_ptr<MeshBase>
DualMeshGenerator::generate()
{
  const auto input_mesh = std::move(_input);

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

    std::size_t split_nonconvex_polyhedron_count = 0;
    std::size_t tetrahedralized_nonconvex_polyhedron_count = 0;
    std::size_t skipped_nonconvex_polyhedron_count = 0;

    const auto tryAddPolyhedron =
        [&](MeshBase & mesh, const std::vector<std::vector<Point>> & polyhedron_side_points) -> bool
    {
      if (polyhedron_side_points.size() < 4)
        return false;

      std::vector<Node *> local_nodes;
      std::vector<std::shared_ptr<libMesh::Polygon>> sides;

      const auto getLocalNode = [&](const Point & point)
      {
        for (auto * const node : local_nodes)
          if (samePoint(*node, point))
            return node;

        Node * const node = mesh.add_point(point);
        local_nodes.push_back(node);

        return node;
      };

      sides.reserve(polyhedron_side_points.size());

      for (const auto & side_points : polyhedron_side_points)
      {
        auto side = std::make_shared<libMesh::C0Polygon>(side_points.size());

        for (const auto i : make_range(side_points.size()))
          side->set_node(i, getLocalNode(side_points[i]));

        sides.push_back(side);
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

    const auto addPolyhedron =
        [&](const std::vector<std::vector<Point>> & polyhedron_side_points) -> bool
    {
      if (polyhedron_side_points.size() < 4)
        return false;

      auto trial_mesh = buildReplicatedMesh(3);

      if (!tryAddPolyhedron(*trial_mesh, polyhedron_side_points))
        return false;

      return tryAddPolyhedron(*dualMesh, polyhedron_side_points);
    };

    const auto addTetrahedralizedPolyhedron =
        [&](const std::vector<std::vector<Point>> & side_points,
            const std::vector<Point> & body_centroid_points) -> bool
    {
      std::vector<std::vector<Point>> surface_triangles;

      for (const auto & side : side_points)
        addTriangulatedSidePoints3D(surface_triangles, side);

      if (surface_triangles.size() < 4)
        return false;

      std::vector<Point> unique_points;

      for (const auto & triangle : surface_triangles)
        for (const auto & point : triangle)
          addUniquePoint(unique_points, point);

      Point interior_point;
      std::vector<Point> interior_point_candidates;

      for (const auto & body_centroid_point : body_centroid_points)
        for (const auto & point : unique_points)
          if (samePoint(body_centroid_point, point))
          {
            addUniquePoint(interior_point_candidates, point);
            break;
          }

      const auto & point_source =
          interior_point_candidates.size() >= 2 ? interior_point_candidates : unique_points;

      for (const auto & point : point_source)
        interior_point += point;

      interior_point /= point_source.size();

      bool guard_concave_plane = false;
      Point concave_plane_normal;
      std::pair<unsigned int, unsigned int> concave_edge;

      if (findConcaveEdge3D(side_points, unique_points, concave_edge))
      {
        const Point edge_point0 = unique_points[concave_edge.first];
        const Point edge_point1 = unique_points[concave_edge.second];
        concave_plane_normal = (edge_point1 - edge_point0).cross(interior_point - edge_point0);

        if (concave_plane_normal.norm() > 1e-12)
        {
          concave_plane_normal /= concave_plane_normal.norm();
          guard_concave_plane = true;
        }
      }

      std::vector<std::array<Point, 4>> tets;

      for (const auto & triangle : surface_triangles)
      {
        const Real volume6 = tetVolume6(interior_point, triangle[0], triangle[1], triangle[2]);

        if (std::abs(volume6) <= 1e-12)
          continue;

        if (guard_concave_plane)
        {
          unsigned int positive_count = 0;
          unsigned int negative_count = 0;

          for (const auto & point : triangle)
          {
            const Real signed_distance = concave_plane_normal * (point - interior_point);

            if (signed_distance > 1e-10)
              ++positive_count;
            else if (signed_distance < -1e-10)
              ++negative_count;
          }

          if (positive_count > 0 && negative_count > 0)
            continue;
        }

        if (volume6 > 0.0)
          tets.push_back({interior_point, triangle[0], triangle[1], triangle[2]});
        else
          tets.push_back({interior_point, triangle[0], triangle[2], triangle[1]});
      }

      if (tets.empty())
        return false;

      std::vector<Node *> local_nodes;

      const auto getLocalNode = [&](const Point & point)
      {
        for (auto * const node : local_nodes)
          if (samePoint(*node, point))
            return node;

        Node * const node = dualMesh->add_point(point);
        local_nodes.push_back(node);

        return node;
      };

      for (const auto & tet_points : tets)
      {
        auto tet = std::make_unique<Tet4>();

        for (const auto i : make_range(tet_points.size()))
          tet->set_node(i, getLocalNode(tet_points[i]));

        dualMesh->add_elem(std::move(tet));
      }

      return true;
    };

    const auto splitConcaveEdgePolyhedron =
        [&](const std::vector<std::vector<Point>> & side_points,
            const std::vector<Point> & body_centroid_points,
            std::vector<std::vector<Point>> & positive_side_points,
            std::vector<std::vector<Point>> & negative_side_points) -> bool
    {
      if (side_points.size() < 4 || body_centroid_points.empty())
        return false;

      positive_side_points.clear();
      negative_side_points.clear();

      std::vector<Point> unique_points;

      for (const auto & side : side_points)
        for (const auto & point : side)
          addUniquePoint(unique_points, point);

      if (unique_points.size() < 4)
        return false;

      Point polyhedron_center;

      for (const auto & point : unique_points)
        polyhedron_center += point;

      polyhedron_center /= unique_points.size();

      const auto pointIndex = [&](const Point & point) -> unsigned int
      {
        for (const auto i : make_range(unique_points.size()))
          if (samePoint(unique_points[i], point))
            return cast_int<unsigned int>(i);

        mooseError("Could not find point while splitting non-convex 3D dual polyhedron.");
        return 0;
      };

      std::vector<Point> face_centers(side_points.size());

      for (const auto side_i : make_range(side_points.size()))
      {
        for (const auto & point : side_points[side_i])
          face_centers[side_i] += point;

        face_centers[side_i] /= side_points[side_i].size();

        Point normal = faceNormal3D(side_points[side_i]);

        if (normal.norm() <= 1e-12)
          return false;
      }

      std::map<std::pair<unsigned int, unsigned int>, std::vector<unsigned int>> edge_to_sides;

      for (const auto side_i : make_range(side_points.size()))
        for (const auto point_i : make_range(side_points[side_i].size()))
        {
          const unsigned int p0 = pointIndex(side_points[side_i][point_i]);
          const unsigned int p1 =
              pointIndex(side_points[side_i][(point_i + 1) % side_points[side_i].size()]);

          edge_to_sides[{std::min(p0, p1), std::max(p0, p1)}].push_back(
              cast_int<unsigned int>(side_i));
        }

      bool found_concave_edge = false;
      std::pair<unsigned int, unsigned int> concave_edge;
      Real best_concave_score = 0.0;

      for (const auto & edge_sides : edge_to_sides)
      {
        const auto & adjacent_sides = edge_sides.second;

        if (adjacent_sides.size() != 2)
          continue;

        const unsigned int side0 = adjacent_sides[0];
        const unsigned int side1 = adjacent_sides[1];
        const Point edge_point0 = unique_points[edge_sides.first.first];
        const Point edge_point1 = unique_points[edge_sides.first.second];
        const Point edge_vector = edge_point1 - edge_point0;

        if (edge_vector.norm() <= 1e-12)
          continue;

        const Point edge_axis = edge_vector / edge_vector.norm();
        const Point edge_midpoint = 0.5 * (edge_point0 + edge_point1);

        const auto edgeRadialDirection = [&](const Point & point) -> Point
        {
          const Point offset = point - edge_midpoint;
          return offset - (offset * edge_axis) * edge_axis;
        };

        Point side_direction0 = edgeRadialDirection(face_centers[side0]);
        Point side_direction1 = edgeRadialDirection(face_centers[side1]);

        if (side_direction0.norm() <= 1e-12 || side_direction1.norm() <= 1e-12)
          continue;

        side_direction0 /= side_direction0.norm();
        side_direction1 /= side_direction1.norm();

        const auto signedAngle = [&](const Point & from, const Point & to) -> Real
        { return std::atan2(edge_axis * from.cross(to), from * to); };

        const Real side_angle = signedAngle(side_direction0, side_direction1);

        if (std::abs(side_angle) <= 1e-10 ||
            std::abs(std::abs(side_angle) - libMesh::pi) <= 1e-10)
          continue;

        const auto inSmallerWedge = [&](Point direction) -> bool
        {
          if (direction.norm() <= 1e-12)
            return true;

          direction /= direction.norm();

          const Real angle = signedAngle(side_direction0, direction);

          if (side_angle > 0.0)
            return angle >= -1e-10 && angle <= side_angle + 1e-10;
          else
            return angle <= 1e-10 && angle >= side_angle - 1e-10;
        };

        unsigned int inside_count = 0;
        unsigned int outside_count = 0;

        for (const auto & point : unique_points)
        {
          if (samePoint(point, edge_point0) || samePoint(point, edge_point1))
            continue;

          const Point direction = edgeRadialDirection(point);

          if (direction.norm() <= 1e-12)
            continue;

          if (inSmallerWedge(direction))
            ++inside_count;
          else
            ++outside_count;
        }

        const Point center_direction = edgeRadialDirection(polyhedron_center);
        const bool center_outside =
            center_direction.norm() > 1e-12 && !inSmallerWedge(center_direction);

        if (!center_outside && outside_count <= inside_count)
          continue;

        const Real concave_score =
            (center_outside ? 1000.0 : 0.0) + static_cast<Real>(outside_count) -
            static_cast<Real>(inside_count);

        if (concave_score > best_concave_score)
        {
          found_concave_edge = true;
          best_concave_score = concave_score;
          concave_edge = edge_sides.first;
        }
      }

      if (!found_concave_edge)
        return false;

      const Point edge_point0 = unique_points[concave_edge.first];
      const Point edge_point1 = unique_points[concave_edge.second];
      const Point edge_vector = edge_point1 - edge_point0;

      if (edge_vector.norm() <= 1e-12)
        return false;

      _console << "DualMeshGenerator 3D concave line nodes:\n"
               << "  " << edge_point0 << "\n"
               << "  " << edge_point1 << "\n";

      const auto distanceToConcaveEdge = [&](const Point & point) -> Real
      {
        return ((point - edge_point0).cross(edge_vector)).norm() / edge_vector.norm();
      };

      std::vector<Point> body_points;

      for (const auto & point : body_centroid_points)
        for (const auto & unique_point : unique_points)
          if (samePoint(point, unique_point) && !samePoint(unique_point, edge_point0) &&
              !samePoint(unique_point, edge_point1))
          {
            addUniquePoint(body_points, unique_point);
            break;
          }

      if (body_points.size() < 2)
        return false;

      std::sort(body_points.begin(),
                body_points.end(),
                [&distanceToConcaveEdge](const Point & a, const Point & b)
                { return distanceToConcaveEdge(a) < distanceToConcaveEdge(b); });

      bool printed_first_body_pair = false;

      for (std::size_t body_i = 0; body_i < body_points.size(); ++body_i)
        for (std::size_t body_j = body_i + 1; body_j < body_points.size(); ++body_j)
        {
          const Point & body_point0 = body_points[body_i];
          const Point & body_point1 = body_points[body_j];
          const unsigned int body_point0_index = pointIndex(body_point0);
          const unsigned int body_point1_index = pointIndex(body_point1);
          const auto body_edge = std::make_pair(std::min(body_point0_index, body_point1_index),
                                                std::max(body_point0_index, body_point1_index));

          if (edge_to_sides.find(body_edge) == edge_to_sides.end())
            continue;

          if (!printed_first_body_pair)
          {
            _console << "DualMeshGenerator first 3D concave split body centroid pair:\n"
                     << "  " << body_point0 << "\n"
                     << "  " << body_point1 << "\n";
            printed_first_body_pair = true;
          }

          const auto concave_edge_sides_it = edge_to_sides.find(concave_edge);
          const auto body_edge_sides_it = edge_to_sides.find(body_edge);

          if (concave_edge_sides_it == edge_to_sides.end() ||
              body_edge_sides_it == edge_to_sides.end() ||
              concave_edge_sides_it->second.size() != 2 || body_edge_sides_it->second.size() != 2)
            continue;

          const auto isCutEdge = [&](const std::pair<unsigned int, unsigned int> & edge)
          { return edge == concave_edge || edge == body_edge; };

          const bool first_diagonal =
              (edge_point0 - body_point1).norm() <= (edge_point1 - body_point0).norm();
          std::vector<std::vector<Point>> split_faces;

          if (first_diagonal)
            split_faces = {{edge_point0, edge_point1, body_point1},
                           {edge_point0, body_point1, body_point0}};
          else
            split_faces = {{edge_point0, edge_point1, body_point0},
                           {edge_point1, body_point1, body_point0}};

          std::vector<bool> positive_piece_side(side_points.size(), false);
          std::vector<unsigned int> side_stack = {concave_edge_sides_it->second[0]};
          positive_piece_side[side_stack.back()] = true;

          while (!side_stack.empty())
          {
            const unsigned int side_i = side_stack.back();
            side_stack.pop_back();

            for (const auto point_i : make_range(side_points[side_i].size()))
            {
              const unsigned int point0 = pointIndex(side_points[side_i][point_i]);
              const unsigned int point1 =
                  pointIndex(side_points[side_i][(point_i + 1) % side_points[side_i].size()]);
              const auto side_edge = std::make_pair(std::min(point0, point1),
                                                     std::max(point0, point1));

              if (isCutEdge(side_edge))
                continue;

              const auto adjacent_sides_it = edge_to_sides.find(side_edge);

              if (adjacent_sides_it == edge_to_sides.end())
                continue;

              for (const auto adjacent_side : adjacent_sides_it->second)
                if (!positive_piece_side[adjacent_side])
                {
                  positive_piece_side[adjacent_side] = true;
                  side_stack.push_back(adjacent_side);
                }
            }
          }

          if (positive_piece_side[concave_edge_sides_it->second[1]] ||
              positive_piece_side[body_edge_sides_it->second[0]] ==
                  positive_piece_side[body_edge_sides_it->second[1]])
            continue;

          std::vector<std::vector<Point>> positive_piece_side_points;
          std::vector<std::vector<Point>> negative_piece_side_points;

          for (const auto side_i : make_range(side_points.size()))
            if (positive_piece_side[side_i])
              addSidePoints3D(positive_piece_side_points, side_points[side_i]);
            else
              addSidePoints3D(negative_piece_side_points, side_points[side_i]);

          for (auto split_face : split_faces)
          {
            addSidePoints3D(positive_piece_side_points, split_face);
            std::reverse(split_face.begin(), split_face.end());
            addSidePoints3D(negative_piece_side_points, split_face);
          }

          if (positive_piece_side_points.size() < 4 || negative_piece_side_points.size() < 4)
            continue;

          positive_side_points = std::move(positive_piece_side_points);
          negative_side_points = std::move(negative_piece_side_points);

          return true;
        }

      return false;
    };

    const auto collectConvexSplitPolyhedra =
        [&](const auto & self,
            const std::vector<std::vector<Point>> & side_points,
            const std::vector<Point> & body_centroid_points,
            std::vector<std::vector<std::vector<Point>>> & convex_pieces,
            const unsigned int depth) -> bool
    {
      if (side_points.size() < 4 || depth > 32)
        return false;

      auto trial_mesh = buildReplicatedMesh(3);

      if (isConvexPolyhedron3D(side_points) && tryAddPolyhedron(*trial_mesh, side_points))
      {
        convex_pieces.push_back(side_points);
        return true;
      }

      std::vector<std::vector<Point>> positive_side_points;
      std::vector<std::vector<Point>> negative_side_points;

      if (!splitConcaveEdgePolyhedron(
              side_points, body_centroid_points, positive_side_points, negative_side_points))
        return false;

      const auto childBodyCentroids =
          [&](const std::vector<std::vector<Point>> & child_side_points)
      {
        std::vector<Point> child_points;
        std::vector<Point> child_body_centroids;

        for (const auto & side : child_side_points)
          for (const auto & point : side)
            addUniquePoint(child_points, point);

        for (const auto & body_centroid_point : body_centroid_points)
          for (const auto & child_point : child_points)
            if (samePoint(body_centroid_point, child_point))
            {
              addUniquePoint(child_body_centroids, body_centroid_point);
              break;
            }

        return child_body_centroids;
      };

      std::vector<std::vector<std::vector<Point>>> positive_convex_pieces;
      std::vector<std::vector<std::vector<Point>>> negative_convex_pieces;

      if (!self(self,
                positive_side_points,
                childBodyCentroids(positive_side_points),
                positive_convex_pieces,
                depth + 1) ||
          !self(self,
                negative_side_points,
                childBodyCentroids(negative_side_points),
                negative_convex_pieces,
                depth + 1))
        return false;

      convex_pieces.insert(
          convex_pieces.end(), positive_convex_pieces.begin(), positive_convex_pieces.end());
      convex_pieces.insert(
          convex_pieces.end(), negative_convex_pieces.begin(), negative_convex_pieces.end());

      ++split_nonconvex_polyhedron_count;
      return true;
    };

    const auto addPolyhedronOrSplit = [&](const std::vector<std::vector<Point>> & side_points,
                                          const std::vector<Point> & body_centroid_points)
    {
      std::vector<std::vector<std::vector<Point>>> convex_pieces;

      if (!collectConvexSplitPolyhedra(
              collectConvexSplitPolyhedra, side_points, body_centroid_points, convex_pieces, 0))
      {
        if (addTetrahedralizedPolyhedron(side_points, body_centroid_points))
        {
          ++tetrahedralized_nonconvex_polyhedron_count;
          return true;
        }

        ++skipped_nonconvex_polyhedron_count;
        return false;
      }

      for (const auto & convex_piece : convex_pieces)
      {
        if (!addPolyhedron(convex_piece))
        {
          if (addTetrahedralizedPolyhedron(convex_piece, body_centroid_points))
          {
            ++tetrahedralized_nonconvex_polyhedron_count;
            continue;
          }

          ++skipped_nonconvex_polyhedron_count;
          return false;
        }
      }

      return true;
    };

    // Build one dual polyhedron around each primal node using element centroids, exterior face
    // centroids, and primal boundary vertices.
    for (const auto & node_elems : source_node_to_elems)
    {
      const dof_id_type source_node_id = node_elems.first;
      const Point & source_point = *input_mesh->node_ptr(source_node_id);
      std::map<std::pair<dof_id_type, dof_id_type>, std::vector<Point>> edge_to_points;
      std::map<std::pair<dof_id_type, dof_id_type>, std::vector<Point>>
          midpoint_boundary_face_centroids;
      std::vector<Point> body_centroid_points;
      std::vector<Point> boundary_face_centroids;
      Point boundary_normal;

      for (const auto & elem : node_elems.second)
      {
        const Point elem_centroid = elem->true_centroid();
        addUniquePoint(body_centroid_points, elem_centroid);

        for (const auto side : elem->side_index_range())
        {
          auto side_elem = elem->build_side_ptr(side);
          const Point face_centroid = side_elem->true_centroid();

          std::vector<dof_id_type> current_side_node_ids;
          std::vector<Point> current_side_points;
          current_side_node_ids.reserve(side_elem->n_vertices());
          current_side_points.reserve(side_elem->n_vertices());

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

        for (std::size_t i = 0; i < midpoint_face_centroids.second.size(); ++i)
          for (std::size_t j = i + 1; j < midpoint_face_centroids.second.size(); ++j)
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
            if ((samePoint(point0, split_boundary_face.first) &&
                 samePoint(point1, split_boundary_face.second)) ||
                (samePoint(point0, split_boundary_face.second) &&
                 samePoint(point1, split_boundary_face.first)))
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
            for (std::size_t i = 0; i < sorted_boundary_points.size(); ++i)
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

      addPolyhedronOrSplit(polyhedron_side_points, body_centroid_points);
    }

    _console << "DualMeshGenerator split " << split_nonconvex_polyhedron_count
             << " non-convex 3D dual polyhedra.\n"
             << "DualMeshGenerator tetrahedralized " << tetrahedralized_nonconvex_polyhedron_count
             << " non-convex 3D dual polyhedra.\n"
             << "DualMeshGenerator skipped " << skipped_nonconvex_polyhedron_count
             << " non-convex 3D dual polyhedra.\n";

    std::map<std::string, std::size_t> elem_type_counts;

    for (const auto & elem : dualMesh->element_ptr_range())
      ++elem_type_counts[Moose::stringify(elem->type())];

    _console << "DualMeshGenerator output element types:";

    if (elem_type_counts.empty())
      _console << " none";
    else
      for (const auto & elem_type_count : elem_type_counts)
        _console << " " << elem_type_count.first << "=" << elem_type_count.second;

    _console << "\n" << std::flush;

    dualMesh->unset_is_prepared();
    return dynamic_pointer_cast<MeshBase>(dualMesh);
  }

  // BEGIN 2D
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

  for (const auto & elem : input_mesh->element_ptr_range())
  {
    for (const auto side : elem->side_index_range())
    {
      if (elem->neighbor_ptr(side) == nullptr)
      {
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

  for (std::size_t i = 0; i < physical_boundary_segments.size(); ++i)
  {
    boundary_node_to_segments[physical_boundary_segments[i].node0].push_back(i);
    boundary_node_to_segments[physical_boundary_segments[i].node1].push_back(i);
  }

  const auto otherBoundaryNode = [&](const BoundarySegment & segment,
                                     const dof_id_type node_id) -> dof_id_type
  { return segment.node0 == node_id ? segment.node1 : segment.node0; };

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
      boundary_vertex_nodes.insert(node_id);
  }

  std::vector<std::pair<Point, Point>> boundary_clip_segments;
  std::vector<bool> used_boundary_segments(physical_boundary_segments.size(), false);

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

  for (const auto boundary_vertex_node : boundary_vertex_nodes)
  {
    const auto node_segment_it = boundary_node_to_segments.find(boundary_vertex_node);

    if (node_segment_it == boundary_node_to_segments.end())
      continue;

    for (const auto segment_id : node_segment_it->second)
      if (!used_boundary_segments[segment_id])
        traceBoundaryClipSegment(boundary_vertex_node, segment_id);
  }

  for (std::size_t segment_id = 0; segment_id < physical_boundary_segments.size(); ++segment_id)
    if (!used_boundary_segments[segment_id])
      traceBoundaryClipSegment(physical_boundary_segments[segment_id].node0, segment_id);

  std::vector<Point> dual_centers;
  std::unordered_map<dof_id_type, dof_id_type> source_elem_to_center_id;
  std::unordered_map<dof_id_type, std::vector<const Elem *>> source_node_to_elems;
  std::unique_ptr<ReplicatedMesh> tri_mesh;

  if (use_voronoi)
  {
    tri_mesh = buildReplicatedMesh(2);

    for (const auto & node : input_mesh->node_ptr_range())
    {
      Node * new_node = tri_mesh->add_point(*node);

      auto node_elem = std::make_unique<NodeElem>();
      node_elem->set_node(0) = new_node;
      tri_mesh->add_elem(std::move(node_elem));
    }

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

    triangulator.insert_extra_points() = false;
    triangulator.smooth_after_generating() = false;

    triangulator.triangulate();

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
  else
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

  auto dualMesh = buildReplicatedMesh(2);

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

  const auto clipDualPolygonToBoundary = [&](const std::vector<Point> & dual_points)
  {
    std::vector<Point> clipped_points;

    if (dual_points.size() < 3)
      return clipped_points;

    for (const auto & point : dual_points)
      if (pointInsideBoundary(point))
        addUniquePoint(clipped_points, point, length_tol);

    for (unsigned int i = 0; i < dual_points.size(); ++i)
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
        samePoint(unique_clipped_points.front(), unique_clipped_points.back(), length_tol))
      unique_clipped_points.pop_back();

    return unique_clipped_points;
  };

  const auto isBoundaryVertexPoint = [&](const Point & point)
  {
    for (const auto boundary_vertex_node : boundary_vertex_nodes)
    {
      const auto point_it = boundary_node_points.find(boundary_vertex_node);

      if (point_it != boundary_node_points.end() && samePoint(point, point_it->second, length_tol))
        return true;
    }

    return false;
  };

  const auto isBoundarySegmentPoint = [&](const Point & point)
  {
    for (const auto & boundary_segment : boundary_clip_segments)
      if (pointOnSegment2D(
              point, boundary_segment.first, boundary_segment.second, length_tol, area_tol))
        return true;

    return false;
  };

  const auto concaveBoundaryVertexIndex = [&](const std::vector<Point> & points) -> std::size_t
  {
    if (points.size() < 4)
      return points.size();

    const Real signed_area = polygonSignedArea2D(points);

    if (std::abs(signed_area) < area_tol)
      return points.size();

    const Real orientation = signed_area > 0.0 ? 1.0 : -1.0;

    for (std::size_t i = 0; i < points.size(); ++i)
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

    for (unsigned int i = 0; i < points.size(); ++i)
      dual_elem->set_node(i, dualMesh->add_point(points[i]));

    // Fixing flipped elems
    if (dual_elem->is_flipped())
    {
      auto reversed_elem = std::make_unique<libMesh::C0Polygon>(points.size());

      for (unsigned int i = 0; i < points.size(); ++i)
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

    for (unsigned int i = 0; i < incident_elems.size(); ++i)
      for (unsigned int j = i + 1; j < incident_elems.size(); ++j)
        if (elementsShareTwoNodes(incident_elems[i], incident_elems[j]))
        {
          adjacency[i].push_back(j);
          adjacency[j].push_back(i);
        }

    // Start at an endpoint for boundary chains, otherwise start anywhere.
    unsigned int start = 0;

    for (unsigned int i = 0; i < adjacency.size(); ++i)
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
          addUniquePoint(dual_points, boundary_point_it->second, length_tol);

        for (const auto & midpoint : boundary_midpoint_it->second)
          addUniquePoint(dual_points, midpoint, length_tol);

        if (boundary_point_it != boundary_node_points.end() && !source_center_points.empty())
        {
          Point center_average;

          for (const auto & center_point : source_center_points)
            center_average += center_point;

          center_average /= source_center_points.size();

          const Point sort_center = 0.5 * (boundary_point_it->second + center_average);

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
      for (std::size_t i = 0; i < dual_points.size(); ++i)
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

      for (std::size_t i = 0; i < sorted_points.size(); ++i)
      {
        if (!isBoundarySegmentPoint(sorted_points[i].first))
          continue;

        const std::size_t next_i = (i + 1) % sorted_points.size();
        const std::size_t prev_i = (i + sorted_points.size() - 1) % sorted_points.size();

        // We pick the direction of fan-triangulating concave polygons such that we never create a
        // triangle that bridges across a boundary
        if (!isBoundarySegmentPoint(sorted_points[next_i].first))
        {
          for (std::size_t k = 0; k < sorted_points.size(); ++k)
            fan_points.push_back(sorted_points[(i + k) % sorted_points.size()].first);

          break;
        }

        if (!isBoundarySegmentPoint(sorted_points[prev_i].first))
        {
          for (std::size_t k = 0; k < sorted_points.size(); ++k)
            fan_points.push_back(
                sorted_points[(i + sorted_points.size() - k) % sorted_points.size()].first);

          break;
        }
      }

      for (std::size_t i = 1; i + 1 < fan_points.size(); ++i)
      {
        const std::vector<Point> triangle_points = {
            fan_points[0], fan_points[i], fan_points[i + 1]};

        if (std::abs(cross2D(triangle_points[0], triangle_points[1], triangle_points[2])) >
            area_tol)
          addDualElement(triangle_points);
      }

      continue;
    }

    addDualElement(dual_points);
  }

  dualMesh->unset_is_prepared();

  return dynamic_pointer_cast<MeshBase>(dualMesh);
}
