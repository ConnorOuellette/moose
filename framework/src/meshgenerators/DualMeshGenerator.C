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
#include "MooseMeshUtils.h"
#include "libmesh/node_elem.h"
#include "libmesh/poly2tri_triangulator.h"
#include "libmesh/elem.h"

registerMooseObject("MooseApp", DualMeshGenerator);

InputParameters
DualMeshGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addClassDescription("Takes a 2D mesh as input and returns a Voronoi dual mesh, i.e.,"
                             "changes each input mode into an element and each input element "
                             "into a node located at its circumcenter.");

  params.addParam<Real>("boundary_node_angular_tol",
                        1,
                        "Tolerance (in degrees) for determining colinearity of boundary sides"
                        "when finding input mesh vertices.");
  params.addParam<Real>("boundary_edge_outside_tol",
                        1e-12,
                        "Tolerance (square of scalar distance) for determining whether polygon "
                        "vertices lie within or outside boundaries.");
  params.addParam<Real>(
      "dual_node_merge_tol", 1e-2, "Tolerance for merging nearly identical circumcenters");
  return params;
}

DualMeshGenerator::DualMeshGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _boundary_node_angular_tol(getParam<Real>("boundary_node_angular_tol")),
    _boundary_edge_outside_tol(getParam<Real>("boundary_edge_outside_tol")),
    _dual_node_merge_tol(getParam<Real>("dual_node_merge_tol"))
{
}

// Circumcenter method
Point
DualMeshGenerator::circumcenter(const Elem * elem)
{
  if (_boundary_edge_outside_tol < 0)
    mooseError("boundary_edge_outside_tol must be a positive value");

  const unsigned int n = elem->n_vertices();
  libmesh_assert_greater(n, 2);

  const Point & p0 = elem->point(0);

  Real A11 = 0.;
  Real A12 = 0.;
  Real A22 = 0.;

  Real b1 = 0.;
  Real b2 = 0.;

  for (unsigned int i = 1; i < n; ++i)
  {
    const Point & pi = elem->point(i);

    const Real dx = pi(0) - p0(0);
    const Real dy = pi(1) - p0(1);

    const Real rhs = 0.5 * (pi(0) * pi(0) + pi(1) * pi(1) - p0(0) * p0(0) - p0(1) * p0(1));

    A11 += dx * dx;
    A12 += dx * dy;
    A22 += dy * dy;

    b1 += dx * rhs;
    b2 += dy * rhs;
  }

  const Real det = A11 * A22 - A12 * A12;

  const Real cx = (A22 * b1 - A12 * b2) / det;
  const Real cy = (A11 * b2 - A12 * b1) / det;

  Point cc(cx, cy, 0.0);

  // Check whether circumcenter is inside polygon using ray casting.
  bool inside = false;

  for (unsigned int i = 0, j = n - 1; i < n; j = i++)
  {
    const Point & pi = elem->point(i);
    const Point & pj = elem->point(j);

    const bool intersects = ((pi(1) > cc(1)) != (pj(1) > cc(1))) &&
                            (cc(0) < (pj(0) - pi(0)) * (cc(1) - pi(1)) / (pj(1) - pi(1)) + pi(0));

    if (intersects)
      inside = !inside;
  }

  return cc;
}

static Real
cross2D(const Point & a, const Point & b, const Point & c)
{
  return (b(0) - a(0)) * (c(1) - a(1)) - (b(1) - a(1)) * (c(0) - a(0));
}

static Point
lineIntersection(const Point & p0, const Point & p1, const Point & q0, const Point & q1)
{
  const Point r = p1 - p0;
  const Point s = q1 - q0;

  const Real denom = r(0) * s(1) - r(1) * s(0);

  if (std::abs(denom) < 1e-14)
    return p1;

  const Point qp = q0 - p0;
  const Real t = (qp(0) * s(1) - qp(1) * s(0)) / denom;

  return p0 + t * r;
}

std::vector<Point>
DualMeshGenerator::clipPolygonToPhysicalBoundary(
    const std::vector<Point> & poly, const std::vector<std::pair<Point, Point>> & boundary_segments)
{
  std::vector<Point> output = poly;

  // Removing polygon vertices that end up outside of the boundary of the primal mesh
  for (const auto & segment : boundary_segments)
  {
    const Point & a = segment.first;
    const Point & b = segment.second;

    std::vector<Point> input = output;
    output.clear();

    if (input.empty())
      break;

    Point prev = input.back();
    bool prev_inside = cross2D(a, b, prev) >= -_boundary_edge_outside_tol;

    for (const Point & curr : input)
    {
      const bool curr_inside = cross2D(a, b, curr) >= -_boundary_edge_outside_tol;

      if (curr_inside)
      {
        if (!prev_inside)
          output.push_back(lineIntersection(prev, curr, a, b));

        output.push_back(curr);
      }
      else if (prev_inside)
        output.push_back(lineIntersection(prev, curr, a, b));

      prev = curr;
      prev_inside = curr_inside;
    }
  }

  return output;
}

// True only when two triangles share a full edge.
static bool
trianglesShareTwoNodes(const Elem * a, const Elem * b)
{
  unsigned int shared_nodes = 0;

  for (unsigned int i = 0; i < a->n_nodes(); ++i)
    for (unsigned int j = 0; j < b->n_nodes(); ++j)
      if (a->node_id(i) == b->node_id(j))
        ++shared_nodes;

  return shared_nodes == 2;
}

/// @brief
/// @return
std::unique_ptr<MeshBase>
DualMeshGenerator::generate()
{
  const auto input_mesh = std::move(_input);

  auto tri_mesh = buildReplicatedMesh(2);

  std::vector<std::pair<Point, Point>> physical_boundary_segments;
  std::unordered_map<dof_id_type, Node *> old_to_new_node;
  std::unordered_set<dof_id_type> real_node_ids;
  std::unordered_set<dof_id_type> boundary_node_ids;

  for (const auto & elem : input_mesh->element_ptr_range())
    for (const auto side : elem->side_index_range())
      if (elem->neighbor_ptr(side) == nullptr)
      {
        auto side_elem = elem->build_side_ptr(side);

        if (side_elem->n_nodes() == 2)
        {
          const auto old0 = side_elem->node_id(0);
          const auto old1 = side_elem->node_id(1);

          physical_boundary_segments.push_back({side_elem->point(0), side_elem->point(1)});
        }
      }

  for (const auto & node : input_mesh->node_ptr_range())
  {
    Node * new_node = tri_mesh->add_point(*node);

    old_to_new_node[node->id()] = new_node;
    real_node_ids.insert(new_node->id());

    auto node_elem = std::make_unique<NodeElem>();
    node_elem->set_node(0) = new_node;
    tri_mesh->add_elem(std::move(node_elem));
  }

  Node * p0 = tri_mesh->add_point(Point(-100.0, -100.0, 0.0));
  Node * p1 = tri_mesh->add_point(Point(100.0, -100.0, 0.0));
  Node * p2 = tri_mesh->add_point(Point(100.0, 100.0, 0.0));
  Node * p3 = tri_mesh->add_point(Point(-100.0, 100.0, 0.0));

  auto big_square = std::make_unique<Quad4>();

  big_square->set_node(0) = p0;
  big_square->set_node(1) = p1;
  big_square->set_node(2) = p2;
  big_square->set_node(3) = p3;

  if (big_square->is_flipped())
  {
    big_square->set_node(1) = p3;
    big_square->set_node(3) = p1;
  }

  tri_mesh->add_elem(std::move(big_square));

  Poly2TriTriangulator triangulator(dynamic_cast<UnstructuredMesh &>(*tri_mesh));
  triangulator.triangulation_type() = libMesh::TriangulatorInterface::PSLG;
  triangulator.minimum_angle() = 0;
  triangulator.desired_area() = 0;

  triangulator.insert_extra_points() = false;
  triangulator.smooth_after_generating() = false;

  triangulator.triangulate();

  struct CircumcenterRecord
  {
    Point point;
    std::vector<const Elem *> triangles;
  };

  std::vector<CircumcenterRecord> circumcenter_records;
  std::unordered_map<dof_id_type, dof_id_type> tri_elem_to_cc_id;
  std::unordered_map<dof_id_type, std::vector<const Elem *>> primal_node_to_triangles;

  // Compute and consolidate triangle circumcenters.
  for (const auto & tri_elem : tri_mesh->element_ptr_range())
  {
    if (tri_elem->n_vertices() != 3)
      continue;

    const Point cc = circumcenter(tri_elem);

    dof_id_type cc_id = libMesh::invalid_uint;

    for (const auto i : index_range(circumcenter_records))
      if ((circumcenter_records[i].point - cc).norm() <= _dual_node_merge_tol)
      {
        cc_id = i;
        break;
      }

    if (cc_id == libMesh::invalid_uint)
    {
      cc_id = circumcenter_records.size();
      circumcenter_records.push_back({cc, {}});
    }

    circumcenter_records[cc_id].triangles.push_back(tri_elem);
    tri_elem_to_cc_id[tri_elem->id()] = cc_id;

    for (const auto n : make_range(tri_elem->n_nodes()))
      primal_node_to_triangles[tri_elem->node_id(n)].push_back(tri_elem);
  }

  auto dualMesh = buildReplicatedMesh(2);

  // Build one dual element around each primal node.
  for (const auto & [primal_node_id, incident_tris] : primal_node_to_triangles)
  {
    // Connect incident triangles only if they share a full edge.
    std::vector<std::vector<unsigned int>> adjacency(incident_tris.size());

    for (unsigned int i = 0; i < incident_tris.size(); ++i)
      for (unsigned int j = i + 1; j < incident_tris.size(); ++j)
        if (trianglesShareTwoNodes(incident_tris[i], incident_tris[j]))
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

    std::vector<bool> used(incident_tris.size(), false);
    std::vector<dof_id_type> ordered_cc_ids;

    unsigned int current = start;
    unsigned int previous = libMesh::invalid_uint;

    // Walk triangle adjacency to collect circumcenters in connected order.
    while (true)
    {
      const Elem * tri = incident_tris[current];

      auto it = tri_elem_to_cc_id.find(tri->id());
      if (it != tri_elem_to_cc_id.end())
      {
        const dof_id_type cc_id = it->second;

        if (std::find(ordered_cc_ids.begin(), ordered_cc_ids.end(), cc_id) == ordered_cc_ids.end())
          ordered_cc_ids.push_back(cc_id);
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

    if (ordered_cc_ids.size() < 3)
      continue;

    Point dual_centroid;

    for (const auto cc_id : ordered_cc_ids)
      dual_centroid += circumcenter_records[cc_id].point;

    dual_centroid /= ordered_cc_ids.size();

    std::sort(ordered_cc_ids.begin(),
              ordered_cc_ids.end(),
              [&](const dof_id_type a, const dof_id_type b)
              {
                const Point & pa = circumcenter_records[a].point;
                const Point & pb = circumcenter_records[b].point;

                return std::atan2(pa(1) - dual_centroid(1), pa(0) - dual_centroid(0)) <
                       std::atan2(pb(1) - dual_centroid(1), pb(0) - dual_centroid(0));
              });

    auto dual_elem = std::make_unique<libMesh::C0Polygon>(ordered_cc_ids.size());

    for (unsigned int i = 0; i < ordered_cc_ids.size(); ++i)
      dual_elem->set_node(i, dualMesh->add_point(circumcenter_records[ordered_cc_ids[i]].point));

    if (dual_elem->is_flipped())
    {
      auto reversed_elem = std::make_unique<libMesh::C0Polygon>(ordered_cc_ids.size());

      for (unsigned int i = 0; i < ordered_cc_ids.size(); ++i)
        reversed_elem->set_node(
            i,
            dualMesh->add_point(
                circumcenter_records[ordered_cc_ids[ordered_cc_ids.size() - 1 - i]].point));

      dual_elem = std::move(reversed_elem);
    }

    if (!dual_elem->is_flipped())
      dualMesh->add_elem(std::move(dual_elem));
  }

  dualMesh->unset_is_prepared();
  tri_mesh->unset_is_prepared();

  return dynamic_pointer_cast<MeshBase>(dualMesh);
  //  return dynamic_pointer_cast<MeshBase>(tri_mesh);
}
