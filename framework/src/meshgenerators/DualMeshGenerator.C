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
  return params;
}

DualMeshGenerator::DualMeshGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _boundary_node_angular_tol(getParam<Real>("boundary_node_angular_tol"))
{
}

// Circumcenter method
Point
DualMeshGenerator::circumcenter(const Elem * elem)
{
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

  auto centroid = [&]() -> Point
  {
    Point c;

    for (unsigned int i = 0; i < n; ++i)
      c += elem->point(i);

    c /= n;

    return c;
  };

  const Real det = A11 * A22 - A12 * A12;

  // if (std::abs(det) < 1e-14)
  //   return centroid();

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

  if (!inside)
    return centroid();

  return cc;
}

std::unique_ptr<MeshBase>
DualMeshGenerator::generate()
{
  const auto mesh = std::move(_input);
  std::vector<libMesh::Point> circumcenters; // vector of all dual nodes
  circumcenters.reserve(mesh->n_elem());

  std::unordered_map<dof_id_type, std::vector<dof_id_type>> _node_to_elem_map;

  // Looping through primal elements, recording circumcenters, and creating map
  for (const auto & in_elem : mesh->element_ptr_range())
  {
    bool is_boundary_elem = false;

    for (const auto side : in_elem->side_index_range())
    {
      if (in_elem->neighbor_ptr(side) == nullptr)
      {
        is_boundary_elem = true;
        break;
      }
    }

    Point dual_point;

    if (is_boundary_elem) // Boundary elements are consistently degenerate/nonradial polygons,
                          // circumcenters are almost always not inside
      dual_point = in_elem->true_centroid();
    else
      dual_point = DualMeshGenerator::circumcenter(in_elem);

    circumcenters.push_back(dual_point);

    for (const auto n : make_range(in_elem->n_nodes()))
      _node_to_elem_map[in_elem->node_id(n)].push_back(in_elem->id());
  }

  auto dualMesh = buildReplicatedMesh(mesh->mesh_dimension());

  if (mesh->mesh_dimension() == 2) // Currently only supports 2D
  {

    for (const auto i : index_range(circumcenters))
    {
      dualMesh->add_point(circumcenters[i]);
    }

    std::unordered_map<dof_id_type, std::vector<dof_id_type>> _elem_to_node_map;
    for (const auto & [node_id, elements] : _node_to_elem_map)
    {
      for (dof_id_type elem_id : elements)
      {
        _elem_to_node_map[elem_id].push_back(node_id);
      }
    }

    // Get what nodes are boundary nodes
    // looping over primal elements....
    std::unordered_map<dof_id_type, std::vector<Point>> node_to_boundary_midpoints;
    std::unordered_map<dof_id_type, std::vector<dof_id_type>> node_to_boundary_neighbors;
    std::unordered_set<dof_id_type> corner_node_ids;
    std::unordered_set<dof_id_type> midpoint_node_ids;
    for (const auto & primalElem : mesh->element_ptr_range())
    {
      // looping over each side
      for (const auto & side : primalElem->side_index_range())
      {
        if (primalElem->neighbor_ptr(side) == nullptr) // if its a boundary side..
        {
          std::unique_ptr<const Elem> side_elem = primalElem->build_side_ptr(side);

          std::vector<const Node *> side_nodes;
          for (unsigned int i = 0; i < side_elem->n_nodes(); ++i)
          {
            const Node * node = side_elem->node_ptr(i);
            side_nodes.push_back(node);
          }

          // find midpoint.
          Point midPoint;
          for (const auto & node : side_nodes)
          {
            midPoint += *node;
          }
          midPoint /= side_nodes.size();

          for (auto * node : side_nodes)
            node_to_boundary_midpoints[node->id()].push_back(midPoint);

          // Recording boundary edge connectivity
          if (side_nodes.size() == 2)
          {
            const auto id0 = side_nodes[0]->id();
            const auto id1 = side_nodes[1]->id();

            node_to_boundary_neighbors[id0].push_back(id1);
            node_to_boundary_neighbors[id1].push_back(id0);
          }
        }
      }
    }
    // boundaryMidPoints now contains all of the node pointers that are boundary nodes that might
    // need to be added to a dual mesh.

    // Helper for determining if a node is a vertex of the mesh
    auto isBoundaryVertex = [&](dof_id_type node_id) -> bool
    {
      auto it = node_to_boundary_neighbors.find(node_id);
      if (it == node_to_boundary_neighbors.end())
        return false;

      auto neighbors = it->second;

      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

      if (neighbors.size() != 2) // boundary vertexes have 2 adjacent sides
        return false;
      // calculating angle between sides to exclude nodes along colinear sides
      const Point & p = *mesh->node_ptr(node_id);
      Point v0 = *mesh->node_ptr(neighbors[0]) - p;
      Point v1 = *mesh->node_ptr(neighbors[1]) - p;

      const Real n0 = v0.norm();
      const Real n1 = v1.norm();

      if (n0 == 0.0 || n1 == 0.0)
        return false;

      Real c = (v0 * v1) / (n0 * n1);
      c = std::max(Real(-1), std::min(Real(1), c));

      const Real angle_deg = std::acos(c) * 180.0 / libMesh::pi;

      return std::abs(angle_deg - 180.0) >
             _boundary_node_angular_tol; // angular tolerance; see params
    };

    // Helper for calculating area --> used to determine concavity
    auto signedArea = [](const std::vector<Node *> & nodes) -> Real
    {
      Real area = 0.0;

      for (unsigned int i = 0; i < nodes.size(); ++i)
      {
        const Point & p = *nodes[i];
        const Point & q = *nodes[(i + 1) % nodes.size()];

        area += p(0) * q(1) - q(0) * p(1);
      }

      return 0.5 * area;
    };

    auto cross2D = [](const Point & a, const Point & b, const Point & c) -> Real
    { return (b(0) - a(0)) * (c(1) - a(1)) - (b(1) - a(1)) * (c(0) - a(0)); };

    // Helper to determine concavity; used to split concave boundaries into TRIs
    auto isConcavePolygon = [&](const std::vector<Node *> & nodes) -> bool
    {
      if (nodes.size() < 4)
        return false;

      const Real orientation = signedArea(nodes);

      for (unsigned int i = 0; i < nodes.size(); ++i)
      {
        const Point & prev = *nodes[(i + nodes.size() - 1) % nodes.size()];
        const Point & curr = *nodes[i];
        const Point & next = *nodes[(i + 1) % nodes.size()];

        const Real cross = cross2D(prev, curr, next);

        if (orientation > 0.0 && cross < -TOLERANCE)
          return true;

        if (orientation < 0.0 && cross > TOLERANCE)
          return true;
      }

      return false;
    };

    // Helper for creating triangular elements
    auto addTriangle = [&](Node * a, Node * b, Node * c)
    {
      auto tri = std::make_unique<libMesh::C0Polygon>(3);
      tri->set_node(0, a);
      tri->set_node(1, b);
      tri->set_node(2, c);
      dualMesh->add_elem(std::move(tri));
    };

    // Helper to determine if node is a corner
    auto isCornerNode = [&](const Node * node) -> bool
    { return corner_node_ids.find(node->id()) != corner_node_ids.end(); };

    // Helper to determine if node is a midpoint
    auto isMidpointNode = [&](const Node * node) -> bool
    { return midpoint_node_ids.find(node->id()) != midpoint_node_ids.end(); };

    // Helper to add finalized polygons (convex or triangular) to the mesh
    auto addConvexOrTriangulatedPolygon = [&](std::vector<Node *> nodes)
    {
      if (nodes.size() < 3)
        return;

      if (!isConcavePolygon(nodes))
      {
        auto elem = std::make_unique<libMesh::C0Polygon>(nodes.size());

        for (unsigned int i = 0; i < nodes.size(); ++i)
          elem->set_node(i, nodes[i]);

        dualMesh->add_elem(std::move(elem));
        return;
      }

      Node * cornerNode = nullptr;

      for (auto * node : nodes)
      {
        if (isCornerNode(node))
        {
          cornerNode = node;
          break;
        }
      }

      if (!cornerNode)
        mooseError("Concave polygon was detected, but no corner node was found.");

      std::vector<std::pair<Node *, Real>> sorted_nodes;

      for (auto * node : nodes)
      {
        if (node == cornerNode)
          continue;

        const Real phi = std::atan2((*node)(1) - (*cornerNode)(1), (*node)(0) - (*cornerNode)(0));

        sorted_nodes.push_back({node, phi});
      }

      std::sort(sorted_nodes.begin(),
                sorted_nodes.end(),
                [](const auto & a, const auto & b) { return a.second < b.second; });
      std::vector<Node *> fan_nodes;

      for (unsigned int i = 0; i < sorted_nodes.size(); ++i)
      {
        if (!isMidpointNode(sorted_nodes[i].first))
          continue;

        const unsigned int next_i = (i + 1) % sorted_nodes.size();
        const unsigned int prev_i = (i + sorted_nodes.size() - 1) % sorted_nodes.size();

        // Prefer the direction where the first step from the boundary midpoint
        // goes to an interior/circumcenter node, not another boundary midpoint.
        if (!isMidpointNode(sorted_nodes[next_i].first))
        {
          for (unsigned int k = 0; k < sorted_nodes.size(); ++k)
            fan_nodes.push_back(sorted_nodes[(i + k) % sorted_nodes.size()].first);

          break;
        }

        if (!isMidpointNode(sorted_nodes[prev_i].first))
        {
          for (unsigned int k = 0; k < sorted_nodes.size(); ++k)
            fan_nodes.push_back(
                sorted_nodes[(i + sorted_nodes.size() - k) % sorted_nodes.size()].first);

          break;
        }
      }

      if (fan_nodes.empty())
        mooseError("Could not find a boundary midpoint to start concave fan triangulation.");

      for (unsigned int i = 0; i + 1 < fan_nodes.size(); ++i)
        addTriangle(cornerNode, fan_nodes[i], fan_nodes[i + 1]);
    };

    // loop over all primal nodes / dual elements
    for (const auto & [primalNodeID, primalElemIDs] : _node_to_elem_map)
    {
      const bool is_boundary_node =
          node_to_boundary_midpoints.find(primalNodeID) != node_to_boundary_midpoints.end();
      std::vector<std::pair<Node *, Real>> dualNodesAndPhis;
      if (!is_boundary_node) // For interior polygons we needn't worry about the boundary
      {
        Node * primalNode = mesh->node_ptr(primalNodeID);
        for (unsigned int j = 0; j < primalElemIDs.size(); ++j)
        {
          const dof_id_type dualNodeOnPElem_id =
              primalElemIDs[j]; // Grab the dual nodes' IDs on each primal element

          auto const dualNodeOnPElem = dualMesh->node_ptr(
              dualNodeOnPElem_id); // Grab the dual nodes corresponding to these IDs

          Real dualNodeX = (dualNodeOnPElem->operator()(0)) - primalNode->operator()(0);
          Real dualNodeY = (dualNodeOnPElem->operator()(1)) - primalNode->operator()(1);
          Real nodePhi = atan2(dualNodeY, dualNodeX);

          dualNodesAndPhis.push_back({dualNodeOnPElem, nodePhi});

          std::sort(dualNodesAndPhis.begin(),
                    dualNodesAndPhis.end(),
                    [](const auto & a, const auto & b) { return a.second < b.second; });
        }
      }
      else // boundary dual elements
      {

        Node * primalNode = mesh->node_ptr(primalNodeID);

        if (isBoundaryVertex(primalNodeID))
        {
          Node * cornerNode = dualMesh->add_point(*mesh->node_ptr(primalNodeID));
          corner_node_ids.insert(cornerNode->id());

          dualNodesAndPhis.push_back({cornerNode, 0.0});
        }

        // Add circumcenter nodes from adjacent primal elements
        for (const auto elem_id : primalElemIDs)
        {
          Node * dualNode = dualMesh->node_ptr(elem_id);
          dualNodesAndPhis.push_back({dualNode, 0.0});
        }

        // Add boundary midpoint nodes
        for (const auto & midpoint : node_to_boundary_midpoints[primalNodeID])
        {
          Node * midpointNode = dualMesh->add_point(midpoint);
          midpoint_node_ids.insert(midpointNode->id());
          dualNodesAndPhis.push_back({midpointNode, 0.0});
        }

        // Recompute angles around the geometric center of boundary elements

        Point circumcenter_avg;

        for (const auto elem_id : primalElemIDs)
          circumcenter_avg += *dualMesh->node_ptr(elem_id);

        circumcenter_avg /= primalElemIDs.size();

        Point sort_center = 0.5 * ((*primalNode) + circumcenter_avg);

        for (auto & [node, phi] : dualNodesAndPhis)
        {
          phi = std::atan2((*node)(1) - sort_center(1), (*node)(0) - sort_center(0));
        }
      }
      std::sort(dualNodesAndPhis.begin(),
                dualNodesAndPhis.end(),
                [](const auto & a, const auto & b) { return a.second < b.second; });

      std::vector<Node *> ordered_nodes;
      ordered_nodes.reserve(dualNodesAndPhis.size());

      for (const auto & node_phi : dualNodesAndPhis)
        ordered_nodes.push_back(node_phi.first);

      addConvexOrTriangulatedPolygon(ordered_nodes);
    }
  }
  // 3D Not implemented
  else if (mesh->mesh_dimension() == 3)
  {
    mooseError("3D mesh support not implemented yet");
    auto makecircumcenterFace = [&](std::vector<Node *> face_nodes)
    {
      std::shared_ptr<libMesh::C0Polygon> face =
          std::make_shared<libMesh::C0Polygon>(face_nodes.size());

      mooseAssert(face->n_nodes() == face_nodes.size(),
                  "C0Polygon face was not initialized with the requested number of nodes.");

      for (unsigned int i = 0; i < face_nodes.size(); ++i)
        face->set_node(i, face_nodes[i]);

      return std::static_pointer_cast<libMesh::Polygon>(face);
    };

    auto ordercircumcentersAroundEdge =
        [&](std::vector<Node *> nodes, const Point & edge_p0, const Point & edge_p1)
    {
      const Point axis = edge_p1 - edge_p0;

      Point center;
      for (const auto * node : nodes)
        center += *node;
      center /= nodes.size();

      Point ref = *nodes[0] - center;
      Point normal = axis.cross(ref);

      std::sort(nodes.begin(),
                nodes.end(),
                [&](const Node * a, const Node * b)
                {
                  Point va = *a - center;
                  Point vb = *b - center;

                  Real aa = std::atan2((axis.cross(ref)) * va, ref * va);
                  Real bb = std::atan2((axis.cross(ref)) * vb, ref * vb);

                  return aa < bb;
                });

      return nodes;
    };

    for (const auto & [primalNodeID, primalElemIDs] : _node_to_elem_map)
    {
      // Skip boundary primal nodes for now
      bool is_boundary_node = false;

      for (const auto elem_id : primalElemIDs)
      {
        const Elem * elem = mesh->elem_ptr(elem_id);

        for (const auto side : elem->side_index_range())
        {
          if (elem->neighbor_ptr(side) == nullptr)
          {
            std::unique_ptr<const Elem> side_elem = elem->build_side_ptr(side);

            for (unsigned int sn = 0; sn < side_elem->n_nodes(); ++sn)
            {
              if (side_elem->node_id(sn) == primalNodeID)
              {
                is_boundary_node = true;
                break;
              }
            }
          }

          if (is_boundary_node)
            break;
        }

        if (is_boundary_node)
          break;
      }

      if (is_boundary_node)
        continue;

      std::vector<std::shared_ptr<libMesh::Polygon>> poly_faces;

      // One dual face per primal edge attached to this primal node.
      std::unordered_map<dof_id_type, std::vector<dof_id_type>> edge_to_elem_ids;

      for (const auto elem_id : primalElemIDs)
      {
        const Elem * elem = mesh->elem_ptr(elem_id);

        for (unsigned int n = 0; n < elem->n_nodes(); ++n)
        {
          const dof_id_type other_node_id = elem->node_id(n);

          if (other_node_id == primalNodeID)
            continue;

          // If both nodes appear in this elem, this elem contributes
          // to the dual face associated with edge primalNodeID-other_node_id.
          edge_to_elem_ids[other_node_id].push_back(elem_id);
        }
      }

      for (const auto & [other_node_id, elem_ids] : edge_to_elem_ids)
      {
        if (elem_ids.size() < 3)
          continue;

        std::vector<Node *> face_nodes;

        for (const auto elem_id : elem_ids)
          face_nodes.push_back(dualMesh->node_ptr(elem_id));

        const Point & p0 = *mesh->node_ptr(primalNodeID);
        const Point & p1 = *mesh->node_ptr(other_node_id);

        face_nodes = ordercircumcentersAroundEdge(face_nodes, p0, p1);

        poly_faces.push_back(makecircumcenterFace(face_nodes));
      }

      if (poly_faces.size() < 4)
        continue;

      std::unique_ptr<Node> mid_elem_node;

      auto poly = std::make_unique<libMesh::C0Polyhedron>(poly_faces, mid_elem_node);

      dualMesh->add_elem(std::move(poly));
    }
  }
  else
  {
    _console << "ERROR: Tried to generate dual of a " << mesh->mesh_dimension()
             << "-dimensional mesh." << std::endl;
    mooseError("DualMeshGenerator only supports 2D or 3D meshes.");
  }

  _console << "Dual mesh nodes: " << dualMesh->n_nodes() << std::endl;
  _console << "Dual mesh elems: " << dualMesh->n_elem() << std::endl;

  dualMesh->unset_is_prepared();
  return dynamic_pointer_cast<MeshBase>(dualMesh);
}