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
  MooseEnum location("INSIDE OUTSIDE", "INSIDE");

  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addClassDescription("Takes a 2D mesh as input and returns a Voronoi dual mesh, i.e.,"
                             "changes each input mode into an element and each input element "
                             "into a node located at its centroid.");
  params.addRequiredParam<RealVectorValue>(
      "bottom_left", "The bottom left point (in x,y,z with spaces in-between).");
  params.addRequiredParam<RealVectorValue>(
      "top_right", "The bottom left point (in x,y,z with spaces in-between).");
  params.addRequiredParam<subdomain_id_type>(
      "block_id", "Subdomain id to set for inside/outside the bounding box");
  params.addParam<SubdomainName>(
      "block_name", "Subdomain name to set for inside/outside the bounding box (optional)");
  params.addParam<MooseEnum>(
      "location", location, "Control of where the subdomain id is to be set");
  params.addParam<std::vector<SubdomainName>>(
      "restricted_subdomains",
      "Only reset subdomain ID for given subdomains within the bounding box");

  params.addParam<std::string>("integer_name",
                               "Element integer to be assigned (default to subdomain ID)");
  return params;
}

DualMeshGenerator::DualMeshGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _location(parameters.get<MooseEnum>("location")),
    _block_id(parameters.get<subdomain_id_type>("block_id")),
    _has_restriction(isParamValid("restricted_subdomains")),
    _bounding_box(MooseUtils::buildBoundingBox(parameters.get<RealVectorValue>("bottom_left"),
                                               parameters.get<RealVectorValue>("top_right")))
{
}

std::unique_ptr<MeshBase>
DualMeshGenerator::generate()
{
  auto mesh = std::move(_input);
  std::vector<libMesh::Point> centroids; // vector of all dual nodes
  centroids.reserve(mesh->n_elem());

  std::unordered_map<dof_id_type, std::vector<dof_id_type>> _node_to_elem_map;

  unsigned int i;

  //////////Adding nodes//////
  for (const auto & in_elem : mesh->element_ptr_range())
  {

    _console << "Found centroid: " << in_elem->true_centroid() << std::endl;

    Point centroid = in_elem->true_centroid();
    centroids.push_back(centroid);
    for (unsigned int n = 0; n < in_elem->n_nodes(); n++)
    {
      _node_to_elem_map[in_elem->node_id(n)].push_back(in_elem->id());
    }
  }
  _console << "Successfully created node-element mapping" << std::endl;

  auto dualMesh = buildReplicatedMesh(mesh->mesh_dimension());

  for (i = 0; i < centroids.size(); ++i)
  {
    dualMesh->add_point(centroids[i]);
  }

  std::unordered_map<dof_id_type, std::vector<dof_id_type>> _elem_to_node_map;
  for (const auto & [node_id, elements] : _node_to_elem_map)
  {
    for (dof_id_type elem_id : elements)
    {
      _elem_to_node_map[elem_id].push_back(node_id);
    }
  } //_elem_to_node_map now has element IDs in the first entry and maps to the nodes that make
    // up that element

  _console << "Mesh populated with dual nodes" << std::endl;

  // loop over all primal nodes / dual elements
  for (const auto & [primalNodeID, primalElemIDs] : _node_to_elem_map)
  {
    _console << "Number of nodes for dual element: " << primalElemIDs.size() << std::endl;

    if (primalElemIDs.size() >= 3)
    {

      _console << "Loading interor polygon!" << std::endl;

      // Define a dual element located at each primal node
      std::unique_ptr<Elem> dualElem = std::make_unique<libMesh::C0Polygon>(primalElemIDs.size());

      // Now loop over the # of nodes on each dual element
      std::vector<std::pair<Node *, Real>> dualNodesAndPhis;
      auto primalNode = mesh->node_ptr(primalNodeID);

      for (unsigned int j = 0; j < primalElemIDs.size();
           ++j) // n.second.size is number of dual nodes to the dual element
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
      for (unsigned int k = 0; k < dualNodesAndPhis.size(); ++k)
      {
        dualElem->set_node(k,
                           dualNodesAndPhis[k].first); // assign these nodes to the the dual element
      }

      // add the element to the mesh, now that it's assigned nodes
      dualMesh->add_elem(std::move(dualElem));
    }
    else
    {
      // Here are the elements with less than 3 dual nodes, that need to have midpoints added to be
      // polygons.
      // _____GETTING BOUNDARY NODES______ //
      // loop over all primal elements
      for (const auto & [elemID, nodeIDs] : _elem_to_node_map)
      {
        Elem * primalElem = mesh->elem_ptr(elemID);
        std::vector<Node *> extDualNodesOnPrimalElem;
        for (const auto & side : primalElem->side_index_range())
        {
          if (primalElem->neighbor_ptr(side) == nullptr)
          {
            std::vector<unsigned int> extPrimalNodeIDs = primalElem->nodes_on_side(side);
            _console << "found bordering side with nodes: " << std::endl;
            for (const auto & nodeID : extPrimalNodeIDs)
            {
              Node * node = mesh->node_ptr(nodeID);
              node->print_info();
            }

            // find midpoint.
            Point midPoint;
            for (const auto & nodeID : extPrimalNodeIDs)
            {
              Node * node = mesh->node_ptr(nodeID);
              midPoint += *node;
            }
            midPoint /= extPrimalNodeIDs.size();

            Node * midpointNode = dualMesh->add_point(midPoint);
            extDualNodesOnPrimalElem.push_back(midpointNode);
            //_console << "found mipoint with info" << std::endl;
            // midpointNode->print_info();
          }
        }

        // We have a vector of points that should be dual nodes. Add them to a dual element!
        std::unique_ptr<Elem> extDualElem = std::make_unique<libMesh::C0Polygon>(
            primalElemIDs.size() + extDualNodesOnPrimalElem.size());
        // Lets compile the interior nodes and the boundary nodes into a single vector for sorting
        std::vector<std::pair<Node *, Real>> nodesAndPhis;
        auto primalNode = mesh->node_ptr(primalNodeID);

        // Add the interior nodes to the vector
        for (unsigned int k = 0; k < primalElemIDs.size(); ++k)
        {
          Node * intNode = dualMesh->node_ptr(primalElemIDs[k]);
          nodesAndPhis.push_back({intNode, 0}); // Phi placeholder
        }

        // Adding boundary nodes to the vector
        for (unsigned int k = 0; k < extDualNodesOnElem.size(); ++k)
        {
          nodesAndPhis.push_back({extDualNodesOnElem[k], 0});
          // extDualElem->set_node((k + primalElemIDs.size()), extDualNodesOnElem[k]);
        }

        for (unsigned int j = 0; j < nodesAndPhis.size(); ++j)
        {
          Node * node = nodesAndPhis[j].first;

          Real nodeX = (node->operator()(0)) - primalNode->operator()(0);
          Real nodeY = (node->operator()(1)) - primalNode->operator()(1);
          Real nodePhi = atan2(nodeY, nodeX);

          nodesAndPhis[j].second = nodePhi;
        }
        std::sort(nodesAndPhis.begin(),
                  nodesAndPhis.end(),
                  [](const auto & a, const auto & b) { return a.second < b.second; });
        // Now they are sorted by their Phi. Now set them to the element in their order

        for (unsigned int k = 0; k < nodesAndPhis.size(); ++k)
        {
          _console << "Adding point with phi: " << nodesAndPhis[k].second << std::endl;
          extDualElem->set_node(k, nodesAndPhis[k].first);
        }
        //_console << "Adding element to mesh with info: " << std::endl;
        // extDualElem->print_info();

        dualMesh->add_elem(std::move(extDualElem));
        _console << "Next Element!" << std::endl;
      }
    }
  }

  //_console << "Printing info" << std::endl;
  dualMesh->print_info();
  //_console << "Finished printing info" << std::endl;

  dualMesh->unset_is_prepared();
  return dynamic_pointer_cast<MeshBase>(dualMesh);
}
