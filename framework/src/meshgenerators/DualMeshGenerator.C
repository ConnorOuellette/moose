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
    _console << "_______INPUT MESH INFO_______" << std::endl;
    in_elem->print_info();

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

  _console << "Mesh populated with dual nodes" << std::endl;

  // loop over all primal nodes / dual elements
  for (const auto & n : _node_to_elem_map)
  {

    // Define a dual element located at each primal node
    std::unique_ptr<Elem> elem = std::make_unique<libMesh::C0Polygon>(n.second.size());
    // Now loop over the # of nodes on each dual element

    _console << "Number of nodes for dual element: " << n.second.size() << std::endl;
    if (n.second.size() < 3)
    {
      _console << "Skipping unbounded polygon with " << n.second.size() << " node(s)" << std::endl;
      continue;
    }
    else
      _console << "Loading interor polygon!" << std::endl;

    for (unsigned int j = 0; j < n.second.size(); ++j)
    {
      const dof_id_type dualNodeOnPElem_id =
          n.second[j]; // Grab the dual nodes' IDs on each primal element

      auto dualNodeOnPElem =
          dualMesh->node_ptr(dualNodeOnPElem_id); // Grab the dual nodes corresponding to these IDs

      // assign these nodes to the polygon element
      elem->set_node(j, dualNodeOnPElem); // assign these nodes to the the dual element
    }
    // add the element to the mesh, now that it's assigned nodes
    dualMesh->add_elem(std::move(elem));
  }

  _console << "Printing info" << std::endl;
  dualMesh->print_info();
  _console << "Finished printing info" << std::endl;

  dualMesh->unset_is_prepared();
  return dynamic_pointer_cast<MeshBase>(dualMesh);
}
