[Mesh]
    [myCCMG]
        type = ConcentricCircleMeshGenerator
        num_sectors = 2
        radii = '0.5 1.0'
        rings = '1 1'
        has_outer_square = false
        preserve_volumes = false
        smoothing_max_it = 0
    []

    [extrude]
        type = MeshExtruderGenerator
        extrusion_vector = '0 0 1'
        input = myCCMG
    []

    [myDualGen]
        type = DualMeshGenerator
        input = extrude
    []

    [check]
        type = MeshDiagnosticsGenerator
        input = myDualGen
        examine_element_overlap = WARNING
        examine_non_matching_edges = WARNING
        examine_element_volumes = WARNING
        minimum_element_volumes = 0
    []

    [SdmPerElemGen]
        type = SubdomainPerElementGenerator
        input = check
    []

    [convert]
        type = ElementsToSimplicesConverter
        input = 'SdmPerElemGen'
    []
[]