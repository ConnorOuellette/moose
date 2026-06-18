[Mesh]
    [mySphere]
        type = SphereMeshGenerator
        nr = 1
        radius = 1
    []

    [SdmPerElemGen]
        type = SubdomainPerElementGenerator
        input = mySphere
    []

    [convert]
        type = ElementsToSimplicesConverter
        input = 'SdmPerElemGen'
    []

    [check]
        type = MeshDiagnosticsGenerator
        input = convert
        examine_element_overlap = WARNING
        examine_non_matching_edges = WARNING
        examine_element_volumes = WARNING
        minimum_element_volumes = 0
    []
[]