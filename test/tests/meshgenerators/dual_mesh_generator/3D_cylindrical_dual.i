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

    [SdmPerElemGen]
        type = SubdomainPerElementGenerator
        input = myDualGen
    []

    [convert]
        type = ElementsToSimplicesConverter
        input = 'SdmPerElemGen'
    []
[]