[Mesh]
    [mySphere]
        type = SphereMeshGenerator
        nr = 1
        radius = 1
    []

    [myDualGen]
        type = DualMeshGenerator
        input = mySphere
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