[Mesh]
    [mySphere]
        type = SphereMeshGenerator
        nr = 1
        radius = 3
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