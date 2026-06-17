[Mesh]
    [mySphere]
        type = SphereMeshGenerator
        nr = 1
        radius = 3
    []

    [SdmPerElemGen]
        type = SubdomainPerElementGenerator
        input = mySphere
    []

    [convert]
        type = ElementsToSimplicesConverter
        input = 'SdmPerElemGen'
    []
[]