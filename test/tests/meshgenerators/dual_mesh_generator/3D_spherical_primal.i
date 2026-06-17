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
[]