[Mesh]
    [myCartMG]
        type = CartesianMeshGenerator
        dim = 3
        dx = "3 3"
        dy = "4 4"
        dz = "5 5"

        subdomain_id = '1 2
                        3 4

                        5 6
                        7 8'
    []

    [cut_1]
        type = BlockDeletionGenerator
        input = myCartMG
        block = '5 6'
    []

    [myDualGen]
        type = DualMeshGenerator
        input = cut_1
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