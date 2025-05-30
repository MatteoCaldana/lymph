# Report changes as follows:
# - **New/Changed/Improved**: DESCRIPTION. (Name Surname YYYY/MM/DD)

## Changes between v1.1.0 and v1.2.0
- **Improved**:`Quadrature`has now (l+1)^2 nodes for integrals over triangles. (Ilario Mazzieri 2024/10/09)
- **Improved**:`Quadrature`files documentation improved. (Mattia Corti 2024/10/09)
- **Changed**:`GetSolutionQuadPoints`->`GetSolutionPostProcessing`: post-processing nodes are independent of quadrature nodes to improve performance. (Mattia Corti 2024/10/09)
- **Improved**: BC imposition in RHS assembly is now more efficient. (Stefano Bonetti 2024/10/09)
- **Improved**: Bugfix in the computation of harmonic coeff for penalty. (Mattia Corti 2024/10/07)
- **New**: Mesh generation for a rectangular bidomain with horizontal interface. (Ilario Mazzieri 2024/10/04)

## Changes between v1.0.0 and v1.1.0
- **Updated**: the PoroElastoAcoustics physics is now updated with the new Core and the Quadrature-Free is implemented. (Mattia Corti  2024/08/13)
- **Improved**: fixed doxygen documentation for all physics. (Ilario Mazzieri 2024/30/07)
- **New**: added documentation for PoroElastoAcoustics. (Ilario Mazzieri 2024/30/07)
- **New**: added a new physcs PoroElastoAcoustics to carry out wave propagation simulation in coupled poroelastic-elastic-acoustic media. (Ilario Mazzieri  2024/30/07)
- **New**: Implemented the quadrature-free assembly for all physics, avoiding sparse allocation. (Mattia Corti 2024/07/29)
- **Changed**: `GetSolutionQuadPoints`: Changed the output of the functions GetSolutionQuadPoints for all the physics, to include the information necessary to the new postprocessing functions. Moreover, the degree of the quadrature is now set from the data file and it is reconstructed also on the edges of each mesh element. (Mattia Corti 2024/07/29)
- **New**: `PostProcessSolution`: Constructed a generic postprocessing function for all physics that manages all the solution plots and saving procedures. Updated time advancement, solution plotting and data structures accordingly. (Mattia Corti 2024/07/29)
- **Improved**: `CreatePolygonalVTK`: The writing of VTK file to save the mesh has been improved to take advantage of fprintf optimization for vectors. (Mattia Corti 2024/07/29)
- **New**: `MeshFemregionSetup`: Constructed a generic initialization procedure for all physics in which the mesh is read and then the femregion structure is created. Updated Error structure and initial conditions (new function `GetInitialConditions`) accordingly. (Mattia Corti 2024/07/29)
- **Changed**: Minor changes in the AssembleNeighEl to improve the scalability. However, the new assembly procedures do not require this function any more. This will be deleted. (Mattia Corti 2024/07/29)
- **Changed**: Updated the main functions to take into account all the last changes of assembly, initialization and postprocessing procedures. (Mattia Corti 2024/07/29)
- **Improved**: Assembly with subtriangulation of all physics now avoid sparse allocation, thus speeding up the computations and the code scalability. (Mattia Corti 2024/07/29)
- **New**: copyright and license information, indicating also the PolyMesher dependency. (Ivan Fumagalli 2024/07/22)
- **New**: New functions for quadrature-free assembly, to evaluate, scale and integrate the monomial expansion of polynomial integrands on polygons and the corresponding Legendre basis functions. (Mattia Corti 2023/11/13)

## Changes between v0.0.0 and v1.0.0
- **New**: added the function ClockWiseElements that control for each element the clockwise storing of the vertices. (Mattia Corti 2023/05/16)
- **New**: added the function MetisMextoRegion that converts the agglomerated meshes in metismex format to lymph one. (Mattia Corti 2023/05/16)
- **New**: added the function CreatePolygonalVTK that allows the saving of the polygonal mesh in VTK format. (Mattia Corti 2023/05/05)
- **Improved**: `Evalshape2D`: The evalshape functions now control the number of outputs and is able to compute only the bases functions of phiq if gradients are not needed. The call needs to be done as phiq = Evalshape2D(...) without [] and without. (Mattia Corti 2023/05/05)
- **Changed**: `Quadrature`: the quadrature function has now a transposed output, which allows to neglect the dx' and ds' in assembly, now dx and ds. (Mattia Corti 2023/05/05)
