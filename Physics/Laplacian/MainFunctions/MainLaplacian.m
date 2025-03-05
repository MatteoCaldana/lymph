%> @file  MainLaplacian.m
%> @author The Lymph Team
%> @date 26 July 2024
%> @brief Solution of the Poisson problem with PolydG
%>
%==========================================================================
%> @section classMainLaplacian Class description
%==========================================================================
%> @brief            Solution of the Poisson problem with PolydG
%
%> @param Data       Struct with problem's data
%> @param Setup      Simulation's setup
%>
%> @retval Error     Struct for convergence test
%>
%==========================================================================

function [Error] = MainLaplacian(Data,Setup)

fprintf('\nSolution of the Poisson problem \n');
fprintf('... here display some info about the simulation ... \n');
fprintf('\n------------------------------------------------------------------\n')

%% Mesh setup
[mesh, femregion, Data.h] = MeshFemregionSetup(Setup, Data, {Data.TagElLap}, {'L'});
timings = struct('deg', Data.degree, 'int_strategy', Data.quadrature , 'mesh_kind', Data.mesh_kind, 'nel', Data.nel, 'solve', 0, 'compute_error', 0, 'rhs', 0, 'assembly', 0);

%% Matrix Assembly
fprintf('\nMatrices computation ... \n');
tic
switch Data.quadrature
    case "QF"
        [Matrices] = MatrixLaplacianQF(Data, mesh.neighbor, femregion);
    case "ST"
        [Matrices] = MatrixLaplacianST(Data, mesh.neighbor, femregion);
end
assembly = toc;
timings.assembly = assembly;
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Right-hand side assembly

fprintf('\nComputing RHS ... \n');
tic
[F] = ForcingLaplacian(Data, mesh.neighbor, femregion);
rhs = toc;
timings.rhs = rhs;
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Solving the linear system
fprintf('\nSolving linear system ... ');
tic
U = Matrices.A \ F;
solve = toc;
timings.solve = solve;
fprintf('\nDone\n')
fprintf('\n------------------------------------------------------------------\n')

%% Postprocess solution
if Setup.isSaveCSV || Setup.isSaveVTK || Setup.isSaveSolution || Setup.isPlotSolution
    PostProcessSolution(Setup, Data, mesh, femregion, 0, U);
end

%% Compute Errors
if Setup.isError
    fprintf('\nComputing Errors ... \n');
    tic
    [Error] = ComputeErrors(Data, femregion, Matrices, U);
    compute_error = toc;
    timings.compute_error = compute_error;
    fprintf('Done\n')
    fprintf('\n------------------------------------------------------------------\n')
else
    Error = [];
end

fprintf('\nBye \n');

fid = fopen(sprintf('timing_d%d_%s_%s_n%d.json',Data.degree, Data.quadrature, Data.mesh_kind, Data.nel),'w');
fprintf(fid,'%s',jsonencode(timings));
fclose(fid);



