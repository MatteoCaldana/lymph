%> @file  RunhConvergenceLaplacian.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Convergence analysis for the Poisson problem (mesh refinements)
%>
%==========================================================================
%> @section classRunhConvergenceLaplacian Class description
%==========================================================================
%> @brief          Sequence of run of MainLaplacian.m
%
%> @param ~
%>
%> @retval ~
%>
%==========================================================================

%% Import lymph and paths of folders related to this problem
run("../ImportLymphPaths.m")
MyPhysicsPath = pwd;
addpath(genpath(fullfile(MyPhysicsPath,'Assembly')));
addpath(genpath(fullfile(MyPhysicsPath,'InputData')));
addpath(genpath(fullfile(MyPhysicsPath,'MainFunctions')));
addpath(genpath(fullfile(MyPhysicsPath,'Error')));
addpath(genpath(fullfile(MyPhysicsPath,'PostProcessing')));

%% Simulation - Setup
run("../RunSetup.m")

for deg = 1:5
for int_strategy = ["QF", "ST"]
for mesh_kind = ['T', 'P']

% Input Data - Boundary conditions - Forcing term
DataConvTestLap;

Data.quadrature = int_strategy; 
Data.degree = deg;

Errors.err_L2 = [];
Errors.err_dG = [];
Errors.h = [];

% Mesh Generation
for ii = 1:numel(Data.N)
    Data.mesh_kind = mesh_kind;
    Data.nel = Data.N{ii};
    kind = mesh_kind;
    rng("default")
    Data.VTKMeshFileName = sprintf('Mesh_%s_%d.vtk', kind, Data.N{ii});
    if Data.MeshFromFile
        % Load an existing mesh
        Data.meshfile = fullfile(Data.FolderName,Data.meshfileseq{ii});
    else
        % Create a new mesh
        if kind ~= 'P'
            nn = round(sqrt(Data.N{ii}));
            NN = [nn, nn];
        else
            NN = Data.N{ii};
        end
        [Data.meshfile] = MakeMeshMonodomain(Data,NN,Data.domain,Data.FolderName,Data.meshfileseq{ii},kind,'laplacian');
    end
    
    
    % Main     
    [Error] = MainLaplacian(Data,Setup);
    Errors.err_L2   = [Errors.err_L2, Error.L2];
    Errors.err_dG   = [Errors.err_dG, Error.dG];
    Errors.h        = [Errors.h, Error.h];

    
end

% Plot of the errors
figure
loglog(Errors.h,Errors.h.^Data.degree,'k:','Linewidth',2)
hold on
loglog(Errors.h,Errors.h.^(Data.degree+1),'k--','Linewidth',2)
loglog(Errors.h,Errors.err_L2,'g','Linewidth',2)
loglog(Errors.h,Errors.err_dG,'r','Linewidth',2)
xlabel('h');
conv1 = ['$h^', num2str(Data.degree), '$'];
conv2 = ['$h^', num2str(Data.degree+1), '$'];
legend(conv1, conv2, "Error $L^2$-norm", "Error DG-norm","Interpreter","latex")
grid on
Errors.order_L2 = log(Errors.err_L2(1:end-1)./Errors.err_L2(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));
Errors.order_dG = log(Errors.err_dG(1:end-1)./Errors.err_dG(2:end))./log(Errors.h(1:end-1)./Errors.h(2:end));

close all
end
end
end