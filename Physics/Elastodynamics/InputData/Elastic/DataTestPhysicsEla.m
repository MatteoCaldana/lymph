% Domain: 
% Meshfile: MeshEmilia
 
Data.name = 'Emilia';

Data.TagElEla   = [1 2 3 4 5 6 7]; % Element tag
Data.TagBcEla   = [8 9 10 11];     % Boundary tag
Data.LabBcEla   = 'DDND';          % Dirichlet/Neumann/Abso/

%% Geometrical properties 

Data.domain          = [0  38463.258537  0 10000];
Data.N               = 0;  %number of elements
Data.MeshFromFile    = true;
Data.FolderName      = 'InputMesh';
Data.VTKMeshFileName = 'MeshEmilia.vtk';
Data.meshfileseq     = ["MeshEmilia.mat"]; %filename for mesh 

%% Discretization properties                            

% Time integration
Data.t0      = 0;
Data.T       = 4;
Data.dt      = 0.001;
Data.timeint = 'newmark';
Data.BetaNM  = 0.25;
Data.GammaNM = 0.5;

% Space discretization
Data.degree        = 5;  % Polynomial degree
Data.penalty_coeff = 10; % Penalty coefficient

%% Quadrature settings

Data.quadrature = "QF";       % Quadrature type: ST/QF

%% Visualization settings
Data.PlotExact         = false;
Data.PlotGridSol       = false;
Data.VisualizationStep = 5;
Data.PlotIniCond       = true;
Data.NPtsVisualization = 5;

%% properties elastic material
Data.rho_el    = {@(x,y) 1800 + 0.*x.*y; @(x,y) 1800 + 0.*x.*y; ...
                  @(x,y) 2050 + 0.*x.*y; @(x,y) 2050 + 0.*x.*y; ...
                  @(x,y) 2050 + 0.*x.*y; @(x,y) 2400 + 0.*x.*y; ...
                  @(x,y) 2450 + 0.*x.*y};

Data.vs_el     = {@(x,y) 294 + 0.*x.*y; @(x,y) 450 + 0.*x.*y; ...
                  @(x,y) 600 + 0.*x.*y; @(x,y) 600 + 0.*x.*y; ...
                  @(x,y) 600 + 0.*x.*y; @(x,y) 1515 + 0.*x.*y; ...
                  @(x,y) 1600 + 0.*x.*y; };

Data.vp_el     = {@(x,y) 1321 + 0.*x.*y; @(x,y) 2024 + 0.*x.*y; ...
                  @(x,y) 1920 + 0.*x.*y; @(x,y) 1920 + 0.*x.*y; ...
                  @(x,y) 1920 + 0.*x.*y; @(x,y) 3030 + 0.*x.*y; ...
                  @(x,y) 3200 + 0.*x.*y; };

Data.zeta      = {@(x,y) 0 + 0.*x.*y; @(x,y) 0 + 0.*x.*y; ...
                  @(x,y) 0 + 0.*x.*y; @(x,y) 0 + 0.*x.*y; ...
                  @(x,y) 0 + 0.*x.*y; @(x,y) 0 + 0.*x.*y; ...
                  @(x,y) 0 + 0.*x.*y; };

Data.mu_el     = {@(x,y) 155584800  + 0.*x.*y; @(x,y) 364500000  + 0.*x.*y;...
                  @(x,y) 738000000  + 0.*x.*y; @(x,y) 738000000  + 0.*x.*y;...
                  @(x,y) 738000000  + 0.*x.*y; @(x,y) 5.5085e+09 + 0.*x.*y;...
                  @(x,y) 6.2720e+09 + 0.*x.*y}; % Data.vs_el^2 * Data.rho_el;

Data.lam_el    = {@(x,y) 0.2830e+10 + 0.*x.*y; @(x,y) 0.6645e+10 + 0.*x.*y;...
                  @(x,y) 0.6081e+10 + 0.*x.*y; @(x,y) 0.6081e+10 + 0.*x.*y;...
                  @(x,y) 0.6081e+10 + 0.*x.*y; @(x,y) 1.1017e+10 + 0.*x.*y;...
                  @(x,y) 1.2544e+10 + 0.*x.*y}; % Data.vp_el^2 * Data.rho_el - 2*Data.mu_el;

% forcing term elastic media
Data.source_f = {@(x,y,t)  0.*t.*x.*y; ...
                 @(x,y,t)  0.*t.*x.*y};

Data.source_g = {@(x,y,t)  0.*t.*x.*y; ...
                 @(x,y,t)  0.*t.*x.*y};

% Moment source tensor
Data.sourceMxx_el   = {@(x,y,t) (1-2*pi^2*4*(t-0.5).^2).*exp(-pi^2*4*(t-0.5).^2).*((x-19432).^2 + (y-7800).^2 <= 100^2)}; % Forcing term Mxx
Data.sourceMyy_el   = {@(x,y,t) (1-2*pi^2*4*(t-0.5).^2).*exp(-pi^2*4*(t-0.5).^2).*((x-19432).^2 + (y-7800).^2 <= 100^2)}; % Forcing term Myy
Data.sourceMxy_el   = {@(x,y,t) 0*t.*x.*y}; % Forcing term Mxy
Data.sourceMyx_el   = {@(x,y,t) 0*t.*x.*y}; % Forcing term Myx

% Dirichlet BC
Data.DirBCEla  = {@(x,y,t)  0*t.*x.*y; ...
                  @(x,y,t)  0*t.*x.*y};

% exact solution --> used to compute the initial conditions
Data.ue_ex = {@(x,y,t) 0.*t.*x.*y; ...
              @(x,y,t) 0.*t.*x.*y};

Data.due_dt_ex = {@(x,y,t) 0.*t.*x.*y; ...
                  @(x,y,t) 0.*t.*x.*y};

% exact gradient --> used for the error analysis
% du1/dx, du1/dy, du2/dx, du2/dy
Data.grad_ue_ex =  {@(x,y,t) 0.*t.*x.*y; ...
                    @(x,y,t) 0.*t.*x.*y; ...
                    @(x,y,t) 0.*t.*x.*y; ...
                    @(x,y,t) 0.*t.*x.*y};
