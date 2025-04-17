%> @file  MakeMeshMonodomain.m
%> @author Ilario Mazzieri
%> @date 8 March 2023 
%> @brief Construction of a polygonal mesh for a rectangular domain.
%>  
%> Creation of a polygonal mesh for a rectangular domain
%>   Omega = [xmin, xmax] x [ymin, ymax]
%>
%>   IDs for volume and boundary elements
%>            _____4_____
%>           |           |
%>          5|     1     |3   
%>           |___________|    
%>                 2
%>
%==========================================================================
%> @section classMakeMeshMonodomain Class description
%==========================================================================
%> @brief  Construction of the polygonal mesh (uses polymesher functions
%>contained in the Polymesher folder 
%>
%> @param Data                    Struct with problem's data
%> @param N                       Number of mesh elements
%> @param DomainLimits            Domain limits
%> @param FolderName              Directory name for saving
%> @param FileName                File name for saving
%> @param MeshType                String 'C' for cartesian grid, 'P' for
%polygonal grid
%> @param SimType                 String simulation type, used for boundary tag       
%>
%> @retval FileNameOut            File name of the *.mat structure 
%>                                containing mesh info 
%>
%==========================================================================
function [FileNameOut] = MakeMeshMonodomain(Data,N,DomainLimits,FolderName,FileName,MeshType,SimType) 
rng("default")
%% Set directories and names

if ~exist(FolderName,'dir')
    mkdir(char(FolderName));
end

%% RECTANGULAR DOMAIN
global Dati

%xmin xmax ymin ymax
Dati.domain = DomainLimits;
Data.domain = DomainLimits;

if nargin < 7
    SimType = 'laplacian';
    Data.TagBcLap = [2, 3, 4, 5];
    Data.LabBcLap = 'DDDD';
end

if (strcmp(MeshType,'P')==1)
    
    % polygonal mesh
    [region] = PolyMesher(@Rectangle,N,100);    
    
elseif (strcmp(MeshType,'C')==1)
    
    % cartesian mesh
    ax = Dati.domain(1); bx = Dati.domain(2);
    ay = Dati.domain(3); by = Dati.domain(4);
    nelx = N(1); nely = N(2);
    Nel = nelx*nely;
    
    dx = (bx-ax)/nelx; dy = (by-ay)/nely;
    [X,Y] = meshgrid(ax+dx/2:dx:bx, ay+dy/2:dy:by);
    P = [X(:) Y(:)];
    
    [region] = PolyMesher(@Rectangle,Nel,1,P);

elseif (strcmp(MeshType,'T')==1)

    % structured triangular mesh
    ax = Dati.domain(1); bx = Dati.domain(2);
    ay = Dati.domain(3); by = Dati.domain(4);
    nelx = N(1); nely = N(2);
    Nel = 2*nelx*nely;
    dx = (bx-ax)/nelx; dy = (by-ay)/nely;
    [X,Y] = meshgrid(ax:dx:bx, ay:dy:by);
    P = [X(:) Y(:)];

    nrows = size(X, 1);
    ncols = size(X, 2);

    % Initialize list of incenters
    incenters = [];

    % Loop through each rectangle
    for i = 1:(nrows-1)
        for j = 1:(ncols-1)
            orientation = mod(mod(i, 2) + mod(j, 2), 2);
            % Define vertices of the rectangle
            p1 = [X(i, j), Y(i, j)];
            p2 = [X(i+1, j), Y(i+1, j)];
            p3 = [X(i, j+1), Y(i, j+1)];
            p4 = [X(i+1, j+1), Y(i+1, j+1)];
            
            if orientation == 0
                tri1 = [p1; p2; p3];
                tri2 = [p2; p4; p3];
            else
                tri1 = [p1; p4; p3];
                tri2 = [p1; p2; p4];
            end
            
            % Triangle 1: p1, p2, p3
            
            incenter1 = computeIncenter(tri1);
            
            % Triangle 2: p2, p4, p3
            
            incenter2 = computeIncenter(tri2);
            
            % Append results
            incenters = [incenters; incenter1; incenter2];
        end
    end

    [region] = PolyMesher(@Rectangle,Nel,1,incenters);
end



%% Compute the neighbor structure
[region,neighbor] = MakeNeighbor(Data,region,SimType);

%% Otuput 
FileNameOut = strcat(FolderName,'/',FileName,'_',num2str(region.ne),'_el.mat');
save(FileNameOut,'region','neighbor');
end

% Helper function to compute incenter
function incenter = computeIncenter(triangle)
    % Vertices
    p1 = triangle(1, :);
    p2 = triangle(2, :);
    p3 = triangle(3, :);
    
    % Side lengths
    a = norm(p2 - p3);
    b = norm(p3 - p1);
    c = norm(p1 - p2);
    
    % Incenter formula
    x = (a * p1(1) + b * p2(1) + c * p3(1)) / (a + b + c);
    y = (a * p1(2) + b * p2(2) + c * p3(2)) / (a + b + c);
    
    % Return incenter
    incenter = [x, y];
end
