%> @file  MatrixLaplacianST.m
%> @author The Lymph Team
%> @date 16 April 2023
%> @brief Assembly of the matrices for the Poisson's problem
%> (subtriangulation)
%>
%==========================================================================
%> @section classMatrixLaplacianST Class description
%==========================================================================
%> @brief            Assembly of the mass and stifness matrices for the
%> Poisson's problem (subtriangulation)
%
%> @param Data       Struct with problem's data
%> @param neighbor   Neighbor struct (see MakeNeighbor.m)
%> @param femregion  Finite Element struct (see CreateDOF.m)
%
%> @retval Matrices  Matrices.Mprj = Mass Matrix;
%>                   Matrices.A    = Stiffness matrix;
%>                   Matrices.dGA  = dG matrix for error analysis;
%>
%==========================================================================

function [Matrices] = MatrixLaplacianST(Data, neighbor, femregion)

%% Quadrature values

[ref_qNodes_1D, w_1D, ref_qNodes_2D, w_2D] = Quadrature(femregion.nqn);

%% Initialization of the matrices
max_nedges = max(neighbor.nedges+1);

ii_index = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel),femregion.nbases,femregion.nbases,ones(1,femregion.nel));
jj_index = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel),femregion.nbases,femregion.nbases,ones(1,femregion.nel));
ii_index_neigh = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel));
jj_index_neigh = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel));

% \int_{\Omega}  ( u.v ) dx
Mprj_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel),femregion.nbases,femregion.nbases,ones(1,femregion.nel));

% \int_{\Omega} ( D grad(u).grad(v) ) dx
A_loc = mat2cell(zeros(femregion.nbases,femregion.nbases,femregion.nel),femregion.nbases,femregion.nbases,ones(1,femregion.nel));

% \int_{E_h} ( {D grad(u)}.[v] ) ds
IA_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel));    % DG-Stiffness Matrix

% \int_{E_h} penalty_E ( h_s^(-1) [u].[v] ) ds
SA_loc = mat2cell(zeros(femregion.nbases*max_nedges,femregion.nbases,femregion.nel),femregion.nbases*max_nedges,femregion.nbases,ones(1,femregion.nel));    % DG-Stiffness Matrix


%% Loop over the elements

% Visualization of computational progress
prog = 0;
fprintf(1,'Computation Progress: %3d%%\n',prog);


IAn_loc_all = cell(femregion.nel, 1);
SAn_loc_all = cell(femregion.nel, 1);

for ie = 1:femregion.nel

    % Visualization of computational progress
    prog = ( 100*(ie/femregion.nel) );
    fprintf(1,'\b\b\b\b%3.0f%%',prog);

    % Selection of the matrix positions associated to element ie
    index = (ie-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';
    ii_index{ie} = repmat(index, 1,femregion.nbases);
    jj_index{ie} = repmat(index',femregion.nbases,1);

    % Extraction of neighbor element and their edges
    neigh_ie      = neighbor.neigh{ie};
    neigh_ie_unq  = unique(neighbor.neigh{ie});
    neighedges_ie = neighbor.neighedges{ie};

    % Extraction of element geometrical information
    coords_ie          = femregion.coords_element{ie};
    [normals,meshsize] = GetNormalsMeshSizeFaces(coords_ie);

    % Creation of the subtriangulation of the element
    edges    = [1:femregion.nedges(ie) ; 2:femregion.nedges(ie), 1]';
    Tria_Del = delaunayTriangulation(coords_ie(:,1),coords_ie(:,2), edges);
    Tria     = Tria_Del( isInterior(Tria_Del) == 1, :);

    % Loop over the subtriangulation
    for iTria = 1:size(Tria,1)

        % Construction of Jacobian and quadrature nodes
        [BJ, qNodes_2D] = GetJacobianPhysicalPoints(coords_ie(Tria(iTria,:),:), ref_qNodes_2D);

        xq  = qNodes_2D(:,1);
        yq  = qNodes_2D(:,2);

        % Scaled weights
        dx = det(BJ) * w_2D;

        % Evaluation of physical parameters
        mu = Data.mu{1}(xq,yq);

        % Construction and evalutation on the quadrature points of the basis functions
        [phiq, gradqx, gradqy] = Evalshape2D(femregion, ie, qNodes_2D);

        %% Matrix assembling

        Mprj_loc{ie} = Mprj_loc{ie} + (dx.*phiq)'*phiq;
        A_loc{ie}    = A_loc{ie} + (dx .* (mu .* gradqx))' * gradqx + (dx .* (mu .* gradqy))' * gradqy;

    end

    %% Boundary integrals and stabilization terms

    % Computation of all the penalty coefficients for the element ie
    [penalty_geom] = PenaltyCoefficient(femregion, Data, ie, neighedges_ie, neigh_ie, meshsize);

    IAn_loc = zeros(femregion.nbases, femregion.nbases, neighbor.nedges(ie));
    SAn_loc = zeros(femregion.nbases, femregion.nbases, neighbor.nedges(ie));

    % Loop over faces
    iedg_real = 1;
    for iedg = 1 : neighbor.nedges(ie)
        % Extraction of the id of the neighboring element in the matrices IAN and SAN
        idneigh = (neigh_ie_unq == neighbor.neigh{ie}(iedg));

        % Extraction of the indexes for assembling face matrices (contribution of the element ie)
        ii_index_neigh{ie}(1:femregion.nbases,:) = repmat(index, 1 ,femregion.nbases);
        jj_index_neigh{ie}(1:femregion.nbases,:) = repmat(index',femregion.nbases,1);

        % Extraction of the edge coordinates
        if iedg == neighbor.nedges(ie)
            p1 = coords_ie(iedg,:);
            p2 = coords_ie(1,:);
        else
            p1 = coords_ie(iedg,:);
            p2 = coords_ie(iedg+1,:);
        end

        % Construction of quadrature nodes on the face
        [qNodes_1D] = GetPhysicalPointsFaces([p1; p2], ref_qNodes_1D);

        xq = qNodes_1D(:,1);
        yq = qNodes_1D(:,2);

        % Scaled weights
        ds = meshsize(iedg) * w_1D;

        % Evaluation of physical parameters
        mu    = Data.mu{1}(xq,yq);

        % Construction and evalutation on the quadrature points of the basis functions
        [phiedgeq, gradedgeqx, gradedgeqy] = Evalshape2D(femregion, ie, qNodes_1D);

        % Extraction of normals to the face
        nx = normals(1,iedg);
        ny = normals(2,iedg);

        %% Matrix assembling

        % Dirichlet boundary faces
        if neigh_ie(iedg) == -1

            IA_loc{ie}(1:femregion.nbases,:) = IA_loc{ie}(1:femregion.nbases,:) + (ds .* mu .* (nx * gradedgeqx + ny * gradedgeqy))' * phiedgeq;
            SA_loc{ie}(1:femregion.nbases,:) = SA_loc{ie}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* mu .* phiedgeq)' * phiedgeq;

            % Internal faces
        elseif neigh_ie(iedg) > 0

            % Element itself
            IA_loc{ie}(1:femregion.nbases,:) = IA_loc{ie}(1:femregion.nbases,:) + 0.5 * (ds .* mu .* (nx * gradedgeqx + ny * gradedgeqy))' * phiedgeq;
            SA_loc{ie}(1:femregion.nbases,:) = SA_loc{ie}(1:femregion.nbases,:) + penalty_geom(iedg) * (ds .* mu .* phiedgeq)' * phiedgeq;

            % Construction and evalutation on the quadrature points of the basis functions for the neighbor
            phiedgeqneigh = Evalshape2D(femregion, neigh_ie(iedg), qNodes_1D);

            % Neighboring element
            neigh_idx = find(idneigh)*femregion.nbases+1:(find(idneigh)+1)*(femregion.nbases);
            index_neigh = (neighbor.neigh{ie}(iedg)-1)*femregion.nbases*ones(femregion.nbases,1) + (1:femregion.nbases)';

            % Extraction of the indexes for assembling face matrices (contribution of the neighboring element)
            ii_index_neigh{ie}(neigh_idx,:) = repmat(index, 1,femregion.nbases);
            jj_index_neigh{ie}(neigh_idx,:) = repmat(index_neigh',femregion.nbases,1);

            % Extracellular component
            IA_loc{ie}(neigh_idx,:)  = IA_loc{ie}(neigh_idx,:) - 0.5 * (ds .* mu .* ( nx * gradedgeqx +  ny * gradedgeqy))' * phiedgeqneigh;
            SA_loc{ie}(neigh_idx,:)  = SA_loc{ie}(neigh_idx,:) - penalty_geom(iedg) * (ds .* mu .* phiedgeq)' * phiedgeqneigh;
            
            IAn_loc(:, :, iedg_real) = IAn_loc(:, :, iedg_real) + (ds .* mu .* ( nx * gradedgeqx +  ny * gradedgeqy))' * phiedgeqneigh;
            SAn_loc(:, :, iedg_real) = SAn_loc(:, :, iedg_real) + penalty_geom(iedg) * (ds .* mu .* phiedgeq)' * phiedgeqneigh;
            iedg_real = iedg_real + 1;
        end

    end
    IAn_loc_all{ie} = IAn_loc; 
    SAn_loc_all{ie} = SAn_loc; 
    writematrix(IAn_loc, sprintf("IAn_loc_%d.txt", ie - 1))
    writematrix(SAn_loc, sprintf("SAn_loc_%d.txt", ie - 1))
end



% Local matrix to global matrix
ii_index  = reshape(cell2mat(ii_index),[femregion.nbases,femregion.nbases*femregion.nel]);
jj_index  = reshape(cell2mat(jj_index),[femregion.nbases,femregion.nbases*femregion.nel]);

Mprj_loc = reshape(cell2mat(Mprj_loc),[femregion.nbases,femregion.nbases*femregion.nel]);
A_loc = reshape(cell2mat(A_loc),[femregion.nbases,femregion.nbases*femregion.nel]);
writematrix(Mprj_loc)
writematrix(A_loc)

IA_loc_self = zeros(size(A_loc));
SA_loc_self = zeros(size(A_loc));
for ie = 1:femregion.nel
    j0 = (ie - 1) * femregion.nbases + 1;
    IA_loc_self(:, j0:j0+femregion.nbases-1) = IA_loc{ie}(1:femregion.nbases,:);
    SA_loc_self(:, j0:j0+femregion.nbases-1) = SA_loc{ie}(1:femregion.nbases,:);
end
writematrix(IA_loc_self)
writematrix(SA_loc_self)

ii_index_neigh = reshape(cell2mat(ii_index_neigh),[femregion.nbases,femregion.nel*max_nedges*femregion.nbases]);
jj_index_neigh = reshape(cell2mat(jj_index_neigh),[femregion.nbases,femregion.nel*max_nedges*femregion.nbases]);

IA_loc  = reshape(cell2mat(IA_loc),[femregion.nbases,femregion.nel*max_nedges*femregion.nbases]);
SA_loc  = reshape(cell2mat(SA_loc),[femregion.nbases,femregion.nel*max_nedges*femregion.nbases]);

del = all(jj_index_neigh == 0,1);
ii_index_neigh(:,del) = [];
jj_index_neigh(:,del) = [];

IA_loc(:,del) = [];
SA_loc(:,del) = [];

Mprj = sparse(ii_index,jj_index,Mprj_loc,femregion.ndof,femregion.ndof);
A    = sparse(ii_index,jj_index,A_loc,femregion.ndof,femregion.ndof);

IA   = sparse(ii_index_neigh,jj_index_neigh,IA_loc,femregion.ndof,femregion.ndof);
SA   = sparse(ii_index_neigh,jj_index_neigh,SA_loc,femregion.ndof,femregion.ndof);

%% Symmetric Interior Penalty Method

dGA = A + SA; % For error computation in dG-norm

A = dGA - IA - transpose(IA); % dG stiffness matrix
writematrix(full(A), 'A_glob.txt')

% Build the Matrices structure
Matrices = struct('Mprj', Mprj, ...
    'A',    A,    ...
    'dGA',  dGA);

end
