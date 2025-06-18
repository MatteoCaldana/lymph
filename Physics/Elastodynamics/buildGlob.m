function [M] = buildGlob(M1, M2, M3, M4, degree)

b = (degree + 1) * (degree + 2) / 2;
Perm = computePermutationMatrix(degree);

[~, ngdof] = size(M1);

starts = 1:2*b:(2*ngdof);  
offsets = (0:b-1)';        
indices = starts + offsets; 
indices = Perm * indices;
indices = indices(:);
M(indices, indices) = M1;
M(indices + b, indices + b) = M4;
M(indices, indices+b) = M2;
M(indices+b, indices) = M3;
end