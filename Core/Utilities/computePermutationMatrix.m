function P = computePermutationMatrix(order)
    ndof_loc = (order + 1) * (order + 2) / 2;
    idx = 1;
    p = zeros(ndof_loc, 1);

    for i = 0:order
        for j = 0:order
            if i + j > order
                continue;
            end
            p(idx) = getLinearIndex2D(i, j);
            idx = idx + 1;
        end
    end

    P = eye(ndof_loc);
    P = P(p + 1, :);  % MATLAB uses 1-based indexing
end

function idx = getLinearIndex2D(i, j)
    % Assuming monomials ordered by total degree then lexicographic
    idx = (i + j) * (i + j + 1) / 2 + j;
end