%> @file  LegendreP.m
%> @author Mattia Corti, Paola F. Antonietti
%> @date 20 February 2023 
%> @brief The function evaluates the Legendre polynomials
%> 
%> The function evaluates the Legendre polynomial of order
%> \f$N\f$ at point \f$x\f$, namely \f$P_N(x)\f$. The polynomials are
%> hard-coded until the degree 7. Then we compute them by means of the recursive formula \cite handbookformulas :
%> \f[ P_{i+1}(x) = \frac{1}{i+1}\{(2i+1)xP_i(x)-iP_{i-1}(x)\} \f]
%>
%======================================================================
%> @section classLegendreP Class description
%======================================================================
%> @brief The function evaluates the Legendre polynomials.
%>
%> @param x     Vector of points \f$x_i\f$ in which the Legendre polynomials are evaluated.
%> @param N     Vector of Legendre polynomial orders \f$N_j\f$ to be evaluated.
%>
%> @retval P    Matrix containing at the position \f$P_{ij}\f$ a Legendre polynomial of order \f$N_j\f$ evaluated at the point \f$x_i\f$.
%======================================================================

function [P] = LegendreP(x, N)
    
    % Turn points into column if needed
    x = x(:); 

    % Initialization of the matrix P
    P = zeros(length(x), length(N));

    Nmax = max(N);
    Ptable = ones(length(x), Nmax + 1);
    Ptable(:, 2) = x;
    for i = 1:Nmax-1
        Ptable(:, i+2) = ((2*i+1)*x.*Ptable(:, i+1) - i*Ptable(:, i))/(i + 1);
    end
    for i = 1:numel(N)
        P(:, i) = Ptable(:, N(i) + 1);
    end

end
