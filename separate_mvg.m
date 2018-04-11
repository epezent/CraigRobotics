% Evan Pezent | evanpezent.com | epezent@rice.edu
% 02/11/2016

function [M,V,G] = separate_mvg(MVG,Qdd,g)
% =========================================================================
% Extracts the mass matrix M, the vector of centrifugal and Coriolis terms
% V, and the vector of gravity terms G from the clumped symbolic expression
% MVG. Qdd is the [n x 1] vector of joint accelerations and g is gravity.
% =========================================================================
n = length(MVG);
M = sym(zeros(n,n));
G = sym(zeros(n,1));
for i = 1:n
    for j = 1:n
        M(i,j) = simplify( diff(MVG(i),Qdd(j)) );
    end
    G(i,1) = simplify( diff(MVG(i),g) * g );
end
V = simplify( MVG - M*Qdd - G );
end
