% Evan Pezent | evanpezent.com | epezent@rice.edu
% 02/04/2017

function linear_sym = linearize_dynamics(nonlinear_sym, x_sym, x0)

linear_sym = subs(nonlinear_sym,x_sym,x0);

for i = 1:length(x_sym)
    linear_sym = linear_sym + subs(diff(nonlinear_sym,x_sym(i)), x_sym, x0)*(x_sym(i)-x0(i));
end

linear_sym = simplify(expand(linear_sym));

end

