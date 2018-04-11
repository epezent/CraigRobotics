% Evan Pezent | evanpezent.com | epezent@rice.edu
% 02/11/2016

function [M,V,G] = separate_mvg(MVG,Qdd,g)
% =========================================================================
% Extracts the mass matrix M, the vector of centrifugal and Coriolis terms
% V, and the vector of gravity terms G from the clumped symbolic expression
% MVG. Qdd is the [n x 1] vector of joint accelerations and g is gravity.
% =========================================================================

% this can be made a hell of a lot easier using diff(MVG(i), Qdd(j))

n = length(MVG);

% Extract M
for i = 1:n
    mvg = MVG(i);
    for j = 1:n
        m_temp = char(collect(mvg,Qdd(j)));
        ind = strfind(m_temp,char(Qdd(j)));
        if length(ind) < 2
            m = m_temp(1:ind-2);
            M(i,j) = simplify(expand(evalin('base',m)));
        else
            error(['Could not collect ' char(Qdd(j)) '.'])
        end
    end
end

% Reduce MVG to VG
VG = simplify(expand(MVG - M*Qdd));
for i = 1:n
    vg = VG(i);
    g_temp = char(collect(vg,g));
    ind = strfind(g_temp,char(g));
    if length(ind) < 2
        g_i = g_temp(1:ind-2);
    else
        error(['Could not collect ' char(g) '.'])
    end
    try
        G(i,1) = simplify(expand(evalin('base',g_i)))*g;
    catch
    end
end

% Reduce VG to V
V = simplify(expand(VG-G));

end
