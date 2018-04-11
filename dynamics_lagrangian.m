function [Tau, L, k, u, K, U] = Lagrangian(m,Pc,Ic,T_array,Q,Qd,Qdd,g0,u_ref)

num = length(m);
vc = cell(num,1);
w  = cell(num,1);

% Local X vector
Z = [0; 0; 1];

for i = 1:num
    R = T_array{i}(1:3,1:3).';
    T0i = eye(4);
    for j = 1:i
        T0i = simplify(T0i*T_array{j});
    end
    P0ci = T0i*[Pc{i}(:);1];
    P0ci = P0ci(1:3,1);
    vc{i} = simplify(jacobian(P0ci,Q)*Qd);
    if i == 1
        w{i} = Qd(i)*Z;
    else
        w{i} = simplify(R*w{i-1} + Qd(i)*Z);
    end
    k(i) = simplify(1/2 * m(i) * vc{i}.' * vc{i} + 1/2 * w{i}.' * Ic{i} * w{i});
    u(i) = simplify(-m(i)*g0.'*P0ci + u_ref);
end

% Total Kinetic Engery
K = simplify(sum(k));
% Total Potential Energy
U = simplify(sum(u));
% Lagrangian
L = simplify(K-U);

% Compute Derivatives
for i = 1:num
    dLdQ(i,1)  = simplify( diff(L, Q(i)) );
    dLdQd(i,1) = simplify( diff(L, Qd(i)) );
    ddtdLdQd(i,1) = sym(0);
    for j = 1:num
        ddtdLdQd(i,1) = ddtdLdQd(i) + diff(dLdQd(i),Qd(j))*Qdd(j) + diff(dLdQd(i),Q(j))*Qd(j);
    end
    ddtdLdQd(i,1) = simplify(ddtdLdQd(i));
end

% compute torques
Tau = simplify(ddtdLdQd - dLdQ);

end
