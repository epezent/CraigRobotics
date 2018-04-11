% Evan Pezent | evanpezent.com | epezent@rice.edu
% 04/10/2018

function [Tau, L, K, U, k, u] = dynamics_lagrangian(m,Pc,Ic,Ti,Q,Qd,Qdd,g0,u_ref)
% =========================================================================
% Computes the dynamic equations of motion for a robotic manipulator using 
% the Lagrangian formulation.
%
% m = [n x 1] vector or link masses 
% Pc = n length cell array of [3 x 1] translations from {i} to {c_i} 
% Ic = n length cell array of [3 x 3] inertia tensors take about {c_i}
% T_array = T_array as obtained from function DH2TF(DH_table)
% Qd = [n x 1] vector of joint angular velocities
% Qdd = [n x 1] vector of joint angular accelerations
% g0 = [3 x 1] gravity vector in {0}. 
%
% where
%
% {i} is the frame attached to link i
% {c_i} is the frame at the COM of link i with the same orientation as {i}
% =========================================================================
% Source: Introduction to Robotics: Mechanics and Control (3e) - Craig, J.
% Eqns: 5.64 (pg. 150), 6.69 - 6.76 (pg. 182 - 182)
% =========================================================================

num = length(m);

% Local X vector
Z = [0; 0; 1];

for i = 1:num
    % compute rotation and transformation needs
    R = Ti{i}(1:3,1:3).'; % ^i+1_i R
    T0i = eye(4);
    for j = 1:i
        T0i = simplify(T0i*Ti{j});  % ^0_i T
    end
    % center of mass of link i in frame 0
    P0ci = T0i*[Pc{i}(:);1];
    P0ci = P0ci(1:3,1);
    % linear velocity of link i at center of mass
    vc{i} = simplify(jacobian(P0ci,Q)*Qd); % 5.64
    % angular velocity of link i
    if i == 1
        w{i} = Qd(i)*Z;
    else
        w{i} = simplify(R*w{i-1} + Qd(i)*Z);
    end
    % kinetic energy of link i
    k(i,1) = simplify(1/2 * m(i) * vc{i}.' * vc{i} + ...
        1/2 * w{i}.' * Ic{i} * w{i}); % 6.69
    % potentianl energy of link i
    u(i,1) = simplify(-m(i)*g0.'*P0ci + u_ref); % 6.73
end

% Total Kinetic Engery
K = simplify(sum(k)); % 6.70
% Total Potential Energy
U = simplify(sum(u)); % 6.74
% Lagrangian
L = simplify(K-U);

% Compute Derivatives
for i = 1:num
    dLdQ(i,1)  = simplify( diff(L, Q(i)) );  % dL/dQ
    dLdQd(i,1) = simplify( diff(L, Qd(i)) ); % dL/dQd
    ddtdLdQd(i,1) = sym(0);                  
    for j = 1:num
        ddtdLdQd(i,1) = ddtdLdQd(i) + diff(dLdQd(i),Qd(j))*Qdd(j) + diff(dLdQd(i),Q(j))*Qd(j);
    end
    ddtdLdQd(i,1) = simplify(ddtdLdQd(i));   % d/dt(dL/dQd)
end

% compute torques / equations of motion
Tau = simplify(ddtdLdQd - dLdQ); % 6.76

end
