% Evan Pezent | evanpezent.com | epezent@rice.edu
% 02/04/2017

function [Tau,w,wd,vd,vcd,F,N,f,n] = dynamics_newtonian(m,Pc,Ic,Ti,Qd,Qdd,g0)
% =========================================================================
% Computes the dynamic equations of motion for a rotational robotic 
% manipulator using the iterative Newton-Euler formulation.
%
% m = [n x 1] vector or link masses 
% Pc = n length cell array of [3 x 1] translations from {i} to {c_i} 
% Ic = n length cell array of [3 x 3] inertia tensors take about {c_i}
% Ti = Ti as obtained from function dh2tf(DH_table)
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
% Eqns: 6.45 - 6.53 (pg. 176)
% =========================================================================

num = length(m);

% Pad vectors for notation consistency
m = [0; m];
Pc = [0 Pc];
Ic = [0 Ic];
Qd = [0; Qd];
Qdd = [0; Qdd];


% Local X vector
Z = [0; 0; 1];

% Frame 0 Variables
w{1} = [0 0 0].';
wd{1}= [0 0 0].';

% Gravity Orientation
% "The effect of gravity loading on the links can be included quite simply 
% by setting vd0 = G, where G has the magnitude of the gravity vector but 
% points in the opposite direction. This is equivalent to saying that the 
% base of the robot is accelerating upward with 1 g acceleration. This 
% fictitious upward acceleration causes exactly the same effect on the 
% links as gravity would." - pg. 176
G = -g0; 
vd{1} = G;

%% Outward Iterations
for i = 1:num % 0 -> n-1
    R        = Ti{i}(1:3,1:3).'; % ^i+1_i R
    P        = Ti{i}(1:3,4); % ^i P_i+1
    w{i+1}   = R*w{i} + Qd(i+1)*Z; % 6.45
    wd{i+1}  = R*wd{i} + cross(R*w{i},Qd(i+1)*Z) + Qdd(i+1)*Z; % 6.46
    vd{i+1}  = R*(cross(wd{i},P) + cross(w{i},cross(w{i},P)) + vd{i}); % 6.47
    vcd{i+1} = cross(wd{i+1},Pc{i+1}) + cross(w{i+1},cross(w{i+1},Pc{i+1})) + vd{i+1}; % 6.48
    F{i+1}   = m(i+1)*vcd{i+1}; % 6.49
    N{i+1}   = Ic{i+1}*wd{i+1} + cross(w{i+1},Ic{i+1}*w{i+1}); % 6.50
end

%% Inward Iterations
for i = num+1:-1:2 % n -> 1
    if i == num+1
        f{i} = F{i}; % 6.51
        n{i} = N{i} + cross(Pc{i},F{i}); % 6.52
    else
        R = Ti{i}(1:3,1:3);
        P = Ti{i}(1:3,4);
        f{i} = R*f{i+1} + F{i}; % 6.51
        n{i} = N{i} + R*n{i+1} + cross(Pc{i},F{i}) + cross(P,R*f{i+1}); % 6.52
    end
    Tau(i,1) = simplify( n{i}.'*Z ); % 6.53
end

%% Clean up elements related to 0th frame
Tau(1) = [];  
w(:,1) = [];
wd(:,1) = [];
vd(:,1) = [];
vcd(:,1) = [];
F(:,1) = [];
N(:,1) = [];
f(:,1) = [];
n(:,1) = [];       
        
end
    
    