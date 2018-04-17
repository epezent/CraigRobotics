clear all
close all
clc

% Evan Pezent
% 4/17/2018

% Exercise 6.5/6.19 of Craig

addpath('../')

%% Problem Setup
syms m1 m2 l1 l2 q1 q2 q1d q2d q1dd q2dd tau1 tau2 g

m  = [m1;m2];
Pc = {[l1;0;0], [l2;0;0]}; % point mass at the distal end of each link
Ic = {zeros(3), zeros(3)}; % zero inertia tensors due to point mass

%           a  alpha   d  theta
DH_table = [0    0     0   q1;   % i = 1
            l1   pi/2  0   q2];  % i = 2
            % note: we don't need the last row
        
[~,Ti] = dh2tf(DH_table); % {T01, T12}

g0 = [0;0;-g]; % gravity acts along negative Z0

Q   = [q1;q2];
Qd  = [q1d;q2d];
Qdd = [q1dd;q2dd];
Tau = [tau1; tau2];

%% Dynamics (Newtonian)
Tau_newtonian = dynamics_newtonian(m,Pc,Ic,Ti,Qd,Qdd,g0);

%% Dynamics (Lagrangian)
Tau_lagrangian = dynamics_lagrangian(m,Pc,Ic,Ti,Q,Qd,Qdd,g0,0)

%% Check Dynamics
Tau_check = simplify(Tau_newtonian - Tau_lagrangian); % [0;0]

%% Separate into M, V, and G
[M,V,G] = separate_mvg(Tau_lagrangian,Qdd,g)

% Note: Results may be factored differently from textbook solutions, but
% they are equal.

M_book = [m1*l1^2 + m2*(l1 + l2*cos(q2))^2      0; 
          0                               m2*l2^2];
      
V_book = [-2*(l1 + l2*cos(q2))*m2*l2*sin(q2)*q1d*q2d; 
          (l1+l2*cos(q2))*m2*l2*sin(q2)*q1d^2];
      
G_book = [0;
          m2*g*l2*cos(q2)];
      
%% Check M,V,G against book solutions
M_check = simplify(M - M_book); % [0 0;0 0]
V_check = simplify(V - V_book); % [0;0]
G_check = simplify(G - G_book); % [0;0]

%% Solve for Qdd (for dynamic simulation)
Qdd_solved = simplify( inv(M) * (Tau - V - G) ) % 6.115