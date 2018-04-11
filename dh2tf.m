% Evan Pezent | evanpezent.com | epezent@rice.edu
% 02/11/2016

function [T0N,Ti,T0i] = dh2tf(DH_table)
% =========================================================================
% Accepts an Nx4 matrix with rows of the form [a alpha d theta] which
% correspond to the DH parameters of subsequent frames i=1:N, where:
%
% a     = the distance from Z_i-1 to Z_i measured along X_i-1
% alpha = the angle from Z_i-1 to Z_i measured about X_i-1
% d     = the distance from X_i-1 to X_i measured along Z_i
% theta = the angle from X_i-1 to X_i measured about Z_i
%
% The first row should be frame 1, the second row frame 2, and so on until
% the final frame N. Returns the transformation from frame N to 0, an
% N-array of transformations each describing frame i relative to frame i-1,
% i.e. mapping i to i-1, and and N-array of transformations each describing
% frame i relative to frame 0.
% =========================================================================
% Source: Introduction to Robotics: Mechanics and Control (3e) - Craig, J.
% Eqns: 3.6 (pg. 75)
% =========================================================================
% EXAMPLE 3.6:
%
% DH_table = [0 0 0 theta1; 0 pi/2 d2 0; 0 0 L2 theta3];
% [T03,Ti] = dh2tf(DH_table)
%
% T03 = [ cos(theta1)*cos(theta3), -cos(theta1)*sin(theta3),  sin(theta1),   L2*sin(theta1) + d2*sin(theta1)]
%       [ cos(theta3)*sin(theta1), -sin(theta1)*sin(theta3), -cos(theta1), - L2*cos(theta1) - d2*cos(theta1)]
%       [             sin(theta3),              cos(theta3),            0,                                 0]
%       [                       0,                        0,            0,                                 1]
%
% T01 = Ti{1}
% [ cos(theta1), -sin(theta1), 0, 0]
% [ sin(theta1),  cos(theta1), 0, 0]
% [           0,            0, 1, 0]
% [           0,            0, 0, 1]
%
% T12 = Ti{2}
% [ 1, 0,  0,   0]
% [ 0, 0, -1, -d2]
% [ 0, 1,  0,   0]
% [ 0, 0,  0,   1]
%
% T23 = Ti{3}
% [ cos(theta3), -sin(theta3), 0,  0]
% [ sin(theta3),  cos(theta3), 0,  0]
% [           0,            0, 1, L2]
% [           0,            0, 0,  1]
% =========================================================================

Ti  = cell(1,size(DH_table,1));
T0i = cell(1,size(DH_table,1));
T0N = eye(4);

for ii = 1:size(DH_table,1)
    a = DH_table(ii,1);
    alpha = DH_table(ii,2);
    d = DH_table(ii,3);
    theta = DH_table(ii,4);
    sinTheta = sin(theta);
    cosTheta = cos(theta);
    sinAlpha = sin(alpha);
    cosAlpha = cos(alpha);
    T = [cosTheta -sinTheta 0 a;
        sinTheta*cosAlpha cosTheta*cosAlpha -sinAlpha -sinAlpha*d;
        sinTheta*sinAlpha cosTheta*sinAlpha cosAlpha cosAlpha*d;
        0 0 0 1];
    Ti{ii} = T;
    T0N = T0N*T;
    T0i{ii} = T0N;
end

end
