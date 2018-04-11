% Evan Pezent
% 02/09/2016
% MECH 498 - Intro to Robotics
% Lab 1 - Homogenous Transforms and Robot Foward Kinematics
% =========================================================================
% screwDH.m
% =========================================================================
% Using your screwtf function, write a function screwDH(a,alpha,d,theta) 
% that accepts four scalar DH parameters and returns a 4x4 homogeneous 
% transformation matrix according to the DH convention. You can check your
% result using the formula for this matrix given in the book and in class.
% Works with numeric and symbolic values.
% =========================================================================
% Eqns: 3.5, 3.6 (pg.75)
% =========================================================================

function T = screwDH(a,alpha,d,theta)

% 3.5 (pg.75)
% Note: [1 0 0] corresponds to X-axis, [0 0 1] corresponds to Z-axis.
T = screwtf(a,alpha,[1 0 0])*screwtf(d,theta,[0 0 1]);

% 3.6 (pg.75) (for checking above only)
% T_check = [cos(theta) -sin(theta) 0 a;
%      sin(theta)*cos(alpha) cos(theta)*cos(alpha) -sin(alpha) -sin(alpha)*d;
%      sin(theta)*sin(alpha) cos(theta)*sin(alpha) cos(alpha) cos(alpha)*d;
%      0 0 0 1];

end