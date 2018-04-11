% Evan Pezent
% 02/09/2016
% MECH 498 - Intro to Robotics
% Lab 1 - Homogenous Transforms and Robot Foward Kinematics
% =========================================================================
% screwtf.m
% =========================================================================
% Accepts two scalars, (translation) and (rotation), and a 3D unit vector 
% [ax] representing the axis of translation/rotation and returns the the 
% corresponding 4x4 homogeneous transformation matrix using axis-angle
% representation. Works with numeric and symbolic values.
% =========================================================================
% Eqns: 2.80 (pg.47), 2.19 (pg.28) 
% =========================================================================

function T = screwtf(translation, rotation, ax)

ax = ax(:);
ax = ax/norm(ax);

% 2.80 (pg.47)
kx = ax(1); ky = ax(2); kz = ax(3);
s = sin(rotation);
c = cos(rotation);
v = 1-c;

R = [kx*kx*v+c kx*ky*v-kz*s kx*kz*v+ky*s;
     kx*ky*v+kz*s ky*ky*v+c ky*kz*v-kx*s;
     kx*kz*v-ky*s ky*kz*v+kx*s kz*kz*v+c];

% 2.19 (pg.28)
T(1:3,1:3) = R;
T(1:3,4) = translation*ax;
T(4,4) = 1;
 
end
