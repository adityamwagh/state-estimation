clear;

% position - x1
syms px py pz;

% euler angles - orientation - x2
syms x y z;

% pdot - x3
syms vx vy vz; % x3

% gyro bias - x4
syms bgx bgy bgz;

% acc bias - x5
syms bax bay baz;

% gyro measurement - body
syms wmx wmy wmz;

%  acc measurement - body
syms amx amy amz;

% gyro noise
syms ngx ngy ngz;

% acc noise
syms nax nay naz;

% gyro drift
syms nbgx nbgy nbgz;

% acc drift
syms nbax nbay nbaz;

% see enclosed jupyter notebook
% https://ethz.ch/content/dam/ethz/special-interest/mavt/robotics-n-intelligent-systems/rsl-dam/documents/RobotDynamics2016/KinematicsSingleBody.pdf
% G = equation 2.75, with with swapped x and z columns. Our orientationn from vicon is x,y,z, and not z,y,x.

G = [
    cos(y) * cos(z), -sin(z), 0;
    sin(z) * cos(y), cos(z), 0;
    -sin(y), 0, 1;
    ];

G_inv = simplify(inv(G));

% see enclosed jupyter notebook
R = [
    cos(y) * cos(z), sin(x) * sin(y) * cos(z) - sin(z) * cos(x), sin(x) * sin(z) + sin(y) * cos(x) * cos(z);
    sin(z) * cos(y), sin(x) * sin(y) * sin(z) + cos(x) * cos(z), -sin(x) * cos(z) + sin(y) * sin(z) * cos(x);
    -sin(y), sin(x) * cos(y), cos(x) * cos(y); ];

% assume that g is in m/s in -z
g = [0; 0; -9.81];

w = [
    (wmx - bgx - ngx);
    (wmy - bgy - ngy);
    (wmz - bgz - ngz);
    ];

a = [
    (amx - bax - nax);
    (amy - bay - nay);
    (amz - baz - naz);
    ];

xdot1 = [vx; vy; vz];
xdot2 = G_inv * R * w;
xdot3 = g + R * a;
xdot4 = [nbgx; nbgy; nbgz];
xdot5 = [nbax; nbay; nbaz];

% x
x1 = [
    px;
    py;
    pz;
    ];
x2 = [
    x;
    y;
    z;
    ];
x3 = [
    vx;
    vy;
    vz;
    ];
x4 = [
    bgx;
    bgy;
    bgz;
    ];
x5 = [
    bax;
    bay;
    baz;
    ];

    % x
x = vertcat(x1, x2, x3, x4, x5);

% computing process model
xdot = vertcat(xdot1, xdot2, xdot3, xdot4, xdot5);
n = [ngx; ngy; ngz; nax; nay; naz; nbgx; nbgy; nbgz; nbax; nbay; nbgz];

%  computing state matrices
A_t = simplify(jacobian(xdot, x));
U_t = simplify(jacobian(xdot, n));
xdot = simplify(xdot);
