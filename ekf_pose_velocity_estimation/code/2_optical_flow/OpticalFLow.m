%% PROJECT 2 VELOCITY ESTIMATION
close all;
clear all;
clc;
addpath('../data')

%Change this for both dataset 1 and dataset 4. Do not use dataset 9.
datasetNum = 4;

[sampledData, sampledVicon, sampledTime] = init(datasetNum);

%% INITIALIZE CAMERA MATRIX AND OTHER NEEDED INFORMATION

K = [311.0520, 0, 201.8724; 0, 311.3885, 113.6210; 0, 0, 1];
tracker = vision.PointTracker('MaxBidirectionalError', 3);
estimatedV = zeros(6, length(sampledData));
e = 0.8; %ransac threshold
useRANSAC = true;

for n = 2:length(sampledData)

    %% Initalize Loop load images
    data_prev = sampledData(n - 1);
    data_cur = sampledData(n);

    %% Detect good points
    corners_prev = detectMinEigenFeatures(data_prev.img);

    %% Initalize the tracker to the last frame.
    initialize(tracker, corners_prev.Location, data_prev.img);

    %% Find the location of the next points;
    [corners_cur, validity_cur] = tracker(data_cur.img);

    %% Calculate velocity
    % Use a for loop
    num_pts = size(corners_cur, 1);
    dt = sampledData(n).t - sampledData(n - 1).t;

    u = [];
    v = [];
    cam_pts = corners_cur;
    cam_pts_prev = corners_prev.Location;

    x = [];
    y = [];
    x_prev = [];
    y_prev = [];

    for i = 1:num_pts

        temp = [cam_pts(i, 1) cam_pts(i, 2) 1];
        cam_xy = K \ temp';
        temp2 = [cam_pts_prev(i, 1) cam_pts_prev(i, 2) 1];
        cam_xy_prev = K \ temp2';

        % transform the pixel coordinates into camera frame
        x(i) = cam_xy(1);
        y(i) = cam_xy(2);
        x_prev(i) = cam_xy_prev(1);
        y_prev(i) = cam_xy_prev(2);
    end

    cam_pts = [x', y'];
    cam_pts_prev = [x_prev', y_prev'];

    for i = 1:num_pts
        u = [u; (cam_pts(i, 1) - cam_pts_prev(i, 1)) / dt];
        v = [v; (cam_pts(i, 2) - cam_pts_prev(i, 2)) / dt];
    end

    %% Calculate Height

    [T_b2w, ang_b2w, R_c2w] = estimatePose(sampledData, n);
    R_b2w = eul2rotm(ang_b2w);
    H_b2w = [R_b2w, T_b2w; 0, 0, 0, 1];
    R_b2c = eul2rotm([-pi / 4, 0, -pi]);
    T_b2c = [-0.04; 0; -0.03];
    H_c2b = [R_b2c', -R_b2c' * T_b2c; 0, 0, 0, 1];
    H_c2w = H_b2w * H_c2b;
    H = H_c2w;
    num_npts = size(cam_pts, 1);
    Z = zeros(num_npts, 1);

    %% Calculating Height
    for i = 1:num_npts
        Z(i) =- (cam_pts(i, 1) * H(3, 1) + cam_pts(i, 2) * H(3, 2) + H(3, 4)) / (H(3, 3));
    end

    optV = [u, v];
    %% RANSAC
    % Write your own RANSAC implementation in the file velocityRANSAC
    if useRANSAC
        Vel = velocityRANSAC(optV, cam_pts, Z, R_c2w, e);
    else
        j = size(cam_pts, 1);
        x = cam_pts(:, 1);
        y = cam_pts(:, 2);
        H = zeros(2 * j, 6);

        for i = 1:j
            A = (1 / Z(i)) * [-1, 0, x(i); 0, -1, y(i)];
            B = [x(i) * y(i), - (1 + x(i)^2), y(i); 1 + y(i)^2, -x(i) * y(i), -x(i)];
            H(2 * i - 1:2 * i, :) = [A, B];
        end

        opt_V = [optV(:, 1); optV(:, 2)];
        Vel = pinv(H) * opt_V;
    end

    %% Fix the linear velocity
    % Change the frame of the computed velocity to world frame
    T_c2b = -R_b2c' * T_b2c;
    r12_c = R_c2w' * T_c2b;
    skew_r12 = [0 -r12_c(3) r12_c(2); r12_c(3) 0 -r12_c(1); -r12_c(2) r12_c(1) 0];
    adjoint = [R_b2c' -R_b2c' * skew_r12; zeros(3) R_b2c'];
    Vel_b = adjoint * Vel;
    adjoint = [R_b2w zeros(3); zeros(3) R_b2w];

    %% Thereshold outputs into a range.
    % Not necessary

    %% ADD SOME LOW PASS FILTER CODE
    % Not neceessary but recommended
    % estimatedV(:,n) = Vel;

    %% STORE THE COMPUTED VELOCITY IN THE VARIABLE estimatedV AS BELOW
    % Structure of the Vector Vel should be as follows:
    % Vel(1) = Linear Velocity in X
    % Vel(2) = Linear Velocity in Y
    % Vel(3) = Linear Velocity in Z
    % Vel(4) = Angular Velocity in X
    % Vel(5) = Angular Velocity in Y
    % Vel(6) = Angular Velocity in Z
    estimatedV(:, n) = adjoint * Vel_b;

    release(tracker);
end

plotData(estimatedV, sampledData, sampledVicon, sampledTime, datasetNum)
