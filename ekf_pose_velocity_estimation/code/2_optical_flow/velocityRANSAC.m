function [Vel] = velocityRANSAC(optV, optPos, Z, R_c2w, e)
    %% CHANGE THE NAME OF THE FUNCTION TO velocityRANSAC
    %% Input Parameter Description
    % optV = The optical Flow
    % optPos = Position of the features in the camera frame
    % Z = Height of the drone
    % R_c2w = Rotation defining camera to world frame
    % e = RANSAC hyper parameter

    %% Output Parameter Description
    % Vel = Linear velocity and angular velocity vector

    %% Get number of iterations
    % sample size
    M = 3;
    % required number of iterations with 3 samples
    k = ceil(log(1 - 0.99) / log(1 - e^M));
    % number of maximum inliers
    max_inliers = 0;
    best_inliers_idx = [];
    n = size(optPos, 1);

    for i = 1:k
        % select random indices
        idx = randperm(n);

        flow = optV(idx(1:M), :);
        x = optPos(idx(1:M), 1);
        y = optPos(idx(1:M), 2);
        depth = Z(idx(1:M), :);

        % construct H matrix
        H = zeros(2 * M, 6);

        for j = 1:M
            A = (1 / depth(j)) * [-1, 0, x(j); 0, -1, y(j)];
            B = [x(j) * y(j), -(1 + x(j)^2), y(j); 1 + y(j)^2, -x(j) * y(j), -x(j)];
            H(2 * j - 1:2 * j, :) = [A, B];
        end

        % find velocities
        currFlow = [flow(:, 1); flow(:, 2)];
        estVel = pinv(H) * currFlow;

        % find inliers
        x = optPos(:, 1);
        y = optPos(:, 2);
        num_inliers = 0;
        inliers_idx = [];

        for l = 1:n
            A = (1 / Z(l)) * [-1, 0, x(l); 0, -1, y(l)];
            B = [x(l) * y(l), -(1 + x(l)^2), y(l); 1 + y(l)^2, -x(l) * y(l), -x(l)];
            estPdot = [A, B] * estVel; % % optical flow prediction using the predicted velocities

            % comparing prection with true values of optical flow
            if norm(estPdot - optV(l, :)) < e
                num_inliers = num_inliers + 1;
                inliers_idx = [inliers_idx; l];
            end

        end

        % get maximum number of inliers
        if num_inliers > max_inliers
            max_inliers = num_inliers;
            best_inliers_idx = inliers_idx;
        end

    end

    % get best estimate of velocity from inliers
    idx = best_inliers_idx;
    idx_size = size(idx, 1);
    flow = optV(idx, :);
    x = optPos(idx, 1);
    y = optPos(idx, 2);
    depth = Z(idx, :);
    H = zeros(2 * idx_size, 6);

    for i = 1:idx_size
        A = (1 / depth(i)) * [-1, 0, x(i); 0, -1, y(i)];
        B = [x(i) * y(i), -(1 + x(i)^2), y(i); 1 + y(i)^2, -x(i) * y(i), -x(i)];
        H(2 * i - 1:2 * i, :) = [A, B];
    end

    currFlow = [flow(:, 1); flow(:, 2)];
    Vel = pinv(H) * currFlow;

end
