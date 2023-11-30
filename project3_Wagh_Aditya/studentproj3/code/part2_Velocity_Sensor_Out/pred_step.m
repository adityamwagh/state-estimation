function [covarEst, uEst] = pred_step(uPrev, covarPrev, angVel, acc, dt)
    %% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step

    %% Parameter Definition
    % uPrev - is the mean of the prev state
    % covarPrev - covar of the prev state
    % angVel - angular velocity input at the time step
    % acc - acceleration at the timestep
    % dt - difference in time

    %% Hyperparameters
    % defining the covariance matrix of the noise
    Q = eye(12) * 0.0001;

    % define noise
    noise = mvnrnd(zeros(12, 1), Q);
    ng = noise(1:3)';
    na = noise(4:6)';
    nbg = noise(7:9)';
    nba = noise(10:12)';
    alpha = 0.001;
    beta = 2;
    k = 1;

    %% Compute augmented state
    uPrev_aug = [uPrev; zeros(12, 1)];
    covarPrev_aug = [covarPrev, zeros(15, 12); zeros(12, 15), Q];

    % finding the sqrt of the covariance which is used for finding sigma points
    sqrt_cov = chol(covarPrev_aug);
    sqrt_cov = sqrt_cov';

    %% Compute sigma points
    n = size(covarPrev_aug, 1);
    lambda = double(alpha * alpha * (n + k) - n);
    sigma_pts = zeros(n, 2 * n + 1);
    sigma_pts(:, 1) = uPrev_aug;

    for i = 1:n
        sigma_pts(:, i + 1) = uPrev_aug + sqrt(n + lambda) * sqrt_cov(:, i);
        sigma_pts(:, n + i + 1) = uPrev_aug - sqrt(n + lambda) * sqrt_cov(:, i);
    end

    Wm = zeros(2 * n + 1, 1);
    Wm(1) = lambda / (n + lambda);
    Wm(2:2 * n + 1) = 1 / (2 * (n + lambda));
    Wc = zeros(2 * n + 1, 1);
    Wc(1) = (lambda / (n + lambda)) + (1 - alpha * alpha + beta);
    Wc(2:2 * n + 1) = 1 / (2 * (n + lambda));
    g = [0; 0; -9.8];

    curr_state = zeros(15, 2 * n + 1);
    curr_sigma_pts = zeros(n, 2 * n + 1);
    uEst = zeros(15, 1);
    covarEst = zeros(15, 15);

    for i = 1:(2 * n + 1)

        % Get angles from sigma point
        roll = sigma_pts(4, i);
        pitch = sigma_pts(5, i);
        yaw = sigma_pts(6, i);

        % Compute G and R
        G = [
            cos(pitch) * cos(yaw), -sin(yaw), 0;
            sin(yaw) * cos(pitch), cos(yaw), 0;
            -sin(pitch), 0, 1;
            ];
        R = eul2rotm([roll, pitch, yaw], 'XYZ');

        %% Propagate the sigma points through the nonlinear process model
        curr_state(:, i) = [
                        sigma_pts(7:9, i);
                        G \ R * (angVel - sigma_pts(10:12, i) - ng);
                        g + R * (acc - sigma_pts(13:15, i) - na);
                        nbg;
                        nba
                        ];

        % discretisation by euler approximation
        curr_sigma_pts(1:15, i) = sigma_pts(1:15, i) + curr_state(:, i) * dt;

        %% Compute mean and covariance matrices
        uEst = uEst + Wm(i, 1) .* curr_sigma_pts(1:15, i);
        covarEst = covarEst + Wc(i, 1) .* ((curr_sigma_pts(1:15, i) - uPrev) * (curr_sigma_pts(1:15, i) - uPrev)');
    end

end
