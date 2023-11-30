function [uCurr, covar_curr] = upd_step(z_t, covarEst, uEst)

    %% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    % defining the covariance matrix of the noise
    R = eye(3) * 0.0002;

    % finding the sqrt of the covariance which is used for finding sigma points
    sqrt_cov = chol(covarEst);
    sqrt_cov = sqrt_cov';

    %% Compute the sigma points
    % hyperparameters
    alpha = 0.0001;
    beta = 2;
    k = 1;
    n = size(covarEst, 1);
    lambda = double(alpha * alpha * (n + k) - n);

    % compute sigma points
    sigma_pts = zeros(n, 2 * n + 1);
    sigma_pts(:, 1) = uEst;

    for i = 1:n
        sigma_pts(:, i + 1) = uEst + sqrt(n + lambda) * sqrt_cov(:, i);
        sigma_pts(:, n + i + 1) = uEst - sqrt(n + lambda) * sqrt_cov(:, i);
    end

    %% Compute weights to find the new mean and covariance
    % compute Wm
    Wm = zeros(2 * n + 1, 1);
    Wm(1) = lambda / (n + lambda);
    Wm(2:2 * n + 1) = 1 / (2 * (n + lambda));

    % Compute Wc
    Wc = zeros(2 * n + 1, 1);
    Wc(1) = lambda / (n + lambda) + (1 - alpha * alpha + beta);
    Wc(2:2 * n + 1) = 1 / (2 * (n + lambda));

    %% Propagate the sigma points through the non-linear function
    velc2w_c = zeros(6, 2 * n + 1);
    z_ut = zeros(3, 1);
    S = zeros(3, 3);
    C = zeros(n, 3);

    %% defining rotations and translations required to chnage the frame of the velocities in the loop
    % rotation from body to camera
    R_b2c = eul2rotm([-pi / 4, 0, -pi], 'ZYX');

    % rotation from camera to body
    R_c2b = R_b2c';

    % translation from body to camera
    Tb2c_b = [0.04 * cos(pi / 4); -0.04 * sin(pi / 4); -0.03];

    for i = 1:2 * n + 1

        % constructing the adjoint to convert the velocities from the body frame to the camera frame
        % Get angles from sigma point
        roll = sigma_pts(4, i);
        pitch = sigma_pts(5, i);
        yaw = sigma_pts(6, i);

        R_b2w = eul2rotm([yaw, pitch, roll], 'ZYX');

        % body velocity from world frame to body frame
        v_b2w_w = sigma_pts(7:9, i);
        v_b2w_b = R_b2w' * v_b2w_w;

        % compute skewsymmetrix form of translation
        skew_Tb2c_b = [0, -Tb2c_b(3), Tb2c_b(2); Tb2c_b(3), 0, -Tb2c_b(1); -Tb2c_b(2), Tb2c_b(1), 0];

        % propogating the sigma points using the adjoint nonlinear equation
        velc2w_c(:, i) = [R_b2c, -R_b2c * skew_Tb2c_b; zeros(3), R_b2c] * [v_b2w_b; R_c2b * z_t(4:6)];

        % Calculating z_ut using only linear velocity
        z_ut = z_ut + velc2w_c(1:3, i) * Wm(i);
    end

    % Calculating the covariance
    for i = 1:2 * n + 1
        C = C + Wc(i) * (sigma_pts(:, i) - uEst) * (velc2w_c(1:3, i) - z_ut)';
        S = S + Wc(i) * (velc2w_c(1:3, i) - z_ut) * (velc2w_c(1:3, i) - z_ut)' + R;

    end

    % updating the mean and variance of the state
    K = C / S;
    uCurr = uEst + K * (z_t(1:3) - z_ut);
    covar_curr = covarEst - K * S * K';

end
