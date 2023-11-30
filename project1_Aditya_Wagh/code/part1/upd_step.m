function [uCurr, covar_curr] = upd_step(z_t, covarEst, uEst)
    % z_t is the measurement
    % covarEst and uEst are the predicted covariance and mean respectively
    % uCurr and covar_curr are the updated mean and covariance respectively
    C_t = [
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        ];

    W_t = eye(6, 6);

    mu = zeros(6, 1);
    % R = 0.0001 * eye(6);
    R = 0.0001 * diag([0.9, 0.5, 0.8, 0.2, 0.5, 0.9]);
    v = mvnrnd(mu, R);
    v = transpose(v);

    % since W_t = I, we will partially precompute the Kalman Gain
    K_t = covarEst * transpose(C_t) / (C_t * covarEst * transpose(C_t) + W_t * R * transpose(W_t));

    % mean and covariance
    uCurr = uEst + K_t * (z_t - (C_t * uEst + W_t * v));
    covar_curr = covarEst - K_t * C_t * covarEst;
end
