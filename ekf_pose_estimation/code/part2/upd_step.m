function [uCurr, covar_curr] = upd_step(z_t, covarEst, uEst)
    % z_t is the measurement
    % covarEst and uEst are the predicted covariance and mean respectively
    % uCurr and covar_curr are the updated mean and covariance respectively
    C_t = [
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
        ];

    W_t = eye(3, 3);
    mu = zeros(3, 1);
    R = 0.0001 * diag([0.9, 0.5, 0.8]);
    v = mvnrnd(mu, R);
    v = transpose(v);

    % since W_t = I, we will partially precompute the Kalman Gain
    K_t = covarEst * transpose(C_t) / (C_t * covarEst * transpose(C_t) + W_t * R * transpose(W_t));

    % mean and covariance
    uCurr = uEst + K_t * (z_t - (C_t * uEst + W_t * v));
    covar_curr = covarEst - K_t * C_t * covarEst;
end

% R = [
%     cos(vcon_thetax) * cos(vcon_thetay), cos(vcon_thetax) * sin(vcon_thetay) * sin(vcon_thetaz) - sin(vcon_thetax) * cos(vcon_thetaz), cos(vcon_thetax) * sin(vcon_thetay) * cos(vcon_thetaz) + sin(vcon_thetax) * sin(vcon_thetaz);
%     sin(vcon_thetax) * cos(vcon_thetay), sin(vcon_thetax) * sin(vcon_thetay) * sin(vcon_thetaz) + cos(vcon_thetax) * cos(vcon_thetaz), sin(vcon_thetax) * sin(vcon_thetay) * cos(vcon_thetaz) - cos(vcon_thetax) * sin(vcon_thetaz);
%     -sin(vcon_thetay), cos(vcon_thetay) * sin(vcon_thetaz), cos(vcon_thetay) * cos(vcon_thetaz)
%     ];
