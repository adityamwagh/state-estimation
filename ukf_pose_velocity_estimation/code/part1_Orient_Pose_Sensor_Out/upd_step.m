function [uCurr, covar_curr] = upd_step(z_t, covarEst, uEst)
    %% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    % z_t - is the sensor data at the time step
    % covarEst - estimated covar of the  state
    % uEst - estimated mean of the state
    C = [eye(3) zeros(3) zeros(3) zeros(3) zeros(3); ...
        zeros(3) eye(3) zeros(3) zeros(3) zeros(3)];

    %% Define the covariance matrix of the noise
    R = eye(6) * 0.01;

    %% Calculate the kalman gain
    K = (covarEst * C') / (C * covarEst * C' + R);

    %% Update the mean and covariance
    uCurr = uEst + K * (z_t - C * uEst);
    covar_curr = covarEst - K * C * covarEst;

end
