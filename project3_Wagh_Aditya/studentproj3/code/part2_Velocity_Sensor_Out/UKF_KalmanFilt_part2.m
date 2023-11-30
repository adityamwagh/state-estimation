clear; % Clear variables
addpath('../data')
datasetNum = 4; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime, proj2Data] = init(datasetNum);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9, 1), zeros(6, 1)); % Copy the Vicon Initial state
covarPrev = 0.01 * eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0;
vel = proj2Data.linearVel;
angVel2 = proj2Data.angVel;

%% Calculate Kalman Filter
for i = 1:length(sampledTime)

    % control inputs
    angVel = sampledData(i).omg;
    acc = sampledData(i).acc;

    % time step computations
    dt = sampledTime(i) - prevTime;
    prevTime = sampledTime(i);

    % mean and covariance estimation
    [covarEst, uEst] = pred_step(uPrev, covarPrev, angVel, acc, dt);
    [uCurr, covar_curr] = upd_step(transpose([vel(i, :), angVel2(i, :)]), covarEst, uEst);

    % save currrent estimated state
    savedStates(:, i) = uCurr;

    % save current step to be used as next state
    uPrev = uCurr;
    covarPrev = covar_curr;
end

plotData(savedStates, sampledTime, sampledVicon, 2, datasetNum);
