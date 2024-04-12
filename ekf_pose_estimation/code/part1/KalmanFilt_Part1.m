clear; % Clear variables
datasetNum = 9; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime] = init(datasetNum);
Z = sampledVicon(1:6, :); % all the measurements that you need for the update

% Set initial condition
uPrev = vertcat(sampledVicon(1:9, 1), zeros(6, 1)); % Copy the Vicon Initial state
covarPrev = eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); % Just for saving state his.
prevTime = 0; % Last time step in real time

% Write your code here calling the pred_step.m and upd_step.m functions
for i = 1:length(sampledTime)

    % control inputs
    angVel = sampledData(i).omg;
    acc = sampledData(i).acc;

    % time step computations
    dt = sampledTime(i) - prevTime;
    prevTime = sampledTime(i);

    % mean and covariance estimation
    [covarEst, uEst] = pred_step(uPrev, covarPrev, angVel, acc, dt);
    [uCurr, covar_curr] = upd_step(Z(:, i), covarEst, uEst);

    % save currrent estimated state
    savedStates(:, i) = uCurr;

    % save current step to be used as next state
    uPrev = uCurr;
    covarPrev = covar_curr;

end

plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);
