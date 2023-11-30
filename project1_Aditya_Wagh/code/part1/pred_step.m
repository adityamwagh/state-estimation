function [covarEst, uEst] = pred_step(uPrev, covarPrev, angVel, acc, dt)
    % covarPrev and uPrev are the previous mean and covariance respectively
    % angVel is the angular velocity
    % acc is the acceleration
    % dt is the sampling time

    % fetch values from previous state
    px = uPrev(1);
    py = uPrev(2);
    pz = uPrev(3);
    x = uPrev(4);
    y = uPrev(5);
    z = uPrev(6);
    vx = uPrev(7);
    vy = uPrev(8);
    vz = uPrev(9);
    bgx = uPrev(10);
    bgy = uPrev(11);
    bgz = uPrev(12);
    bax = uPrev(13);
    bay = uPrev(14);
    baz = uPrev(15);

    % control inputs from gyroscope and accelerometer
    wmx = angVel(1);
    wmy = angVel(2);
    wmz = angVel(3);
    amx = acc(1);
    amy = acc(2);
    amz = acc(3);

    % make some noise
    mu = zeros(12, 1);
    Q = 0.00015 * eye(12);
    Q_d = Q * dt;
    n = mvnrnd(mu, Q_d);
    n = transpose(n);

    % nullify noise to compute state matrices, process model and predicted state
    ngx = n(1);
    ngy = n(2);
    ngz = n(3);
    nax = n(4);
    nay = n(5);
    naz = n(6);
    nbgx = n(7);
    nbgy = n(8);
    nbgz = n(9);
    nbax = n(10);
    nbay = n(11);
    nbaz = n(12);

    A_t = [
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
        0, 0, 0, - (sin(y) * (bgy * cos(x) + ngy * cos(x) - bgz * sin(x) - wmy * cos(x) - ngz * sin(x) + wmz * sin(x))) / cos(y), - (bgz * cos(x) + ngz * cos(x) + bgy * sin(x) - wmz * cos(x) + ngy * sin(x) - wmy * sin(x)) / cos(y)^2, 0, 0, 0, 0, -1, - (sin(x) * sin(y)) / cos(y), - (cos(x) * sin(y)) / cos(y), 0, 0, 0;
        0, 0, 0, bgz * cos(x) + ngz * cos(x) + bgy * sin(x) - wmz * cos(x) + ngy * sin(x) - wmy * sin(x), 0, 0, 0, 0, 0, 0, -cos(x), sin(x), 0, 0, 0;
        0, 0, 0, - (bgy * cos(x) + ngy * cos(x) - bgz * sin(x) - wmy * cos(x) - ngz * sin(x) + wmz * sin(x)) / cos(y), - (sin(y) * (bgz * cos(x) + ngz * cos(x) + bgy * sin(x) - wmz * cos(x) + ngy * sin(x) - wmy * sin(x))) / cos(y)^2, 0, 0, 0, 0, 0, -sin(x) / cos(y), -cos(x) / cos(y), 0, 0, 0;
        0, 0, 0, - (sin(x) * sin(z) + cos(x) * cos(z) * sin(y)) * (bay - amy + nay) - (cos(x) * sin(z) - cos(z) * sin(x) * sin(y)) * (baz - amz + naz), cos(z) * sin(y) * (bax - amx + nax) - cos(x) * cos(y) * cos(z) * (baz - amz + naz) - cos(y) * cos(z) * sin(x) * (bay - amy + nay), (cos(x) * cos(z) + sin(x) * sin(y) * sin(z)) * (bay - amy + nay) - (cos(z) * sin(x) - cos(x) * sin(y) * sin(z)) * (baz - amz + naz) + cos(y) * sin(z) * (bax - amx + nax), 0, 0, 0, 0, 0, 0, -cos(y) * cos(z), cos(x) * sin(z) - cos(z) * sin(x) * sin(y), - sin(x) * sin(z) - cos(x) * cos(z) * sin(y);
        0, 0, 0, (cos(z) * sin(x) - cos(x) * sin(y) * sin(z)) * (bay - amy + nay) + (cos(x) * cos(z) + sin(x) * sin(y) * sin(z)) * (baz - amz + naz), sin(y) * sin(z) * (bax - amx + nax) - cos(y) * sin(x) * sin(z) * (bay - amy + nay) - cos(x) * cos(y) * sin(z) * (baz - amz + naz), (cos(x) * sin(z) - cos(z) * sin(x) * sin(y)) * (bay - amy + nay) - (sin(x) * sin(z) + cos(x) * cos(z) * sin(y)) * (baz - amz + naz) - cos(y) * cos(z) * (bax - amx + nax), 0, 0, 0, 0, 0, 0, -cos(y) * sin(z), - cos(x) * cos(z) - sin(x) * sin(y) * sin(z), cos(z) * sin(x) - cos(x) * sin(y) * sin(z);
        0, 0, 0, cos(y) * sin(x) * (baz - amz + naz) - cos(x) * cos(y) * (bay - amy + nay), cos(y) * (bax - amx + nax) + cos(x) * sin(y) * (baz - amz + naz) + sin(x) * sin(y) * (bay - amy + nay), 0, 0, 0, 0, 0, 0, 0, sin(y), -cos(y) * sin(x), -cos(x) * cos(y);
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

        ];

    U_t = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        - 1, - (sin(x) * sin(y)) / cos(y), - (cos(x) * sin(y)) / cos(y), 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, -cos(x), sin(x), 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, -sin(x) / cos(y), -cos(x) / cos(y), 0, 0, 0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, -cos(y) * cos(z), cos(x) * sin(z) - cos(z) * sin(x) * sin(y), - sin(x) * sin(z) - cos(x) * cos(z) * sin(y), 0, 0, 0, 0, 0, 0;
        0, 0, 0, -cos(y) * sin(z), - cos(x) * cos(z) - sin(x) * sin(y) * sin(z), cos(z) * sin(x) - cos(x) * sin(y) * sin(z), 0, 0, 0, 0, 0, 0;
        0, 0, 0, sin(y), -cos(y) * sin(x), -cos(x) * cos(y), 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

        ];

    % f(mu_(t-1), u_t, 0)
    f = [
        vx;
        vy;
        vz;
        - (bgx * cos(y) + ngx * cos(y) - wmx * cos(y) + bgz * cos(x) * sin(y) + ngz * cos(x) * sin(y) + bgy * sin(x) * sin(y) - wmz * cos(x) * sin(y) + ngy * sin(x) * sin(y) - wmy * sin(x) * sin(y)) / cos(y);
        bgz * sin(x) - ngy * cos(x) - bgy * cos(x) + wmy * cos(x) + ngz * sin(x) - wmz * sin(x);
        - (bgz * cos(x) + ngz * cos(x) + bgy * sin(x) - wmz * cos(x) + ngy * sin(x) - wmy * sin(x)) / cos(y);
        (cos(x) * sin(z) - cos(z) * sin(x) * sin(y)) * (bay - amy + nay) - (sin(x) * sin(z) + cos(x) * cos(z) * sin(y)) * (baz - amz + naz) - cos(y) * cos(z) * (bax - amx + nax);
        (cos(z) * sin(x) - cos(x) * sin(y) * sin(z)) * (baz - amz + naz) - (cos(x) * cos(z) + sin(x) * sin(y) * sin(z)) * (bay - amy + nay) - cos(y) * sin(z) * (bax - amx + nax);
        sin(y) * (bax - amx + nax) - cos(x) * cos(y) * (baz - amz + naz) - cos(y) * sin(x) * (bay - amy + nay) - 981/100;
        nbgx;
        nbgy;
        nbgz;
        nbax;
        nbay;
        nbaz;

        ];

    % state matrices
    F_t = eye(15) + dt * A_t;
    V_t = U_t;

    % Final state prediction equations
    uEst = uPrev + dt * f; % can replace * with matrix multiplication; + with matrix addition
    covarEst = F_t * covarPrev * transpose(F_t) + V_t * Q_d * transpose(V_t);
end
