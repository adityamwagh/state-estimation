function [position, orientation, R_c2w] = estimatePose(data, t)
    %% CHANGE THE NAME OF THE FUNCTION TO estimatePose
    % Please note that the coordinates for each corner of each AprilTag are
    % defined in the world frame, as per the information provided in the
    % handout. Ideally a call to the function getCorner with ids of all the
    % detected AprilTags should be made. This function should return the X and
    % Y coordinate of each corner, or each corner and the centre, of all the
    % detected AprilTags in the image. You can implement that anyway you want
    % as long as the correct output is received. A call to that function
    % should made from this function.
    %% Input Parameter Definition
    % data = the entire data loaded in the current dataset
    % t = index of the current data in the dataset

    %% Output Parameter Definition

    % position = translation vector representing the position of the
    % drone(body) in the world frame in the current time, in the order XYZ

    % orientation = euler angles representing the orientation of the
    % drone(body) in the world frame in the current time, in the order XYZ

    %R_c2w = Rotation which defines camera to world frame

    A = [];

    % calculate world coordinates of centers & corners of a tag
    corners = getCorner(data(t).id);
    %
    for i = 1:(length(data(t).id))

        % compute A matrix for each
        temp = [
            [corners(1).p0(1, i), corners(1).p0(2, i), 1, 0, 0, 0, -data(t).p0(1, i) * corners(1).p0(1, i), -data(t).p0(1, i) * corners(1).p0(2, i), -data(t).p0(1, i)];
            [0, 0, 0, corners(1).p0(1, i), corners(1).p0(2, i), 1, -data(t).p0(2, i) * corners(1).p0(1, i), -data(t).p0(2, i) * corners(1).p0(2, i), -data(t).p0(2, i)];
            [corners(1).p1(1, i), corners(1).p1(2, i), 1, 0, 0, 0, -data(t).p1(1, i) * corners(1).p1(1, i), -data(t).p1(1, i) * corners(1).p1(2, i), -data(t).p1(1, i)];
            [0, 0, 0, corners(1).p1(1, i), corners(1).p1(2, i), 1, -data(t).p1(2, i) * corners(1).p1(1, i), -data(t).p1(2, i) * corners(1).p1(2, i), -data(t).p1(2, i)];
            [corners(1).p2(1, i), corners(1).p2(2, i), 1, 0, 0, 0, -data(t).p2(1, i) * corners(1).p2(1, i), -data(t).p2(1, i) * corners(1).p2(2, i), -data(t).p2(1, i)];
            [0, 0, 0, corners(1).p2(1, i), corners(1).p2(2, i), 1, -data(t).p2(2, i) * corners(1).p2(1, i), -data(t).p2(2, i) * corners(1).p2(2, i), -data(t).p2(2, i)];
            [corners(1).p3(1, i), corners(1).p3(2, i), 1, 0, 0, 0, -data(t).p3(1, i) * corners(1).p3(1, i), -data(t).p3(1, i) * corners(1).p3(2, i), -data(t).p3(1, i)];
            [0, 0, 0, corners(1).p3(1, i), corners(1).p3(2, i), 1, -data(t).p3(2, i) * corners(1).p3(1, i), -data(t).p3(2, i) * corners(1).p3(2, i), -data(t).p3(2, i)];
            [corners(1).p4(1, i), corners(1).p4(2, i), 1, 0, 0, 0, -data(t).p4(1, i) * corners(1).p4(1, i), -data(t).p4(1, i) * corners(1).p4(2, i), -data(t).p4(1, i)];
            [0, 0, 0, corners(1).p4(1, i), corners(1).p4(2, i), 1, -data(t).p4(2, i) * corners(1).p4(1, i), -data(t).p4(2, i) * corners(1).p4(2, i), -data(t).p4(2, i)];
            ];

        A = cat(1, A, temp);

    end

    % compute svd of A
    [~, ~, Vh] = svd(A);
    Vh = Vh * sign(Vh(9, 9));

    % get h from v
    h = Vh(:, 9);
    h = [h(1), h(2), h(3); h(4), h(5), h(6); h(7), h(8), h(9)];

    % calibration matrix from parameters.txt
    K = [311.0520, 0, 201.8724; 0, 311.3885, 113.6210; 0, 0, 1];

    Rt = K \ h;
    Rt = Rt * sign(Rt(3, 3));

    % translation vector in the camera frame
    T_c2w = Rt(:, 3) / norm(Rt(:, 1));
    newRt = [Rt(:, 1), Rt(:, 2), cross(Rt(:, 1), Rt(:, 2))];

    % minimize R
    [Ur, ~, Vr] = svd(newRt);
    Sr = [1, 0, 0; 0, 1, 0; 0, 0, det(Ur * Vr)];
    R_c2w = Ur * Sr * Vr; % Rotation in the camera frame

    % change from camera to body
    H_c2w = [R_c2w', -R_c2w' * T_c2w; 0, 0, 0, 1];
    H_b2c = [eul2rotm([-pi / 4, 0, -pi]), [-0.04; 0; -0.03]; 0, 0, 0, 1];
    H_b2w = H_c2w * H_b2c;
    R_b2w = H_b2w(1:3, 1:3);
    T_b2w = H_b2w(1:3, 4);

    orientation = (rotm2eul(R_b2w));
    position = (T_b2w);
end
