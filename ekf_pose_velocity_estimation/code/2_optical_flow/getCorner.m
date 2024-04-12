function res = getCorner(id)
    %% CHANGE THE NAME OF THE FUNCTION TO getCorner
    %% Input Parameter Description
    % id = List of all the AprilTag ids detected in the current image(data)

    %% Output Parameter Description
    % res = List of the coordinates of the 4 corners (or 4 corners and the
    % centre) of each detected AprilTag in the image in a systematic method
    % store all points

    res.p0 = zeros(2, length(id));
    res.p1 = zeros(2, length(id));
    res.p2 = zeros(2, length(id));
    res.p3 = zeros(2, length(id));
    res.p4 = zeros(2, length(id));

    for i = 1:length(id)

        % start at the first tag
        corners = [[0.076, 0.076]; [0.152, 0]; [0.152, 0.152]; [0, 0.152]; [0, 0]; ];

        % get rows and columns
        j = rem(id(i), 12);
        k = floor(id(i) / 12);

        % compute x coordinate of corners
        for a = 1:j
            corners(:, 1) = corners(:, 1) + (0.152 + 0.152);
        end

        % compute y coordinate of corners
        for b = 1:k

            if rem(b, 3) == 0
                corners(:, 2) = corners(:, 2) + (0.152 + 0.178);
            else
                corners(:, 2) = corners(:, 2) + (0.152 + 0.152);
            end

        end

        % interchange rows and columns of corners
        corners = transpose(corners);

        % store x values
        res.p0(1, i) = corners(1, 1);
        res.p1(1, i) = corners(1, 2);
        res.p2(1, i) = corners(1, 3);
        res.p3(1, i) = corners(1, 4);
        res.p4(1, i) = corners(1, 5);

        % store y values
        res.p0(2, i) = corners(2, 1);
        res.p1(2, i) = corners(2, 2);
        res.p2(2, i) = corners(2, 3);
        res.p3(2, i) = corners(2, 4);
        res.p4(2, i) = corners(2, 5);

    end
