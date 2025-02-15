function z = elliptical_2D_gaussian_withd(xy,params)
    % Two-dimensional Gaussian function with an extra parameter d (to allow us to fit upside down gaussian like data)
    x = xy(:,1);            % x coordinates
    y = xy(:,2);            % y coordinates
    amp = params(1);        % amplitude of the Gaussian (will be added to d - the minimum)
    x0 = params(2);         % x localization of Gaussian peak
    y0 = params(3);         % y localization of Gaussian peak
    sigma_X = params(4);    % std along the main axis
    sigma_Y = params(5);    % std along the axis perpendicular to the main axis
    theta = params(6);      % main axis in radians
    d = params(7);          % minimum of the Gaussian

    a = cos(theta)^2 / (2 * sigma_X^2) + sin(theta)^2 / (2 * sigma_Y^2);
    b = sin(2 * theta) / (4 * sigma_X^2) - sin(2 * theta) / (4 * sigma_Y^2);
    c = sin(theta)^2 / (2 * sigma_X^2) + cos(theta)^2 / (2 * sigma_Y^2);

    first_piece=a.*(x-x0).^2;
    second_piece=2.*b.*(x-x0).*(y-y0);
    third_piece=c.*(y-y0).^2;
    
    z = amp.*exp(-(first_piece+second_piece+third_piece))+d;
end
