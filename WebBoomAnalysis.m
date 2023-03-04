clear; close all; clc;
% This code performs a web boom analysis.
% It assumes that webs are always straight between two booms.

% Axis follow Right hand rule
% Positive x is from wing leading edge to trailing edge
% Positive y is up
% Positive z is along left wing

%
% Defining Loads
%

% The foce arrays are broken down as:
% 1st Element: Force Value
% 2nd Element: Force x position
% 3rd Element: Force y position
Fx = [10000, -0.623, 40.6923];
Fy = [1044440, -322.6923, 104.6923];
Fz = [1, 0, 0];

Mx = -2800;
My = -5000;
Mz = 2000;


%
% Web & Boom Parameter Inputs
%

% Boom Positions
x_b = [-1, 1, 1, -1];
y_b = [1, 1, -1, -1];


% Boom Youngs Modulus
E = 1e9*[70, 70, 70, 70];

% Boom area
A = 1e-3*[10, 12, 1, 1];

% Web Shear Modulus
G = 1e9*[28, 28, 28, 28];

% Web Thickness
t = [0.002, 0.002, 0.002, 0.002];

%
% Size of Array
%

[i_h, i_w] = size(x_b);


%
% Finding the Stiffness Centre
%

x_c = sum(E.*A.*x_b) / sum(E.*A);
y_c = sum(E.*A.*y_b) / sum(E.*A);

%
% Converting Loads into Stiffeness Centre Axis
%

Mx_c = Mx + Fz(1)*Fz(3);
My_c = My - Fz(1)*(Fz(2) - x_c);
Mz_c = Mz + Fx(1)*(Fx(3) - y_c) - Fy(1)*(Fy(2) - x_c);

Fx_c = Fx(1);
Fy_c = Fy(1);
Fz_c = Fz(1);

%
% Moving Coordinate System to Stiffeness Centre as Origin
%

x_b_c = x_b - x_c*ones(size(x_b));
y_b_c = y_b - y_c*ones(size(y_b));

%
% Boom Stiffeness Coefficients
%

F = sum(E.*A);
Kxy = sum(E.*A.*x_b_c.*y_b_c);
Kxx = sum(E.*A.*y_b_c.^2);
Kyy = sum(E.*A.*x_b_c.^2);


%
% Boom Stress Calcaultion
%

stress_b = -(My_c*Kxx + Mx_c*Kxy)*(E.*x_b_c)/(Kxx*Kyy - Kxy^2) + (Mx_c*Kyy + My_c*Kxy)*(E.*y_b_c)/(Kxx*Kyy - Kxy^2) + Fz_c/F*E;

%
% Web Stiffeness Coefficients
%

Rx = zeros(1, i_w);
Ry = zeros(1, i_w);
for i = 1:i_w
    Rx(i) = sum(E(1:i).*A(1:i).*y_b_c(1:i));
    Ry(i) = sum(E(1:i).*A(1:i).*x_b_c(1:i));
end

%
% Open Shear Flow Calculation
%

q_open = +Fx_c*(Kxx/(Kxx*Kyy - Kxy^2))*Ry - Fy_c*(Kxy/(Kxx*Kyy - Kxy^2))*Ry ...
         +Fy_c*(Kyy/(Kxx*Kyy - Kxy^2))*Rx - Fx_c*(Kxy/(Kxx*Kyy - Kxy^2))*Rx;

%
% Finding Open Shear Centre of Structure
%

% Finds the horizontal and vertical components of the force due to the
% shear flow in the webs. Then uses the respective x and y distance to find
% the moment caused by each of the shear flows.
%
% By finding the vertical and horizontal forces you can verify that the
% code works as the internal forces should equal the external forces.
s = zeros(1, i_w);
for i = 1:i_w-1
    % Finds the angle from the x-axis to the panel. Assumes that shear flow
    % acts along the web
    theta(i) = atan2(y_b_c(i+1) - y_b_c(i), x_b_c(i+1) - x_b_c(i));
    
    % Finds length of the Web. Note that the for loop doesn't calculate the
    % last web until outside the for loop.
    s(i) = sqrt((x_b_c(i+1)-x_b_c(i))^2 + (y_b_c(i+1)-y_b_c(i))^2);

    % This is because the last web is connected to the last and first boom.
end

% Finds the length of the final web.
s(i_w) = sqrt((x_b_c(1)-x_b_c(i_w))^2 + (y_b_c(1)-y_b_c(i_w))^2);
theta(i_w) = atan2(y_b_c(1) - y_b_c(i_w), x_b_c(1) - x_b_c(i_w));

M_q_open_x = q_open.*s.*cos(theta).*y_b_c;
M_q_open_y = q_open.*s.*sin(theta).*x_b_c;
M_q_open = - M_q_open_x + M_q_open_y;

% Finds the total moment caused by the shear flows
M_q_open_tot = sum(M_q_open);



%
% Finding Shear Centre
%

shear_centre_open = - (M_q_open_tot)/sqrt(Fy_c^2 + Fx_c^2);



%
% Finding the shear flow of the closed section
%

% Finding area of the closed cell
area = polyarea(x_b, y_b);
omega = 2*area;

% Finding the torque experienced by the structure around the shear
% centre caused by the force loads.
To = shear_centre_open*sqrt(Fy_c^2 + Fx_c^2);

% Calcaulting the closed cell shear flow due to 
qo = -(To + Mz_c) / omega;

%
% Calcaulting the total Shear flow
%

% Find the total shear flow, including effects due to torques.
q = q_open + qo;

%
% Finding anlge of twist
%

Twist = (1/omega)*sum(q.*s./G./t) * 180/pi;

%
% Finding Shear Centre of the Web-Boom Strucutre
%

% It is assumed that the distance between the origin and the shear centre
% is perpendicular to the resultant force vector. By finding the angle from
% the positive x axis between the x and y force, you can find the angle
% between the shear centre line and the vertical y axis.
%
% Using the angle and length of the shear centre, an x and y coordinate can
% be found.

%shear_centre_closed = -(omega * sum(q./G./t.*s) / sum(s./G./t)) / sqrt(Fx_c^2 + Fy_c^2);
x_shear_centre = (-(omega * sum(q./G./t.*s) / sum(s./G./t) + Mz_c) / Fy_c + x_c) * Fy_c/Fx_c;
y_shear_centre = (+(omega * sum(q./G./t.*s) / sum(s./G./t) + Mz_c) / Fx_c + y_c) * Fx_c/Fy_c;




X_eq = sum(q.*s.*cos(theta)) + Fx_c;
Y_eq = sum(q.*s.*sin(theta)) + Fy_c;
Z_eq = sum(stress_b.*A) - Fz_c;
M_eq = sum(q.*s.*cos(theta).*y_b_c) - sum(q.*s.*sin(theta).*x_b_c) + Mz_c;

fprintf('Boom Stress (MPa): \n\t')
disp(stress_b*1e-6)

fprintf('Web Shear Flow (M N/m)): \n\t')
disp(q*1e-6)

fprintf('Twist (deg): \n\t')
disp(Twist)

fprintf('x shear centre (m): \n\t')
disp(x_shear_centre)

fprintf('y shear centre (m): \n\t')
disp(y_shear_centre)


