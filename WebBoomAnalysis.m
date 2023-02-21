close all; clear; clc;
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
Fx = [2000, 0.0925, 0.7767];
Fy = [0, 0, 0];
Fz = [0, 0, 0];

Mx = 0;
My = 0;
Mz = -138;


%
% Web & Boom Parameter Inputs
%

% Boom Positions
x_b = [-0.21, 0.31, 0.41, -0.14];
y_b = [0.1, 2.1, -0.31, -0.14];


% Boom Youngs Modulus
E = 1e9*[200, 200, 200, 200];

% Boom area
A = [0.5, 0.5, 0.5, 0.5];

% Web Shear Modulus
G = 1e6*[20, 20, 20, 20];

% Web Thickness
t = [0.15, 0.15, 0.15, 0.15];

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
My_c = My - Fz(1)*Fz(2);
Mz_c = Mz - Fx(1)*(Fx(3) - y_c) + Fy(1)*(Fy(2) - x_c);

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

shear_centre_open = - M_q_open_tot/sqrt(Fy_c^2 + Fx_c^2);


%
% Finding the shear flow of the closed section
%

% Finding area of the closed cell
area = polyarea(x_b, y_b);
omega = 2*area;

% Finding the torque experienced by the structure around the shear
% centre caused by the force loads.
To = shear_centre_open*sqrt(Fy_c^2 + Fx_c^2) + Mz_c;

% Calcaulting the closed cell shear flow due to 
qo = -To / omega;

%
% Calcaulting the total Shear flow
%

% Find the total shear flow, including effects due to torques.
q = q_open + qo;

%
% Finding anlge of twist
%

Twist = (1/omega)*sum(q.*s./G./t) * 180/pi

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

shear_centre_closed = -(omega * sum(q./G./t.*s) / sum(s./G./t) + Mz_c) / sqrt(Fx_c^2 + Fy_c^2);
x_shear_centre = shear_centre_closed*sin(atan2(Fy_c , Fx_c)) + x_c
y_shear_centre = shear_centre_closed*cos(atan2(Fy_c , Fx_c)) + y_c

% x_shear_centre_closed = shear_center_distance*sin(force_angle)
% y_shear_centre_closed = shear_center_distance*cos(force_angle)


X_eq = sum(q.*s.*cos(theta)) + Fx_c
Y_eq = sum(q.*s.*sin(theta)) + Fy_c
M_eq = sum(q.*s.*cos(theta).*y_b_c) - sum(q.*s.*sin(theta).*x_b_c) + Mz_c


