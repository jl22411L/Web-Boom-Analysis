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
Fx = [9080, 0, 0];
Fy = [10, 0, 0];
Fz = [0, 0, 0];

Mx = 0;
My = 0;
Mz = 5078;


%
% Web & Boom Parameter Inputs
%

% Boom Positions
x_b = [0, 0.2, 0.4, 0.6, 0.8, 1.4, 1.6];
y_b = [0, 0.3, 0.7, 1.2, -3 ,-5, -9];


% Boom Youngs Modulus
E = 1e9*[70, 70, 70, 70, 70, 70, 70];

% Boom area
A = 12e-3*[1, 1, 1, 1, 1, 1, 1];

% Web Shear Modulus
G = 1e9*[28, 28, 28, 28, 28, 28, 28];

% Web Thickness
t = [0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002];

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

% These are used for shear centre calcaultions
Mz_c_x = Mz + Fx(1)*(Fx(3) - y_c);
Mz_c_y = Mz - Fy(1)*(Fy(2) - x_c);

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

% Find the open cell shear flows due to x and y force. This is done for
% shear flow calculations.
qx_open = +Fx_c*(Kxx/(Kxx*Kyy - Kxy^2))*Ry - Fx_c*(Kxy/(Kxx*Kyy - Kxy^2))*Rx;
qy_open = +Fy_c*(Kyy/(Kxx*Kyy - Kxy^2))*Rx - Fy_c*(Kxy/(Kxx*Kyy - Kxy^2))*Ry;
q_open = qx_open + qy_open;
         

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

% This is finding the total moment
M_q_open_x = q_open.*s.*cos(theta).*y_b_c;
M_q_open_y = q_open.*s.*sin(theta).*x_b_c;
M_q_open = - M_q_open_x + M_q_open_y; 

% This is finding the moments from the open cell shear flow due to the x
% force.
M_qx_open_x = qx_open.*s.*cos(theta).*y_b_c;
M_qx_open_y = qx_open.*s.*sin(theta).*x_b_c;
M_qx_open = - M_qx_open_x + M_qx_open_y;
M_qx_open_tot = sum(M_qx_open);

% This is finding the moments from the open cell shear flow due to the x
% force.
M_qy_open_x = qy_open.*s.*cos(theta).*y_b_c;
M_qy_open_y = qy_open.*s.*sin(theta).*x_b_c;
M_qy_open = - M_qy_open_x + M_qy_open_y; 
M_qy_open_tot = sum(M_qy_open);

% Finds the total moment caused by the shear flows
M_q_open_tot = sum(M_q_open);
% REMEMBER THAT SHEAR FLOW IS OFF BACK. NEED TO BRING TO FRONT FOR
% EQUILIBRIUM CALCULATIONS


%
% Finding Shear Centre
%

% Bassically, the open cell shear centre os calculated using the shear flow
% due to x and y forces independently. Using the total open cell shear flow
% just caused to many problems and led to an indeterminate problem.
x_shear_centre_open = -M_qy_open_tot / Fy_c;
y_shear_centre_open = M_qx_open_tot / Fx_c;

%
% Finding the shear flow of the closed section
%

% Finding area of the closed cell
area = polyarea(x_b, y_b);
omega = 2*area;

% Finding the torque experienced by the structure around the shear
% centre caused by the force loads.

% Torque due to x force
To_y = x_shear_centre_open*Fy_c;
% Torque due to y force
To_x = - y_shear_centre_open*Fx_c;

% Calcaulting the closed cell shear flow due to their respective load cases
qo_x = (To_x - Mz_c_x) / omega;
qo_y = (To_y - Mz_c_y) / omega;
qo_tot = qo_x + qo_y;

%
% Calcaulting the total Shear flow
%

% Find the total shear flow, including effects due to torques.
qx = qx_open + qo_x;
qy = qy_open + qo_y;
q = q_open + qo_tot;

%
% Finding anlge of twist
%

Twist = (1/omega)*sum(q.*s./G./t) * 180/pi;

%
% Finding Shear Centre of the Web-Boom Strucutre
%

% So using the total shear flows for the individual forces, the shear
% centre is foind. Similar to the open cell shear flow calcaultions.

x_shear_centre = -(omega * sum(qy./G./t.*s) / sum(s./G./t) + Mz_c_y) / Fy_c + x_c;
y_shear_centre = (omega * sum(qx./G./t.*s) / sum(s./G./t) + Mz_c_x) / Fx_c + y_c;


%
% Equilibrium Calcualtions
%

% Checking X Force Equlibrium
% Xx_eq = sum(qx.*s.*cos(theta)) + Fx_c
% Yx_eq = sum(qx.*s.*sin(theta))
% Mx_eq = -sum(qx.*s.*cos(theta).*y_b_c) + sum(qx.*s.*sin(theta).*x_b_c) + Mz_c_x

% Checking Y Force Equlibrium
% Xy_eq = sum(qy.*s.*cos(theta))
% Yy_eq = sum(qy.*s.*sin(theta)) + Fy_c
% My_eq = -sum(qy.*s.*cos(theta).*y_b_c) + sum(qy.*s.*sin(theta).*x_b_c) + Mz_c_y

% Checking equlibrium for whole problem
X_eq = sum(q.*s.*cos(theta)) + Fx_c
Y_eq = sum(q.*s.*sin(theta)) + Fy_c
Z_eq = sum(stress_b.*A) - Fz_c
M_eq = -sum(q.*s.*cos(theta).*y_b_c) + sum(q.*s.*sin(theta).*x_b_c) + Mz_c

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


