close all; clear; clc;
% This Scrip Finds the x and y position of the booms from the leading edge.

% Axis follow Right hand rule
% Positive x is from wing leading edge to trailing edge
% Positive y is up
% Positive z is along left wing

% Code finds position for booms on the top surface of the aerofoil and the
% borrom surface of the aerofoil.

%
% Initial Inputs
%

t_c = 0.15; % Thickness to Chord Ratio (-)
p = 0.2; % Position of Max Camber (-)
m = 0.02; % Camber (-)

c = 2.2; % Chord Length (m)
b_span = 16; % Wing Span (m)
n = 500000; % Number of Points (-)

% x-coordinate of booms. Top row is for top surface and Bottom row is for
% bottom surface.
x_booms_top = [0.35, 0.75, 0.8, 1.1];
x_booms_bottom = [0.35, 0.75];


%
% Creating Symmetric Aerofoil
%

% x-coorinates of symmetric aerofoil
x_o = linspace(0, c, n);

% Thickness of aerofoil with respect to x
y_t = 5*t_c*c*(0.2969*(x_o/c).^0.5 - 0.126*(x_o/c) - 0.3516*(x_o/c).^2 + 0.2843*(x_o/c).^3 - 0.1036*(x_o/c).^4);

% This plots the camber line of the aerofoil, y_c, along with the gradient
% of the camber line of the aerofoil, dy_c. Using the gradient of the
% camber line you can find theta which is used to find the coordinates of
% the NACA aerofoil.
for i=1:n
    if 0 <= x_o(i) & x_o(i) < p*c
        y_c(i) = m*(2*p-x_o(i)/c)*x_o(i)/p^2;
        dy_c(i) = 2*m*(p - x_o(i)/c)/p^2;
    elseif p*c <= x_o(i) & x_o(i) <= c
        y_c(i) = m*(c-x_o(i))*(1 + x_o(i)/c - 2*p)/(1-p)^2;
        dy_c(i) = 2*m*(p-x_o(i)/c)/(1-p)^2;
    end

    theta(i) = atan(dy_c(i));

end

%
% Calcaulting the coordinates of the aerofoil
%

% Upper Surface
X_U = x_o - y_t.*sin(theta);
Y_U = y_c + y_t.*cos(theta);

% Lower Surface
X_L = x_o + y_t.*sin(theta);
Y_L = y_c - y_t.*cos(theta);

% Putting the coordinates together
X = [X_U, flip(X_L)];
Y = [Y_U, flip(Y_L)];


%
% Finding the y position of the booms
%

% Gets the size of the boom array so for loop will go through all elements
[h_t, w_t] = size(x_booms_top);
[h_b, w_b] = size(x_booms_bottom);

% Boom number array
n = linspace(1, w_t+w_b, w_t+w_b);
count = 0;
for j = 1:2
    for i = 1:w_t      
        if j == 1
            count = count + 1;
            % Finds the x of the boom that is closes to the aerofoil
            [val, ind] = min(abs(X_U - x_booms_top(i)));
            % Uses the index, ind, to find y value and stores in array
            x_b_top(count) = X_U(ind);
            y_b_top(count) = Y_U(ind);
        end
    end
    for i = 1:w_b
        if j == 2
            count = count + 1;
            % Finds the x of the boom that is closes to the aerofoil
            [val, ind] = min(abs(X_L - x_booms_bottom(i)));
            % Uses the index, ind, to find y value and stores in array
            x_b_bottom(count - w_t) = X_L(ind);
            y_b_bottom(count - w_t) = Y_L(ind);
        end
    end
end

% Combines respective positions of booms into one array
x_b = [x_b_top, flip(x_b_bottom)];
y_b = [y_b_top, flip(y_b_bottom)];

% Makes one big array with all the data
boom = [n', x_b', y_b'];

disp(boom)


f = figure(1);
f.Position = [50, 50, 1000, 300];
axis([0, 2.2, -0.15, 0.25]);
figure(f)
hold on
plot(X, Y, 'g')
plot([x_b, x_b(1)], [y_b, y_b(1)], 'c', 'MarkerFaceColor', 'b')
plot([x_b, x_b(1)], [y_b, y_b(1)], 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold off
grid on
set(gca, 'XTick', 0:0.1:2.2)
set(gca, 'YTick', -0.15:0.1:0.25)
daspect([1, 1, 1]);









