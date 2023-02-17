clear;
clc;
close all;
format long;

%% Initializing Constants
Cm = 0.01;
Cp = 2500;
delta = 2;
k = 0.65;
ro = 370;
Dm = 2.2e-8;
Tb = 10;
Ti = 60;
T0 = 110;
l = 0.024;
hl = 2.5e-9;
mb = 86;
mi = 45;
m0 = 4;
gamma = 0;
epsilon = 0.3;

%% User Defined Variables
totalT = 6000;

xStep = 0.001;
tStep = 120;

xNum = (l / xStep);
tNum = totalT / tStep;

%% Initializing Result Matrixes
% We use these two matrixes to hold our results,
% Each column has 24 values for different X values in the same time
% Each row has 50 values for the same X with different times
% Boundary values is stored in the first column of each matrix
temp = zeros(xNum + 1, tNum);
mois = zeros(xNum + 1, tNum);

for i = 1 : tNum
    temp(1, i) = Ti;
    mois(1, i) = mi;
    temp(xNum + 1, i) = T0;
    mois(xNum + 1, i) = m0;
end

for i = 1 : xNum + 1
    temp(i, 1) = Tb;
    mois(i, 1) = mb;
end


%% Initializing Equations
% L, D, v, lambda values is calculated using the equations in the article
L = k / ro * Cp;
D = (k * Dm) / (ro * Cm * (k + Dm * delta * (epsilon * hl + gamma)));
v = Cm * (epsilon * hl + gamma) / Cp;
lambda = (Cp * Dm * delta) / (Cm * (k + Dm * delta * (epsilon * hl + gamma)));

%% Alpha and Beta are two Constants
alpha = D * tStep / xStep * xStep;
alpha = alpha * 1e6;

beta = L * tStep / xStep * xStep;

%% Initializing A and B Matrixes in AX = B Equation
%createMatrixA;
A = zeros((xNum - 1) * 2, (xNum - 1) * 2);

% Using the First Equation 
for i = 2 : xNum
	row = i - 1;
	if i == 2                   % Lower Boundary
		% Temp Variables
		A(row, i - 1) = -1 * ((2 * beta) + 1);
		A(row, i) = beta;
		
		% Mois Variables
		A(row, i + xNum - 2) = v;     
	elseif i == xNum            % Upper Boundary
		% Temp Variables
		A(row, i - 1) = -1 * ((2 * beta) + 1);
		A(row, i - 2) = beta;
		
		% Mois Variables
		A(row, i + xNum - 2) = v;
    else                        % Otherwise
		% Temp Variables
		A(row, i - 1) = -1 * ((2 * beta) + 1);
		A(row, i) = beta;
		A(row, i - 2) = beta;
		
		% Mois Variables
		A(row, i + xNum - 2) = v;
	end
end

% Using the Second Equation
for i = 2 : xNum
	row = i + xNum - 2;
	if i == 2                   % Lower Boundary            
		% Mois Variables
		A(row, i - 1 + xNum - 1) = -1 * ((2 * alpha) + 1);
		A(row, i + xNum - 1) = alpha;
		
		% Temp Variables
		A(row, i - 1) = lambda;                        
	elseif i == xNum            % Upper Boundary            
		% Mois Variables
		A(row, i - 1 + xNum - 1) = -1 * ((2 * alpha) + 1);
		A(row, i - 2 + xNum - 1) = alpha;
		
		% Temp Variables
		A(row, i - 1) = lambda;
    else                        % Otherwise          
		% Mois Variables
		A(row, i - 1 + xNum - 1) = -1 * ((2 * alpha) + 1);
		A(row, i + xNum - 1) = alpha;
		A(row, i - 2 + xNum - 1) = alpha;
		
		% Temp Variables
		A(row, i - 1) = lambda;
	end
end

B = zeros((xNum - 1) * 2, 1);

%% Calculation
% Here we use these two FORs to irativley calculate different matrixes
for j = 1 : tNum - 1
    %updateMatrixB;
    % Using the First Equation 
for i = 2 : xNum
	row = i - 1;
	if i == 2
		B(row) = (-beta) * temp(i - 1, j + 1) - temp(i, j) + v * mois(i, j);
    elseif i == xNum
		B(row) = v * mois(i, j) - beta * temp(i + 1, j + 1) - temp(i, j);
	else
		B(row) = v * mois(i, j) - temp(i, j);
	end
end

% Using the Second Equation
for i = 2 : xNum
	row = i + xNum - 2;
	if i == 2
		B(row) = (-alpha) * mois(i - 1, j + 1) - mois(i, j) + lambda * temp(i, j);
	elseif i == xNum
		B(row) = lambda * temp(i, j) - alpha * mois(i + 1, j + 1) - mois(i, j);
	else
		B(row) = lambda * temp(i, j) - mois(i, j);
    end
end
    x = A \ B;
   for i = 1 : xNum - 1
   temp(i + 1, j + 1) = x(i, 1);
   mois(i + 1, j + 1) = x(i + xNum - 1, 1);
end
end
%%%%%%%%%%%%%%%%%%%
figure;
hold on;
x = 1 : tNum;
x = x * tStep;
y = temp(12, :);
set(gca, 'xlim', [0 30*tStep]);
set(gca, 'ylim', [0 100]);
plot(x, y, 'LineWidth', 2)
title('Temperature - Time (X = 0.012)');
xlabel('Time');
ylabel('Temp.');
%%%%%%%%%%%%%%%%
figure;
hold on;
x = 1 : tNum;
x = x * tStep;
y = mois(12, :);
set(gca, 'xlim', [0 30*tStep]);
set(gca, 'ylim', [80 90]);
plot(x, y, 'LineWidth', 2)
title('Moisture - Time (X = 0.012)');
xlabel('Time');
ylabel('Mois.');
%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
cmap = hsv(6);
for i = 2 : 6
    x = 0 : xNum;
    x = x / 1000;
    y = temp(:, i);
    set(gca, 'xlim', [0 0.025]);
    set(gca, 'ylim', [0 110]);
    plot(x, y, 'Color',cmap(i, :), 'MarkerSize', 10, 'LineWidth', 2);
end
title('Temperature - X (Different Times)');
xlabel('X');
ylabel('Temp.');
%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
cmap = hsv(tNum);
for i = 2 : tNum
    x = 0 : xNum;
    x = x / 1000;
    y = mois(:, i);
    set(gca, 'xlim', [0 0.03]);
    set(gca, 'ylim', [0 120]);
    plot(x, y, 'Color',cmap(i, :), 'MarkerSize', 10, 'LineWidth', 0.5);
end
title('Moisture - X (Different Times)');
xlabel('X');
ylabel('Mois.');
