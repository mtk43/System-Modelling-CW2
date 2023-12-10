close all; clear;
%% Initialise variables, material parameters and boundary conditions
method = 'Crank-Nicholson';             % Method
dt = 0.001;                              % Time step
N = 1/dt;                               % No. of timesteps
Ne = 10;                                % Number of elements
gq = CreateGQScheme(2);                 % Initialise GQ scheme
msh = OneDimLinearMeshGen(0, 1, Ne);    % 1D mesh
dx = msh.nvec(2);                       % Position Step
param.material.f = 0;                   % Source term
param.material.D = 1;                   % Diffusion coefficient
param.material.lambda = 0;              % Linear reaction coefficient
param.BC.dirichlet.start = NaN;           % Dirichlet BC at x=0
param.BC.dirichlet.end = 1;             % Dirichlet BC at x=1
param.BC.neumann.start = 1;           % Neumann BC at x=0
param.BC.neumann.end = NaN;             % Neumann BC at x=1


% Check stability criterion
criterion = (param.material.D*dt)/(dx^2);
if criterion > 1/2
    fprintf('WARNING: Solution may be unstable, consider decreasing time step or increasing position step\n');
end
%% Solve problem 
C = TransientFEMSolver(method, N, dt, msh, gq, param);

%% Plot at specified times
% Finding column indices of time points in question
time_points_of_interest = [0.1];
column_indices = time_points_of_interest / dt;

% Plot the solution at the specified time points
figure (1);
hold on
for i = 1:length(column_indices)
    plot_handles_t(i) = plot(msh.nvec, C(:, column_indices(i)).');
end
hold off
% Add axes labels
xlabel('x');
ylabel('c(x, t)');
xlim([0, 1]); xticks(0:0.2:1);
ylim([0, 1]); yticks(0:0.2:1);
title(method);
grid on 

% Add legend
lgd_t = legend(plot_handles_t, cellfun(@(x) [num2str(x)], num2cell(time_points_of_interest), 'UniformOutput', false), 'Location', 'northwest');
title(lgd_t, 't');

% Adjust figure size
figure_size = [300, 300, 700, 400];  % [left, bottom, width, height]
set(gcf, 'Position', figure_size);

%% Plot at specified positions
% Finding row indices of position points in question
pos_points_of_interest = [0.8];
row_indices = pos_points_of_interest / dx;

% Define time array
t = linspace(0, 1, N);

% Plot the solution at the specified time points
figure (2);
hold on
for i = 1:length(row_indices)
    plot_handles_x(i) = plot(t, C(row_indices(i), :));
end
hold off
% Add axes labels
xlabel('t');
ylabel('c(x, t)');
xlim([0, 1]); xticks(0:0.2:1);
ylim([0, 1]); yticks(0:0.2:1);
title(method);
grid on

% Add legend
lgd_x = legend(plot_handles_x, cellfun(@(x) [num2str(x)], num2cell(pos_points_of_interest), 'UniformOutput', false), 'Location', 'northwest');
title(lgd_x, 'x');

% Adjust figure size
figure_size = [300, 300, 700, 400];  % [left, bottom, width, height]
set(gcf, 'Position', figure_size);

%% Plot Analytical solution
t = 0.1;
x = linspace(0, 1, 100);

for i = 1:length(x)
   c_analytical(i) =  TransientAnalyticSoln(x(i), t);
end

figure (1)
hold on
plot(x, c_analytical);
hold off