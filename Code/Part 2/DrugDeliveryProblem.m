close all; clear;
%% Initialise variables, material parameters and boundary conditions
% Important x coordinates
E = 0.00166667;
D = 0.005;
B = 0.01;

method = 'Crank-Nicholson';                 % Method
dt = 0.01;                                  % Time step
T = 30;                                     % Total runtime (s)
N = T/dt;                                   % No. of timesteps
Ne = 12;                                    % Number of elements
Ngq = 3;                                    % GQ level to use
gq = CreateGQScheme(Ngq);                   % Initialise GQ scheme
msh = OneDimMeshGen(0, B, Ne);              % 1D mesh
dx = msh.nvec(2);                           % Position Step
param.material.f = 0;                       % Source term
param.material.D = [4e-6, 5e-6, 2e-6];      % Diffusion coefficient
param.material.beta = [0, 0.01, 0.01];      % Extra-vascular diffusivity
param.material.gamma = [0.02, 0.02, 0.02];  % Drug degradation rate
param.material.lambda = 0;                  % Linear reaction coefficient
param.BC.dirichlet.start = 30;              % Dirichlet BC at x=0
param.BC.dirichlet.end = 0;                 % Dirichlet BC at x=1
param.BC.neumann.start = NaN;               % Neumann BC at x=0
param.BC.neumann.end = NaN;                 % Neumann BC at x=1
type = 'Quadratic';                         % Basis function to use

% Find last node corresponding to tissue layer boundaries
param.layers.E = find(msh.nvec <= E, 1, 'last') - 1;
param.layers.D = find(msh.nvec <= D, 1, 'last') - 1;

% Check stability criterion
criterion = (param.material.D*dt)/(dx^2);
if criterion > 1/2
    fprintf('WARNING: Solution may be unstable, consider decreasing time step or increasing position step\n');
end
%% Solve problem 
C = TransientFEMSolver(method, N, dt, msh, gq, param, type);

%% Plot at specified times
% Finding column indices of time points in question
time_points_of_interest = [1:30];
column_indices = time_points_of_interest / dt;

% Determine what to plot the result against
switch type
    case 'Linear'
       x = msh.nvec;
    case 'Quadratic'
       x = linspace(0, 0.01, (2*Ne)+1);
end
       

% Plot the solution at the specified time points
figure (1);
hold on
for i = 1:length(column_indices)
%     plot_handles_t(i) = plot(msh.nvec, C(:, column_indices(i)).');
    plot_handles_t(i) = plot(x, C(:, column_indices(i)).');
end
hold off
% Add axes labels
xlabel('x');
ylabel('c(x, t)');
%xlim([0, 1]); xticks(0:0.2:1);
%ylim([0, 1]); yticks(0:0.2:1);
title(method);
grid on 

% Add legend
lgd_t = legend(plot_handles_t, cellfun(@(x) [num2str(x)], num2cell(time_points_of_interest), 'UniformOutput', false), 'Location', 'northeast');
title(lgd_t, 't');

% Adjust figure size
figure_size = [300, 300, 700, 400];  % [left, bottom, width, height]
set(gcf, 'Position', figure_size);

%% Plot at specified positions
% % Finding row indices of position points in question
% pos_points_of_interest = [0.8];
% row_indices = pos_points_of_interest / dt;
% 
% % Define time array
% t = linspace(0, 1, N);
% 
% % Plot the solution at the specified time points
% figure (2);
% hold on
% for i = 1:length(row_indices)
%     plot_handles_x(i) = plot(t, C(row_indices(i), :).');
% end
% hold off
% % Add axes labels
% xlabel('t');
% ylabel('c(x, t)');
% %xlim([0, 1]); xticks(0:0.2:1);
% %ylim([0, 1]); yticks(0:0.2:1);
% title(method);
% grid on
% 
% % Add legend
% lgd_x = legend(plot_handles_x, cellfun(@(x) [num2str(x)], num2cell(pos_points_of_interest), 'UniformOutput', false), 'Location', 'northwest');
% title(lgd_x, 'x');
% 
% % Adjust figure size
% figure_size = [300, 300, 700, 400];  % [left, bottom, width, height]
% set(gcf, 'Position', figure_size);