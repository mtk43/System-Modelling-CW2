theta = 1/2;                            % Crank-Nicholson
dt = 0.001;                             % Time step
N = 1/dt;                               % No. of timesteps
f = 0;                                  % Source term
D = 1;                                  % Diffusion term
lambda = 0;                             % Linear reaction term
Ne = 10;                                % Number of elements
gq = CreateGQScheme(2);                 % Initialise GQ scheme
msh = OneDimLinearMeshGen(0, 1, Ne);    % 1D mesh

Ccurr = TransientFEMSolver(theta, N, dt, f, D, lambda, msh, gq);

% Plot the solution Ccurrent
plot(msh.nvec, Ccurr);
xlabel('X-axis');
ylabel('Solution');
title(['Solution at t = ', num2str(N * dt)]);
pause(0.1); % Pause to visualize the plot



% Plot analytical solution
hold on
x = 0.5;
c = TransientAnalyticSoln(x, 0.1);
scatter(x, c);
hold off