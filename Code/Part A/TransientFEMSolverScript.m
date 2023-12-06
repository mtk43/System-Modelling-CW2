% Step 1: Initialise mesh
% Assuming you have a function OneDimLinearMeshGen that generates the mesh
N = 10; % Number of elements
msh = OneDimLinearMeshGen(0, 1, N);

% Step 2: Initialise time integration scheme
theta = 0.5; % Theta |value for time integration (0.5 is for Crank-Nicolson)
dt = 0.1;   % Time step
N_steps = 50; % Number of time steps
T = N_steps * dt; % Total time

% Step 3: Define material coefficients
D = 1.0;     % Diffusion coefficient
lambda = 0.1;% Some parameter
f = @(x) sin(pi * x); % Source term function

% Step 4: Initialise Global Matrices and Vectors
GM = zeros(msh.nNodes, msh.nNodes);
M = zeros(msh.nNodes, msh.nNodes);
K = zeros(msh.nNodes, msh.nNodes);
GV = zeros(msh.nNodes, 1);

% Step 5: Define solution variable vectors
Ccurr = zeros(msh.nNodes, 1);
Cnext = zeros(msh.nNodes, 1);

% Step 6: Set initial conditions on Ccurrent
Ccurr = sin(pi * msh.nvec);

% Step 7: Time loop
for tstep = 1:N_steps
    % Step 7.1: Loop over elements
    for e = 1:msh.nElements
        % Calculate local element mass matrix and stiffness matrix
        Me = LocalMassMatrix(msh, e);
        Ke = LocalStiffnessMatrix(msh, e, D);
        
        % Add to global mass and stiffness matrices
        M(msh.index(e), msh.index(e)) = M(msh.index(e), msh.index(e)) + Me;
        K(msh.index(e), msh.index(e)) = K(msh.index(e), msh.index(e)) + Ke;
    end
    
    % Step 7.2: Calculate global matrix (GM)
    GM = M + theta * dt * K;
    
    % Step 7.3: Calculate global vector (GV)
    GV = M * Ccurr + dt * theta * K * Ccurr;
    
    % Step 7.4: Loop over elements for additional contributions
    for e = 1:msh.nElements
        % Calculate local element source vector
        Fe = LocalSourceVector(msh, e, f);
        
        % Add to global vector
        GV(msh.index(e)) = GV(msh.index(e)) + dt * Fe;
        
        % Add Neumann boundary conditions if specified
        % ...
    end
    
    % Step 7.5: Apply Dirichlet Boundary Conditions
    % ...
    
    % Step 7.6: Solve the matrix system to obtain Cnext
    Cnext = GM \ GV;
    
    % Step 8: Set Ccurrent equal to Cnext
    Ccurr = Cnext;
    
    % Step 9: Re-initialise global matrices and vectors to zero
    GM = zeros(msh.nNodes, msh.nNodes);
    M = zeros(msh.nNodes, msh.nNodes);
    K = zeros(msh.nNodes, msh.nNodes);
    GV = zeros(msh.nNodes, 1);
    
    % Step 10: Plot/write to file the solution Ccurrent
    plot(msh.nvec, Ccurr);
    xlabel('X-axis');
    ylabel('Solution');
    title(['Solution at t = ', num2str(tstep * dt)]);
    pause(0.1); % Pause to visualize the plot
end
