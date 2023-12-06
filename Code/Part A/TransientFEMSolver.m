function [Ccurr] = TransientFEMSolver(theta, N, dt, f, D, lambda, msh, gq)

%% Initialise all the necessary variables

% Initialise the necessary matrices
GM = zeros(msh.ngn);     % Global Matrix
M = zeros(msh.ngn);             % Global Mass Matrix
K = zeros(msh.ngn);              % Global Stiffness Matrix
GV = zeros(msh.ne, 1);     % Global vector

% Initialise solution vectors
Ccurr = zeros(msh.ngn, 1);       % Current solution
Cnext = zeros(msh.ngn, 1);        % Next solution

%% Determine global local mass & stiffness matrices
% Step Loop over elements
for elem = 1:msh.ne
    % Calculate local element mass matrix and stiffness matrix
    Melem = LocalMassMatrix(elem, msh, gq);
    Kelem = LocalStiffnessMatrix(elem, msh, D, lambda, gq);

    % Create an offset based on the element number
    offset = (elem - 1) * (length(Melem) - 1);

    % Update the global mass matrix
    M((elem:length(Melem) + offset), (elem:length(Melem) + offset)) = ...
        M((elem:length(Melem) + offset), (elem:length(Melem) + offset)) + Melem;

    % Update the global stiffness matrix
    K((elem:length(Kelem) + offset), (elem:length(Kelem) + offset)) = ...
        K((elem:length(Kelem) + offset), (elem:length(Kelem) + offset)) + Kelem;
end

% Calculate global matrix (GM)
GM = M + (theta * dt * K);

%% Loop through elements at each time interval
for tstep = 1:N 
       
    % Calculate global vector (GV)
    GV = (M - (1 - theta) * dt * K) * Ccurr;
    
    % Loop over elements for additional contributions - for now there
    % is no source term so not including it
    %         for elem = 1:Ne
    %             % Calculate local element source vector
    %             Fe = LocalSourceVector(msh, elem, f);
    %
    %             % Add to global vector
    %             GV(msh.index(elem)) = GV(msh.index(elem)) + dt * Fe;
    %
    %             % Add Neumann boundary conditions if specified
    %             % ...
    %         end
    
    % Apply Dirichlet Boundary Conditions
    [GM, GV] = DirichletBC(0, 1, GM, GV);
    
    % Solve the matrix system to obtain Cnext
    Cnext = GM \ GV;
    
    % Set Ccurrent equal to Cnext
    Ccurr = Cnext;
    
    % Re-initialise global vector to zero
    GV = zeros(msh.ne, 1);   
    Cnext = zeros(msh.ngn, 1);
    
    % Plot the solution every 5 time steps 
% if(mod(tstep,5)==0) 
%     x = msh.nvec; 
%     hold on
%     figure (1); 
%     plot(x,Cnext); 
% end
end
hold off

end