function [C] = TransientFEMSolver(method, N, dt, msh, gq, param, type)

%% Initialise all the necessary variables
% Define material parameter values
D = param.material.D;            % Diffusion coefficient;
f = param.material.f;            % Source term
lambda = param.material.lambda;  % Linear reaciton term

% Initialise the necessary matrices
switch type
    case 'Linear'
        GM = zeros(msh.ngn);     % Global Matrix
        M = zeros(msh.ngn);      % Global Mass Matrix
        K = zeros(msh.ngn);      % Global Stiffness Matrix
        GV = zeros(msh.ngn, 1);  % Global vector
        Ccurr = zeros(msh.ngn, 1);   % Current solution
        Cnext = zeros(msh.ngn, 1);   % Next solution
        C = zeros(msh.ngn, N);       % Solution matrix to store all Ccurr values
    
    case 'Quadratic'
        GM = zeros((2*msh.ne)+1);     % Global Matrix
        M = zeros((2*msh.ne)+1);      % Global Mass Matrix
        K = zeros((2*msh.ne)+1);      % Global Stiffness Matrix
        GV = zeros((2*msh.ne)+1, 1);  % Global vector
        Ccurr = zeros((2*msh.ne)+1, 1);   % Current solution
        Cnext = zeros((2*msh.ne)+1, 1);   % Next solution
        C = zeros((2*msh.ne)+1, N);       % Solution matrix to store all Ccurr values
end

% Initialise solution vectors
% Ccurr = zeros(msh.ngn, 1);   % Current solution
% Cnext = zeros(msh.ngn, 1);   % Next solution
% C = zeros(msh.ngn, N);       % Solution matrix to store all Ccurr values

%% Determine the method to use
switch method
    case 'Crank-Nicholson'
        theta = 1/2;
    case 'Forward Euler'
        theta = 0;
    case 'Backward Euler'
        theta = 1;
    otherwise
        fprintf('Invalid method specified\n');
end

%% Determine global local mass & stiffness matrices
% Step Loop over elements
start_index = 1;
for elem = 1:msh.ne
    % Calculate local element mass matrix and stiffness matrix
    Melem = LocalMassMatrix(elem, msh, gq, type);
    Kelem = LocalStiffnessMatrix(elem, msh, D, lambda, gq, type);
     
    % Update the global mass matrix
    M(start_index:(start_index + length(Melem))-1, start_index:(start_index + length(Melem))-1) = ...
        M(start_index:(start_index + length(Melem))-1, start_index:(start_index + length(Melem))-1) + Melem;
    
    % Update the global stiffness matrix
    K(start_index:(start_index + length(Kelem))-1, start_index:(start_index + length(Kelem))-1) = ...
        K(start_index:(start_index + length(Kelem))-1, start_index:(start_index + length(Kelem))-1) + Kelem;
    
    start_index = start_index + (length(Melem)-1);
end

% Calculate global matrix (GM)
GM = M + (theta * dt * K);

%% Loop through elements at each time interval
for tstep = 1:N 
    
    % Calculate global vector (GV)
    GV = (M - (1 - theta) * dt * K) * Ccurr;
    
    % Apply Neumann Boundary Conditions - assumes they are constant in time
    neumann_vector = NeumannBC(param.BC.neumann.start, param.BC.neumann.end, D, GV);
    GV = GV + dt*(neumann_vector);
    
    % Apply Dirichlet Boundary Conditions
    [GM, GV] = DirichletBC(0, 1, GM, GV);
    
    % Solve the matrix system to obtain Cnext
    Cnext = GM \ GV;
    
    % Set Ccurrent equal to Cnext
    Ccurr = Cnext;
    
    % Save all the Ccurr value to a matrix, C
    C(:, tstep) = Ccurr;
    
    % Re-initialise global vector to zero
    GV = zeros(msh.ne, 1);   
    Cnext = zeros(msh.ngn, 1);
    
    % Plot the solution every n time steps 
%     if(mod(tstep,10)==0) 
%         x = msh.nvec; 
%         hold on
%         figure (1); 
%         plot(x, Ccurr, 'DisplayName', ['tstep = ' num2str(tstep)]);
%     end
end
% legend('show');
% hold off

% C(end, 1) = 1;
% for tstep = 1:N
% 
%     % Calculate global vector (GV)
%     GV = (M - (1 - theta) * dt * K) * C(:, tstep);
% 
%     % Loop over elements for additional contributions - for now there
%     % is no source term so not including it
% %     for elem = 1:Ne
% %         % Calculate local element source vector
% %         Fe = LocalSourceVector(msh, elem, f);
% % 
% %         % Add to global vector
% %         GV(msh.index(elem)) = GV(msh.index(elem)) + dt * Fe;
% % 
% %         % Add Neumann boundary conditions if specified
% %         % ...
% %     end
%     % Apply Neumann Boundary Conditions - assumes they are constant in time
%     neumann_vector = NeumannBC(param.BC.neumann.start, param.BC.neumann.end, D, GV);
%     GV = GV + dt*(neumann_vector);
% 
%     % Apply Dirichlet Boundary Conditions
%     [GM, GV] = DirichletBC(param.BC.dirichlet.start, param.BC.dirichlet.end, GM, GV);
% 
%     % Solve the matrix system to obtain C at the next time step 
%     C(:, tstep+1) = GM \ GV;
% 
% end
end