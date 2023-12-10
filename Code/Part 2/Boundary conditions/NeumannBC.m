function [neumann_vector] = NeumannBC(cond0, cond1, D, GV_curr)
%%Function to apply Neumann BC (NBC) to an FEM problem
% Takes in the current global vector (GV_curr), the NBC at the position
% boundaries of the problem, (cond0 and cond1) and the diffusion
% coefficient (D) and outputs the new global vector (GV_new)
% values) and

%% Initialise Neumann vector
neumann_vector = zeros(length(GV_curr), 1);

%% Apply condition at x = 0
% Determine whether NBC has been set
if ~isnan(cond0)
    neumann_vector(1) = -D*cond0;
end

%% Apply condition at x = 1
% Determine whether NBC has been set
if ~isnan(cond1)
    neumann_vector(end) = D*cond1;
end

% %% Calculate new global vector accounting for NBC
% GV_new = GV_curr + neumann_vector;

end