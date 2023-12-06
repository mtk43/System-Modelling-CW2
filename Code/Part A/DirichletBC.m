function [GM, GV] = DirichletBC(cond0, cond1, GM, GV)
%%Needs to be heaviliy generalised
%% Apply condition at x = 0
GM(1, :) = 0; % Set the last row to zeros
GM(1, 1) = 1; % Set the diagonal element to 1
GV(1, :) = cond0;

%% Apply condition at x = 1
GM(end, :) = 0; % Set the last row to zeros
GM(end, end) = 1; % Set the diagonal element to 1
GV(end, :) = cond1;

end