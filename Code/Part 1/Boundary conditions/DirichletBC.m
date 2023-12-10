function [GM, GV] = DirichletBC(cond0, cond1, GM, GV)
%%Needs to be heaviliy generalised
%% Apply condition at x = 0
if ~isnan(cond0)
    GM(1, :) = 0; % Set the last row to zeros
    GM(1, 1) = 1; % Set the diagonal element to 1
    GV(1, :) = cond0;
end

%% Apply condition at x = 1
if ~isnan(cond1)
    GM(end, :) = 0; % Set the last row to zeros
    GM(end, end) = 1; % Set the diagonal element to 1
    GV(end, :) = cond1;
end

end