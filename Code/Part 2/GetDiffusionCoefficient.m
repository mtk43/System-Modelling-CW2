function D = GetDiffusionCoefficient(D_coords, x)
    % GetDiffusionCoefficient returns the diffusion coefficient for the given coordinates
    % based on the specified values and coordinates.
    %
    % Inputs:
    %   D_coords - matrix of diffusion coefficient values and corresponding coordinates
    %   x - x-coordinate of the current element
    %
    % Output:
    %   D - diffusion coefficient for the current element

    % Assuming D_coords is sorted based on coordinates
    [~, index] = max(x >= D_coords(:, 2));
    D = D_coords(index, 1);
end
