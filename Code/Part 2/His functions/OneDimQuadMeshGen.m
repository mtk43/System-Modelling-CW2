function [mesh] = OneDimQuadMeshGen(xmin, xmax, Ne)
    % This function generates a one-dimensional, equispaced, quadratic finite
    % element mesh with Ne number of elements between the points at x
    % position xmin and xmax.

    mesh.ne = Ne;                   % Set the number of elements
    mesh.ngn = 2 * Ne + 1;          % Set the number of global nodes for quadratic elements
    mesh.nvec = zeros(mesh.ngn, 1); % Allocate vector to store global node values
    dx = (xmax - xmin) / Ne;        % Calculate element size

    % Calculate the global node values for quadratic elements
    mesh.nvec = xmin:dx/2:xmax;

    % Loop over elements and set the element properties
    for i = 1:Ne
        % Set spatial positions of nodes
        mesh.elem(i).x(1) = xmin + (i - 1) * dx;
        mesh.elem(i).x(2) = xmin + (2 * i - 1) * dx / 2;
        mesh.elem(i).x(3) = xmin + i * dx;

        % Set global IDs of the nodes
        mesh.elem(i).n(1) = 2 * i - 1;
        mesh.elem(i).n(2) = 2 * i;
        mesh.elem(i).n(3) = 2 * i + 1;

        % Set element Jacobian based on mapping to a standard element
        mesh.elem(i).J = dx / 2; % Assuming a standard element from -1 to 1
    end
end
