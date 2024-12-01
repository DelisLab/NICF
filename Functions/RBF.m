function K = RBF(X, Y, gamma)
    % X and Y are input matrices where each row is a data point.
    % gamma is the RBF kernel parameter.

    % Compute the squared Euclidean distance between each pair of points
    % using the formula ||x - y||^2 = ||x||^2 + ||y||^2 - 2*x*y'
    
    % X and Y are (n x d) and (m x d) matrices, respectively
    X_norm = sum(X.^2, 2); % Column vector of squared norms of X
    Y_norm = sum(Y.^2, 2); % Column vector of squared norms of Y
    
    % Compute the squared Euclidean distance matrix
    dist = bsxfun(@plus, X_norm, Y_norm') - 2 * (X * Y');
    
    % Compute the RBF kernel matrix
    K = exp(-gamma * dist);
end