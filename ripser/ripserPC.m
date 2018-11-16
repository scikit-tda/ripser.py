function [Is, Cocycles] = ripserPC( X, varargin )
    % Run ripser on a Euclidean point cloud
    %:param X: N X d point cloud with N points in d dimensions
    %:param maxdim: Maximum dimension of homology
    %:param thresh: Threshold up to which to add edges (default max dist*2)
    %:param coeff: Field coefficient with which to run TDA (default 2)
    p = inputParser;
    addParameter(p, 'maxdim', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'thresh', inf, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'coeff', 2, @(x) isprime(x));
    parse(p, varargin{:});
    XSqr = sum(X.^2, 2);
    D = bsxfun(@plus, XSqr(:), XSqr(:)') - 2*(X*X');
    D = 0.5*(D + D');
    D(D < 0) = 0;
    D = sqrt(D);
    [Is, Cocycles] = ripserDM(D, varargin{:});
end

