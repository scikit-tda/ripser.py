function [Is, Cocycles] = ripserDM( D, varargin )
    % Run ripser using a distance matrix
    %:param D: N X N distance matrix (dense or sparse)
    %:param maxdim: Maximum dimension of homology
    %:param thresh: Threshold up to which to add edges (default max dist*2)
    %:param coeff: Field coefficient with which to run TDA (default 2)
    p = inputParser;
    addParameter(p, 'maxdim', 1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'thresh', inf, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'coeff', 2, @(x) isprime(x));
    parse(p, varargin{:});
    if issparse(D)
        [I, J, S] = find(D);
        I = int32(I)-1; %Matlab is 1-indexed
        J = int32(J)-1; %Matlab is 1-indexed
        S = single(S);
        N = size(D);
        [Is, Cocycles] = ripser(p.Results.coeff, p.Results.maxdim, p.Results.thresh, I, J, S, N);
    else
        %pdist surrogate
        N = size(D, 1);
        [I, J] = meshgrid(1:N, 1:N);
        d = D(I < J);
        d = single(d(:));
        [Is, Cocycles] = ripser(p.Results.coeff, p.Results.maxdim, p.Results.thresh, d);
    end
end

