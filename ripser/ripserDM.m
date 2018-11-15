function [Is, Cocycles] = ripserDM( D, coeff, maxdim, thresh )
    %:param D: (N X N distance matrix
    %:param coeff: Field coefficient with which to run TDA
    %:param maxdim: Maximum dimension of homology
    %:param thresh: Threshold up to which to add edges
    if nargin < 4
        thresh = max(D(:))*2;
    end

    %pdist surrogate
    N = size(D, 1);
    [I, J] = meshgrid(1:N, 1:N);
    d = D(I < J);
    d = single(d(:));
    [Is, Cocycles] = ripser(d, coeff, maxdim, thresh);
    
    %Delete entries with numerically insignificant persistence
%     for ii = 1:length(Is)
%         I = Is{ii};
%         if numel(I) == 0
%             continue;
%         end
%         P = I(:, 2) - I(:, 1);
%         idx = P > eps;
%         Is{ii} = I(idx, :);
%         Cocycles{ii} = Cocycles{ii}(idx);
%     end
end

