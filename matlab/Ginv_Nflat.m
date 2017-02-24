function G = Ginv_Nflat(x, data_N)
    M = size(x,1);
    G = repmat([1 0; 0 0.5]/data_N, [1 1 M]);
    G = bsxfun(@times, G, reshape(x(:,2).^2, [1 1 M]));
end
