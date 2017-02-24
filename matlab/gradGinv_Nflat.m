function grad = gradGinv_Nflat(x, data_N)
    grad = zeros(size(x));
    grad(:,2) = x(:,2) / data_N;
end
