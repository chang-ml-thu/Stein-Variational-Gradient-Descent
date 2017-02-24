function grad = gradGinv_warpN(x, data_N, sig0, sigx)
    coeff = data_N * (sigx^2 + data_N * sig0^2);
    grad = zeros(size(x));
    grad(:,1) = 1 ./ (((sigx/sig0)^2 + data_N + 4*data_N*(x(:,2).^2)).^2);
    grad(:,2) = grad(:,1);
    grad(:,1) = grad(:,1) .* ( -2*coeff + 8 * data_N^2 * sig0^2 * (x(:,2).^2) );
    grad(:,2) = grad(:,2) .* ( -8*coeff * x(:,2) );
end
