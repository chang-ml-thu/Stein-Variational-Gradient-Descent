function grad = gradDet_warpN(x, data_N, sig0, sigx)
	grad = zeros(size(x));
    grad(:,2) = 4*x(:,2) ./ ( sigx^2/sig0^2/data_N + 1 + 4*(x(:,2).^2) );
end
