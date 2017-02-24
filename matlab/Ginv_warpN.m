function G = Ginv_warpN(x, data_N, sigx, sig0)
	M = size(x,1);
    G = zeros(2,2,M);
	w2sq4 = 4*(x(:,2).^2);
    G(1,1,:) = w2sq4 * data_N * sig0^2 + sigx^2;
    G(1,2,:) = -2 * data_N * (sig0^2) * x(:,2);
    G(2,1,:) = G(1,2,:);
    G(2,2,:) = sig0^2 * data_N + sigx^2;
	G = bsxfun(@rdivide, G, reshape( (sigx/sig0)^2 + data_N * (w2sq4 + 1), [1, 1, M] ));
end
