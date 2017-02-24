function grad = dlog_p_warpN(x, m1, data_N, sig0, sigx)
	grad = zeros(size(x));
    grad(:,1) = (data_N * (x(:,2).^2 + x(:,1)) - m1)  / sigx^2;
    grad(:,2) = - x(:,2)/sig0^2 - 2*grad(:,1).*x(:,2);
    grad(:,1) = - x(:,1)/sig0^2 - grad(:,1);
end
