function grad = dlog_p_Nflat(x, m1, m2, data_N)
	grad = zeros(size(x));
    grad(:,1) = (m1 - data_N*x(:,1)) ./ (x(:,2).^2);
    grad(:,2) = (m2 - 2*m1*x(:,1) + data_N*(x(:,1).^2)) ./ (x(:,2).^3) - data_N ./ x(:,2);
end
