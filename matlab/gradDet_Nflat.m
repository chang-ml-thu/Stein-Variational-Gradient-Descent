function grad = gradDet_Nflat(x)
	grad = zeros(size(x));
    grad(:,2) = -2 ./ x(:,2);
end
