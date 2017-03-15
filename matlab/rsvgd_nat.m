function  theta = rsvgd_nat(theta0, dlog_p, gradDet, Ginv, max_iter, master_stepsize, h, auto_corr, method)

%%%%%%%%
% Bayesian Inference via Stein Variational Gradient Descent

% input:
%   -- theta0: initialization of particles, m * d matrix (m is the number of particles, d is the dimension)
%   -- dlog_p: function handle of first order derivative of log p(x)
%   -- max_iter: maximum iterations
%   -- master_stepsize: the general learning rate for adagrad
%   -- h/bandwidth: bandwidth for rbf kernel. Using median trick as default
%   -- auto_corr: momentum term
%   -- method: use adagrad to select the best \epsilon

% output:
%   -- theta: a set of particles that approximates p(x)
%%%%%%%%

if ~exist('master_stepsize','var') || isempty(master_stepsize); master_stepsize = 0.1; end;
if ~exist('h','var') || isempty(h); h = -1; end;
if ~exist('auto_corr','var') || isempty(auto_corr); auto_corr = 0.9; end;
if ~exist('method','var') || isempty(method); method = 'adagrad'; end;

switch lower(method)
    
    case 'adagrad'
        %% AdaGrad with momentum
        theta = theta0;
        
        fudge_factor = 1e-6;
        historial_grad = 0;
        
        for iter = 1:max_iter
            grad = OptVecField(theta, dlog_p, gradDet, Ginv, h);   %\Phi(theta)
            if historial_grad == 0
                historial_grad = historial_grad + grad.^2;
            else
                historial_grad = auto_corr * historial_grad + (1 - auto_corr) * grad.^2;
            end
            adj_grad = grad ./ (fudge_factor + sqrt(historial_grad));
            theta = theta + master_stepsize * adj_grad; % update
        end
        
	case 'fixed'
		theta = theta0;
		for iter = 1:max_iter
            grad = OptVecField(theta, dlog_p, gradDet, Ginv, h);   %\Phi(theta)
            theta = theta + master_stepsize * grad; % update
        end
        
    otherwise
        error('wrong method');
end
end


function [Akxy, info] = OptVecField(x, dlog_p, gradDet, Ginv, h)
%%%%%%%%%%%%%%%%%%%%%%
% Input:
%    -- x: particles, n*d matrix, where n is the number of particles and d is the dimension of x 
%    -- dlog_p: a function handle, which returns the first order derivative of log p(x), n*d matrix
%    -- h: bandwidth. If h == -1, h is selected by the median trick

% Output:
%    --Akxy: n*d matrix, \Phi(x) is our algorithm, which is a smooth
%    function that characterizes the perturbation direction
%    --info: kernel bandwidth
%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5; h = -1; end % median trick as default

[n, d] = size(x);
% Using rbf kernel as default
XY = x*x';
x2= sum(x.^2, 2);
X2e = repmat(x2, 1, n);
H = (X2e + X2e' - 2*XY); % calculate pairwise distance

% median trick for bandwidth
if h == -1
    h = sqrt(0.5*median(H(:)) / log(n+1));   %rbf_dot has factor two in kernel
end

Kxy = exp(-H/(2*h^2));   % calculate rbf kernel
ginv = Ginv(x);
trans_score = dlog_p(x) + gradDet(x);
GinvDxKxy = zeros(n,d);
for j = 1:n
	trans_score(j,:) = trans_score(j,:) * ginv(:,:,j);
	dxKxy = bsxfun(@times, Kxy(:,j)/(h^2), bsxfun(@minus, x(j,:), x));
	for i = 1:n
		GinvDxKxy(j,:) = GinvDxKxy(j,:) + dxKxy(i,:) * ginv(:,:,i);
	end
end

Akxy = (Kxy * trans_score + GinvDxKxy)/n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info.bandwidth = h;
end
