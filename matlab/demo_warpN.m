clear all

% warpN
fid = fopen('data_warpN_50_0_1.txt', 'r');
if fid == -1
    disp('Cannot open file!'); return;
end
data_N = fscanf(fid, '%d', 1);
data_mu = fscanf(fid, '%f', 1);
data_sig = fscanf(fid, '%f', 1);
data_x = fscanf(fid, '%f', [data_N 1]);
fclose(fid);
m1 = sum(data_x);
m2 = data_x' * data_x;

sig0 = 2;
sigx = data_sig;
M = 10; eps = 0.001; % h = 0.1;
xbnd = [-5,3]; ybnd = [-3, 3];
truew = [0.1, 0.1];
xy_init = bsxfun(@plus, randn(M,2)*1e-1, [-2,0]);
% xy_init = randn(M,2) * sig0; % from the prior
num_iter = 10; num_cont = 10; ssize = 150;

%%%%%%%%%%%%%%
truex = xbnd(1):truew(1):xbnd(2);
truey = ybnd(1):truew(2):ybnd(2);
truep = pdf_warpN(truex, truey, m1, m2, data_N, sig0, sigx);
truep = truep / sum(sum(truep)) / prod(truew);

% svgd
dlog_p = @(x)dlog_p_warpN(x, m1, data_N, sig0, sigx);
xy = xy_init;
figure;
for i = 1:9
	xy = svgd(xy, dlog_p, num_iter, eps);
	subplot(3,3,i);
	contour(truex, truey, truep, num_cont);
    hold on;
    scatter(xy(:,1),xy(:,2),ssize,'.');
end

% rsvgd_div
dlog_p = @(x)dlog_p_warpN(x, m1, data_N, sig0, sigx);
gradDet = @(x)gradDet_warpN(x, data_N, sig0, sigx);
xy = xy_init;
figure;
for i = 1:9
	xy = rsvgd_div(xy, dlog_p, gradDet, num_iter, eps);
	subplot(3,3,i);
	contour(truex, truey, truep, num_cont);
    hold on;
    scatter(xy(:,1),xy(:,2),ssize,'.');
end

% rsvgd_nat
dlog_p = @(x)dlog_p_warpN(x, m1, data_N, sig0, sigx);
gradDet = @(x)gradDet_warpN(x, data_N, sig0, sigx);
Ginv = @(x)Ginv_warpN(x, data_N, sigx, sig0);
xy = xy_init;
figure;
for i = 1:9
	xy = rsvgd_nat(xy, dlog_p, gradDet, Ginv, num_iter, eps);
	subplot(3,3,i);
	contour(truex, truey, truep, num_cont);
    hold on;
    scatter(xy(:,1),xy(:,2),ssize,'.');
end

% rsvgd_spe
dlog_p = @(x)dlog_p_warpN(x, m1, data_N, sig0, sigx);
gradDet = @(x)gradDet_warpN(x, data_N, sig0, sigx);
Ginv = @(x)Ginv_warpN(x, data_N, sigx, sig0);
gradGinv = @(x)gradGinv_warpN(x, data_N, sig0, sigx);
xy = xy_init;
figure;
for i = 1:9
	xy = rsvgd_spe(xy, dlog_p, gradDet, Ginv, gradGinv, num_iter, eps);
	subplot(3,3,i);
	contour(truex, truey, truep, num_cont);
    hold on;
    scatter(xy(:,1),xy(:,2),ssize,'.');
end
