clear all

% Nflat
fid = fopen('data_Nflat_30_0_10.txt', 'r');
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

M = 50; eps = 0.05; % h = 2;
% xbnd = [-10,10]; ybnd = [5, 41];
xbnd = [-13,20]; ybnd = [1, 20];
truew = [0.2, 0.1];
% xy_init = bsxfun(@plus, randn(M,2)*1e-1, [5,40]);
xy_init = bsxfun(@plus, randn(M,2)*1e-1, [15,2]);
num_iter = 100; num_cont = 10; ssize = 150;

%%%%%%%%%%%%%%
truex = xbnd(1):truew(1):xbnd(2);
truey = ybnd(1):truew(2):ybnd(2);
truep = pdf_Nflat(truex, truey, m1, m2, data_N);
truep = truep / sum(sum(truep)) / prod(truew);

% svgd
dlog_p = @(x)dlog_p_Nflat(x, m1, m2, data_N);
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
dlog_p = @(x)dlog_p_Nflat(x, m1, m2, data_N);
gradDet = @gradDet_Nflat;
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
dlog_p = @(x)dlog_p_Nflat(x, m1, m2, data_N);
gradDet = @gradDet_Nflat;
Ginv = @(x)Ginv_Nflat(x, data_N);
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
dlog_p = @(x)dlog_p_Nflat(x, m1, m2, data_N);
gradDet = @gradDet_Nflat;
Ginv = @(x)Ginv_Nflat(x, data_N);
gradGinv = @(x)gradGinv_Nflat(x, data_N);
xy = xy_init;
figure;
for i = 1:9
	xy = rsvgd_spe(xy, dlog_p, gradDet, Ginv, gradGinv, num_iter, eps);
	subplot(3,3,i);
	contour(truex, truey, truep, num_cont);
    hold on;
    scatter(xy(:,1),xy(:,2),ssize,'.');
end
