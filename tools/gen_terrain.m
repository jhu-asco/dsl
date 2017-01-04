clc;

%Basic way to generate image files with the pixel values varying as 
% mixture of gaussians

% user input grid size and dimension
cs = 0.1; %meters
xlb_center = [-15,-15]; % lower bound of cell center
xub_center = [ 15, 15]; % upper bound of cell center

% user input gaussian hills
idx_g=0;
gauss_sigma = zeros([2,2,1]);

idx_g = idx_g+1;
gauss_mu(idx_g,:) = [ 0, 0];
gauss_sigma(1:2,1:2,idx_g) = diag([10,10]);
height(idx_g) = 2; % meters

idx_g = idx_g+1;
gauss_mu(idx_g,:) = [ -15, -9];
gauss_sigma(1:2,1:2,idx_g) = diag([10,10]);
height(idx_g) = 1; % meters

idx_g = idx_g+1;
gauss_mu(idx_g,:) = [ 0, 12];
gauss_sigma(1:2,1:2,idx_g) = diag([4,10]);
height(idx_g) = 1; % meters

%% height channel

xlb = xlb_center -cs/2; %lower bound of the grid
xub = xub_center +cs/2; %lower bound of the grid
ds = xub - xlb;    %dimension size
gs = round(ds/cs); %grid size

n = size(gauss_mu,1);
img = zeros([gs,3]);

for idx_x = 1:gs(1)
   for idx_y = 1:gs(2)
       xy = xlb_center + ([idx_x, idx_y] - 1)*cs;
       X = repmat(xy,[n,1]);
       for idx_g = 1:n
           h = height(idx_g)*exp(-0.5*(xy - gauss_mu(idx_g,:))/gauss_sigma(:,:,idx_g)*(xy - gauss_mu(idx_g,:))');
           img(idx_y,idx_x,1) = img(idx_y,idx_x,1) + h;
       end
   end
end
hmax = max(img(:));
img = img/hmax; %normalize the image
%% traversibility channel

% for idx_x = 1:gs(1)
%    for idx_y = 1:gs(2)
%      img(idx_y,idx_x,2) = img(idx_y,idx_x,1) + h;
%    end
% end

%% display image and write to file
imshow(flipud(img));
imwrite(img,'terrain.ppm');

%% write the .tmap file

hscale = hmax/255;
tscale = 0.1;

fid = fopen('terrain.tmap','w');
fprintf(fid,'# Data that makes a terrain map\n');
fprintf(fid,'#   image: image file location relative to this file(R: height, G: traversibility)\n');
fprintf(fid,'#   cellsize: in meters along x and y direction\n');
fprintf(fid,'#   hscale: height = hscale*pixel value\n');
fprintf(fid,'#   tscale: traversibility = tscale*pixel value\n');
fprintf(fid,'image = terrain.ppm\n');
fprintf(fid,'cellsize= %f,%f\n',cs,cs);
fprintf(fid,'hscale= %f\n', hscale);
fprintf(fid,'tscale= %f\n', tscale);

