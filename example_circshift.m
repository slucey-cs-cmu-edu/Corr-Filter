% Example code for Assignment 4 in 16423
rgb = imread('lena.jpg'); 
img = rgb2gray(rgb); 

% 2D projected points on the book
pts = [248, 292, 248, 292;
       252, 252, 280, 280]; 
 
% Height and width of the template
dsize = [pts(2,4)-pts(2,1)+1,pts(1,2)-pts(1,1)+1]; 

% Set the template points (in order that points appear in image)
tmplt_pts = [0, dsize(2)-1, 0, dsize(2)-1; 
             0, 0, dsize(1)-1, dsize(1)-1]; 

% Step 1. Get the positive example
t = Translation; % Define a Translation object
gnd_p = t.fit(tmplt_pts, pts); % Get the ground-truth warp
x = t.imwarp(img, gnd_p, dsize); % Get the image
figure(1); clf; subplot(1,2,1); colormap('gray'); 
h1 = imagesc(x); axis off; axis image; % Display the image
title('Circular-Shifted Cropped Image'); 

% Step 2. Loop through and gather data
dx = -floor(dsize(2)/2):floor(dsize(2)/2);
dy = -floor(dsize(1)/2):floor(dsize(1)/2);
[dp1,dp2] = meshgrid(dx, dy);
N = length(dp1(:)); % Get the number of warps to obtain
dP = [dp1(:),dp2(:)]'; % Set the 
subplot(1,2,2); all_patches = zeros(N*dsize(1),dsize(2)); 
% Display all the images we are collecting
h2 = imagesc(all_patches); axis off; axis square; 
title('Concatenation of Sub-Images (X)'); 

% Allocate space for the sub-images and labels
X = zeros(N,N); % vectorized sub-images
Xf = zeros(N,N); % vectorized fft2 sub-images 
y = zeros(N,1); % output 

sigma = 5; 
for n = 1:N
    dpn = dP(:,n); 
    xn = circshift(x,[dpn(2),dpn(1)]); % Circular shift 
    X(n,:) = xn(:)'; % Store the sub-patch
    xfn = fft2(xn); Xf(n,:) = xfn(:)'; % Store the fft2 of the patch
    y(n) = exp(-dpn'*dpn/sigma); % Store the labels
    all_patches((n-1)*dsize(1)+1:n*dsize(1),:) = xn; % Set the image
    set(h1,'CData',xn); % Set the data showing all the patches
    set(h2,'CData',all_patches); % Set the data showing all the patches
    drawnow; 
end

% Take the fft of the response

% Display the desired output response
figure(2); clf; 
mesh(dx,dy,reshape(y,dsize)); title('Desired output response'); 

% Finally form the auto-scatter and cross-scatter matrices
Sxx = X'*X/N; Sxy = X'*y/N; I = eye(N);
figure(3); h = (Sxx + 1e3*I)\Sxy; 
r = imfilter(img,reshape(h,dsize)); imagesc(r)




