% Example code for visualizing spatial circulant multi-channel correlation filters
rgb = imread('lena.jpg'); 
img = im2double(rgb); 

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
subplot(1,2,2); all_patches = zeros(N*dsize(1),dsize(2),3); 

% Display all the images we are collecting
h2 = imagesc(all_patches); axis off; axis square; 
title('Concatenation of Sub-Images (X)'); 

% Allocate space for the sub-images and labels
X = zeros(N,N,3); % vectorized sub-images
y = zeros(N,1); % output 

sigma = 5; 
for n = 1:N
    dpn = dP(:,n); 
    xn = circshift(x,[dpn(2),dpn(1)]); % Circular shift 
    X(n,:,:) = reshape(xn,[1,size(xn,1)*size(xn,2),3]); % Store the sub-patch
    y(n) = exp(-dpn'*dpn/sigma); % Store the labels
    all_patches((n-1)*dsize(1)+1:n*dsize(1),:,:) = xn; % Set the image
    set(h1,'CData',xn); % Set the data showing all the patches
    set(h2,'CData',all_patches); % Set the data showing all the patches
    drawnow; 
end

% Shift the y-response so the center is at position (1,1)
% fftshift function does this, or could be done using circshift
yf = fft2(ifftshift(reshape(y,dsize)));

% Display the desired output response
figure(2); clf; 
mesh(dx,dy,reshape(y,dsize)); title('Desired output response'); 

% Reshape to solve
X = reshape(X,[size(X,1),size(X,2)*size(X,3)])'; 

% Finally form the auto-scatter and cross-scatter matrices
S = X*X'; I = eye(3*N); % Setup the S matrix and the identity matrix
figure(3); h = (S + 1e1*I)\(X*y); % Use backslash instead of inv (more stable)
g = reshape(h,[dsize(1),dsize(2),3]); % Reshape the weight vector

% Place your code here for solving for g in the Fourier domain (remember to
% use the yf (the fft of y after being fftshifted) provided. 

% Obtain the final single-channel response using multi-channel weight vectors and image
r = imfilter(img(:,:,1),g(:,:,1)) + imfilter(img(:,:,2),g(:,:,2)) + imfilter(img(:,:,3),g(:,:,3));
imagesc(r); 




