% Simple toy problem to show off properties of an FFT

% Randomly generate a 5x5 2D sub-image
D = 5; 
x = randn(D,D); 

% Randomly generate some weights for N possible sub-images
N = D*D; 
y = randn(N,1); 

% Apply the 2D fft to the patch
xf = fft2(x);

% Now obtain the matrix F => equivalent to F*x = fft(x)
F = fftmtx(D); 

% Get the 2D version using the property
% vec(F'*X*F) = kron(F',F')*vec(X)
F2 = kron(conj(F'),conj(F'));

% Output the two results side by side (should be identitical)
fprintf('Result 1: F2*x(:) is equivalent to vec(fft2(x))\n'); 
fprintf('Two complex column vectors should be identitcal\n'); 
[F2*x(:),xf(:)]
fprintf('Press any key to continue.....\n'); 
pause;

% Now generate a covariance matrix S over ALL possible circular shifts
Sxx = zeros(D*D,D*D); 
Sxy = zeros(D,1); 
Cxx = zeros(D*D,D*D); 
for n = 1:D
    for m = 1:D
        xnm = circshift(x,[n,m]);
        xfnm = fft2(xnm); % Take the 2D fft
        Sxx = Sxx + xnm(:)*xnm(:)';
        Cxx = Cxx + xfnm(:)*xfnm(:)'; 
    end
end
Cxx = Cxx/(D*D); % Normalize by D^2 samples

% Finally apply the F2 matrix
Sf = F2*Sxx*F2';

% Calculate the spectrum
sf = xf.*conj(xf); 

% Highlight result 2
fprintf('Result 2: the spectrum of x and diagonal of Cxx are the same\n'); 
[diag(Cxx),sf(:)]
fprintf('Press any key to continue.....\n'); 
pause;

% Highlight result 3
fprintf('Result 3: the spectrum of x and eigenvalues of Sxx are the same\n'); 
[eig(Sxx),sort(sf(:))]
fprintf('Press any key to continue.....\n'); 
pause;


