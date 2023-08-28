%%	Dynamic Mode Decommposition to 2D IR spectra
%   Use this script to evaluate FFCF from matlab 2D IR spectra file by 
%   dynamic mode decomposition.
%   Modified Aug. 22nd, 2023   Cong Xu, Baiz group

%% Input 
% The input data should be provided in a three-dimensional matrix stored in
% a .mat file format. In this matrix, the x-axis corresponds to the 
% excitation frequency, while the y-axis represents the detection 
% frequency. The spectra within the matrix should be arranged in ascending 
% order based on the t2 parameter, which is aligned along the z-direction.

% If a datastore or 2D IR figures are used, please rewrite the following 
% part to make sure a three-dimensional matrix is loaded in the varible 
% 'spectra'.

close all; clc; clear;
% for example
filename='Example_OnePeakSpectra_2p5';
savepath=pwd;
load([filename '.mat']);

% Please change the variable name on the right to make sure 'spectra' is 
% the loaded three dimesional data matrix like described.
spectra=spectra(:,:,:);

% UserInput
% waiting time step of spectra (ps)
dt = .1; 
% rank of truncation (see Discussion in paper)
r = 4; 
% targeted Input Pixels
targetedPixel=256; 
% the number of displayed original and reconstructed spectra
ImageToDisplay=55; 
% number of Augmented spectra (see equation 7)
s=1;

%% Output
% Phi   -> DMD projected eigenmodes (see equation 5)
% omega -> continuous-time eigenvalues of DMD modes (see step 6 in Algorithm 1)
% ReconstructedImage -> Reconstructed spectra by DMD
% b -> amplitudes of DMD modes

%% data processing

% Data Interpolation
resizedFactor=targetedPixel/size(spectra,1);
ResizedSpectra = imresize(spectra, resizedFactor);
[PlotX,PlotY]=meshgrid(freqAx(1:1/resizedFactor:end),freqAx(1:1/resizedFactor:end));
tToRead=size(ResizedSpectra,3);

% Reshape data from 2D to 1D
Ori_X=reshape(ResizedSpectra,targetedPixel*targetedPixel,tToRead);
X=[];
for i=1:s
    X=[X;Ori_X(:,i:end-s+i)];
end 

%% DMD
% SVD
X1= X(:, 1:end-1);
X2 = X(:, 2:end);
[U, S, V] = svd(X1, 'econ'); % equation 3

% matrix size truncation to rank R (step 2)
r = min(r, size(U,2));
U_r = U(:, 1:r);
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);

% Propogation matrix A (step 3)
PropogationMatrix=U_r' * X2 * V_r / S_r;

% eigenvalues μ_i and eigenvectors w_i (step 4)
[W_r, D] = eig(PropogationMatrix);

% low-rank DMD modes Phi (step 5)
Phi = X2 * V_r / S_r * W_r; 

% discrete-time eigenvalues
lambda = diag(D); 

% continuous-time eigenvalues
omega = log(lambda)/dt; 

% find purely decay modes
decayModeIndex=find(imag(omega)==0);

% Compute DMD mode amplitudes b
x1 = X1(:, 1);
b = Phi\x1;

% DMD spectra reconstruction
Tsteps = size(X1, 2); % Tsteps = m - 1
time_dynamics = zeros(r, Tsteps);
t = (0:Tsteps-1)*dt; % time vector
for i = 1:Tsteps
    time_dynamics(:,i) = (b.*exp(omega*t(i)));
end

Xdmd = Phi(:,decayModeIndex) * time_dynamics(decayModeIndex,:);

%% Plot reconstructed Data
figure();
subplot(1,2,1);
ReconstructedImage=reshape(Xdmd(:,ImageToDisplay),targetedPixel,targetedPixel,[]);
im1=imagesc(ResizedSpectra(:,:,ImageToDisplay));
title('Original Image');
subplot(1,2,2);
im2=imagesc(real(ReconstructedImage(:,:,1)));
title('DMD Reconstrcuted Image');

save([savepath '/' filename '_R=' num2str(r)],'omega','r','b','decayModeIndex');
