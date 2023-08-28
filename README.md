# Dynamics-Mode-Decomposition analysis

Code for the paper [Cutting Through the Noise: Extracting Dynamics from Ultrafast Spectra using Dynamic Mode Decomposition]() by Cong Xu, and Carlos R. Baiz.
Use this script to evaluate FFCF from matlab 2D IR spectra file by dynamic mode decomposition.

# Input
* `dt` -> Waiting time step between adjecent spectra (ps)
* `r` -> Rank of truncation (see Discussion in paper)
* `targetedPixel` -> Targeted input pixels
* `ImageToDisplay` -> The number of displayed original and reconstructed spectra
* `s` -> Number of Augmented spectra (see equation 7)

# Output
* `Phi`   -> DMD projected eigenmodes (see equation 5)
* `omega` -> Continuous-time eigenvalues of DMD modes (see step 6 in Algorithm 1)
* `ReconstructedImage` -> Reconstructed spectra by DMD
* `b` -> amplitudes of DMD modes

