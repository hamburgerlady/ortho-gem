# ortho-gem
MATLAB implementations of minimal and optimal algorithms for Orthographic essential matrix estimation

Please cite:

Oskarsson, M. J Math Imaging Vis (2017). https://doi.org/10.1007/s10851-017-0753-1

Contains: 
* test_orthogem.m - testscript for all solvers
* minorthoF_mat.m - minimal solver using 3 point correspondence
* ortho_ransac_mat.m - RANSAC wrapper using minorthoF_mat.m
* ortho_optimal.m - Optimal solver that maximizes the number of inliers
* mlorthoF_cm.m - Least squares solver 
