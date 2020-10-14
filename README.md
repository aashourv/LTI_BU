This script estimates the internal system and unknown external inputs parameters of a linear time-invariant dynamical system from the user-provided time series (i.e., outputs of the LTI system), using the algorithm proposed in Ashourvan et al. 2020 ("A dynamical systems framework to uncover the drivers of large-scale cortical activity"). Specifically, LTI_BU_Est.m function estimates parameters of a first-order multivariate AR model: 

  V(:,k) = AV(:,k-1) + Bu(:,k-1) + noise(C),

Matrix V (n x t) must contain the time series data, with rows of V representing variables and columns of V representing observations. This script provides the least squares estimates of the coefficient matrices A (n x n) and the noise covariance matrix C (n x n). Input matrix B (n x p) encodes the spatial profiles of external inputs u (p x t), estimated using Expectation–Maximization-like algorithm and Lasso regularization from the model residuals. The intercept vector is zero, which assumes that the AR process has zero mean. 

# Function parameters:
 'sensInd' Index of modeled channels/sensors (Default=1:size(V,1)) 
 
 'numInp' Dimensions of input matrix B (Default=1)
 
 'iterations' Number of Lasso regularization repetitions (Default=25)
 
 'Reg_Fact' Regularization factor (Default=0.5)

# Function dependencies: 
 System matrix A is estimated using ARfit package:(https://github.com/tapios/arfit)
T. Schneider and A. Neumaier, 2001: Algorithm 808: ARfit A Matlab package for the estimation of parameters and eigenmodes of multivariate autoregressive models. ACM Trans. Math. Softw., 27, 58-65.

Please contact Arian Ashourvan (ashourv@seas.upenn.edu) with any questions regarding this code.
  
Arian Ashourvan, University of Pennsylvania
 
Sérgio Pequito, TU Delft
 
October 2020 


