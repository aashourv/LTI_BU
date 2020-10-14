function [A,B,u,C] = LTI_BU_Est(V,varargin)
%[A,B,u,C] = LTI_BU_Est(varargin) estimates the parameters of a first-order
%  multivariate AR model,
%
%      V(:,k) = A*V(:,k-1) + B*u(:,k-1) + noise(C),
%
%  Matrix V (n*t) must contain the time series data, with rows of V 
%  representing variables and columns of V representing observations. This 
%  script provides the least squares estimates of the coefficient matrices 
%  A (n*n) and the noise covariance matrix C (n*n). Input matrix B (n*p) 
%  encodes the spatial profiles of external inputs u (p*t), estimated using 
%  Expectation–Maximization-like algorithm and Lasso regularization from 
%  the model residuals. The intercept vector is zero, which assumes that 
%  the AR process has zero mean. System matrix A is estimated using ARfit 
%  package (https://github.com/tapios/arfit).
%
%  T. Schneider and A. Neumaier, 2001: Algorithm 808: ARfit â€“ A Matlab 
%  package for the estimation of parameters and eigenmodes of multivariate 
%  autoregressive models. ACM Trans. Math. Softw., 27, 58-65.
%
%  Function paramters:
%  'sensInd'        Index of modeled channels/sensors (Default=1:size(V,1)) 
%  'numInp'                        Dimensions of input matrix B (Default=1)
%  'iterations'     Number of Lasso regularization repetitions (Default=25)
%  'Reg_Fact'                           Regularization factor (Default=0.5)
%
%  
%  Arian Ashourvan, University of Pennsylvania,
%  Sérgio Pequito, TU Delft,
%  October 2020 
%  ========================================================================

global numCh
A_flag=0;
nargin = length(varargin);
if nargin > 1
    for i = 1:2:nargin
        switch varargin{i}
            case 'sensInd'
                sensInd = varargin{i+1};
            case 'numInp'
                numInp = varargin{i+1};
            case 'iterations'
                niter = varargin{i+1};
            case 'Reg_Fact'
                Reg_Fact = varargin{i+1};
            case 'A_Linear'
                A_flag=1;
                A = varargin{i+1};
                C=nan;
        end
    end
end

K=size(V,2);
if exist('sensInd') ==0
    sensInd = 1:size(V,1);
end
if exist('numInp') ==0
    numInp = 1;
end
if exist('Reg_Fact') ==0
    Reg_Fact = 0.5;
end
if exist('niter') ==0
    iterations = 25;
end

u = zeros(numInp, K); 
numCh = length(sensInd);
sampleID = [0:K-1]+1;
X = V(sensInd,sampleID);

% Estimate system matrix A if not provided 
if A_flag==0  
    [w, A, C, sbc, fpe, th]=arfit(X', 1, 1,'zero');
end
% Initialize B
B = zeros(size(A,1),size(A,2));
Flag_B_Loop=0;
while Flag_B_Loop==0
    SS=rand(size(B,1));
    SS=SS+SS'/2;
    SS=SS>0.7;
    B=SS;
    B=double(B);
    [~,r] = qr(B);
    colInd = find(abs(diag(r))>1e-7);
    if   length(colInd)<numInp
        display('Loop B');
    else
        
        colInd = colInd(1:numInp);
        B = B(:,colInd);
        Flag_B_Loop=1;
    end
end

if rank(B) < numInp
    error('rank deficient B');
end

% Temporal (u) and spatial (B) Lasso regularization step
for iterInd = 1:niter
    for kInd = 2:K
        yUse(:,kInd) = X(:,kInd) - A*X(:,kInd-1);
        u(:,kInd) = getLassoSoln(B, yUse(:,kInd), Reg_Fact); 
    end
    B= getLassoSoln( u(:,:)',yUse(:,:)', Reg_Fact);
    B=B';
end

function out = getLassoSoln(A, b, lambda)
% method :
% 1: CVX
% 2: lasso MATLAB
% 3: ADMM

method = 3;
switch method
    
    case 1
        cvx_begin quiet
        variable uGG(32)
        minimize(sum_square(A*uGG - b) + lambda*norm(uGG,1))
        cvx_end
        out = uGG;
    case 2
        [out, fInfo] = lasso(A, b, 'Lambda', lambda);
        
        % problem: it doesn't allow to remove the intercept
        
    case 3
        %         downloaded from: https://web.stanford.edu/~boyd/papers/admm/lasso/lasso.html
        % QUIET = 1;
        MAX_ITER = 100;
        ABSTOL   = 1e-4;
        RELTOL   = 1e-2;
        
        [m, n] = size(A);
        Atb = A'*b;
        % lambda = 1;
        rho = 1/lambda;
        alpha = 1;
        
        x = zeros(n,1);
        z = zeros(n,1);
        u = zeros(n,1);
        
        [L, U] = factor(A, rho);
        
        for k = 1:MAX_ITER
            
            % x-update
            q = Atb + rho*(z - u);    % temporary value
            if( m >= n )    % if skinny
                x = U \ (L \ q);
            else            % if fat
                x = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
            end
            
            % z-update with relaxation
            zold = z;
            x_hat = alpha*x + (1 - alpha)*zold;
            z = shrinkage(x_hat + u, lambda/rho);
            
            % u-update
            u = u + (x_hat - z);
            
            % diagnostics, reporting, termination checks
            %  history.objval(k)  = objective(A, b, lambda, x, z);
            
            history.r_norm(k)  = norm(x - z);
            history.s_norm(k)  = norm(-rho*(z - zold));
            
            history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
            history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
            
            %     if ~QUIET
            %         fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            %             history.r_norm(k), history.eps_pri(k), ...
            %             history.s_norm(k), history.eps_dual(k), history.objval(k));
            %     end
            
            if (history.r_norm(k) < history.eps_pri(k) && ...
                    history.s_norm(k) < history.eps_dual(k))
                break;
            end
            
        end
        out = z;
end


function p = objective(A, b, lambda, x, z)

p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) );


function z = shrinkage(x, kappa)

z = max( 0, x - kappa ) - max( 0, -x - kappa );

function [L, U] = factor(A, rho)

[m, n] = size(A);
if ( m >= n )    % if skinny
    L = chol( A'*A + rho*speye(n), 'lower' );
else            % if fat
    L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
end

% force matlab to recognize the upper / lower triangular structure
L = sparse(L);
U = sparse(L');


