function [xkp1 fxkp1 hxkp1 gxkp1 fxkp1l info hxkp1l gxkp1l xlist] = ...
    tvreg_upn(geom,b,alpha,tau,dims,constraint,opt)
% 
% Changes made to original algorithm to calculate A on the fly.  The matrix
% A is too large to store.  Changes incorporate forward and backward projection
% algorithms based on Jacob's ray tracing code.
%
% N. Wadeson 26/06/2012
%
% Solve the problem 
%
% min f(x) = h(x) + g(x) = alpha* TV(x,tau) + 1/2 ||Ax-b||_2^2
% s.t. x in Q
%
% with implicit definitions of h and g.
%
% The function h(x)=TV(X,tau) is a smooth approximation of the TV.
% In fact TV(x,tau) is the Huber functional
% 
% TV(x,tau) = sum_{i,j,l}^{m,n,l} Phi_tau(D_ijl x)
% with
%
%              { ||y||_2 - tau/2   if  ||y||_2>= tau 
% Phi_tau(y) = {
%              { ||y||_2^2 / (2 tau) else
%    
% and D_ijl a finite difference matrix computing the approximated 
% gradient at coordinate (i,j,l).
%
% Input definitions:
% geom:  A struct holding source and detector positions for use in forward
%        and backward projections
%       
% b:     Observed data
%
% alpha: Regularization parameter
%
% tau:   Smoothing parameter of the TV. Suggestion tau =
%        1e-4*norm(X0(:),'inf') where XO is the true X. 
%         
% dims:  dims=[m,n,l] or dims=[m,n]. A vector describing the
%        dimensionality of x, e.g. dims =  [256,256] in the case 
%        of two dimensional problems.
%
% constraint: A struct with the following fields,
%             type = Constraint type
%             c, d = lower and upper bound on x if type=2, where c and d 
%             has the same dimensions as x.
%
%     type == 1
%          Q = R the reals (no constraints)
% 
%     type == 2
%          Q = {x | d_i =>x_i>= c_i, i=0..prod(dims)}
%
% opt:    (Optional) A struct with one or more of the following fields
%          
%         epsb_rel = stopping criteria (default 1e-4)
%
%           The algorithm stops when the iterate y satisfies
%
%              ||G_t(y)||_2/(m*n*l) \leq epsb_rel 
%    
%           where 
%
%               G_t(x) =  1/t (x - P_Q(x - t \nabla f(x)) ) 
%
%           and returns P_Q(x - t \nabla f(x))
%
%         k_max = the maximum number of iterates (default 10000)
% 
%         x0 = initial estimate (default uses 5 cgls iterations using b 
%              and A)
%
%          bL = Initial setting of L_k. Default value is based on an 
%             an estimate of L.
%
%          bmu = Initial setting of mu_k. Default vaule is
%             bmu = min(bL/50,eigmax), where eigmax is the largest
%             eigenvalue of A.
%
% Output definitions
%         Default return values
%         xkp1 = the last iterate
%         fxkp1 = f(xkp1), objective of the last iterate
%         hxkp1 = h(xkp1), the smooth TV of the last iterate
%         gxkp1 = g(xkp1), the fidelity term of the last iterate
%         fxkp1l = a vector containing f(xkp1) for all iterates
%   
%         Optional
%         info = A struct containing additional information on the
%                the behaviour of the algorithm.
%                numFunc = # the objective function is evaluated
%                numGrad = # the gradient function is evaluated
%                numBack = # of backtrackings
%                numRest = # of restarts to reduce \mu_k
%                Lklist = contains the iterates L_k
%                muklist =  constains the iterates mu_k
%                rklist = contains the restart iteration positions.
%
%         hxkp1l = a vector containing h(xkp1) for all iterates
%         gxkp1l = a vector containing g(xkp1) for all iterates 
%         xlist = a matrix containing the iterates xkp1 in all rows.
%                 Make sure only to operate with small x and 
%                 small k_max.
%
%

% Default values.
epsb_rel = 1e-4;
k_max     = 10000;
verbose = 0;

% If options is given as input, use these values
if nargin > 5
    if isfield(opt,'epsb_rel')
        epsb_rel = opt.epsb_rel;
    end
    if isfield(opt,'k_max')
        k_max = opt.k_max;
    end
    if isfield(opt,'verbose')
        verbose = opt.verbose; 
    end


    if isfield(opt,'x0')
        x = opt.x0;   
    else
        addpath /home/nwadeson/MATERIALS/cgls_XTek/
        [x, rho, eta] = cgls_XTek_single(b,5,geom,dims);
    end
else
    addpath /home/nwadeson/MATERIALS/cgls_XTek/
    [x rho eta] = cgls_XTek_single(b,5,geom,dims);
end


if (length(fieldnames(geom)) == 5) % data from file (only one file loaded)
    maxDets_y = 2000;
    geom.voxel_size = ((2*geom.mask_radius)/maxDets_y)*[1 1 1];
    geom.offset = [0.0 0.0 0.0];
elseif (length(fieldnames(geom)) == 6) % using phantom data
    geom.offset = [0.0 0.0 0.0];
elseif (length(fieldnames(geom)) == 8) % data from file (plus 2nd FDK reconstruction file)
    if (dims(3) == 1)
        geom.nVoxels(3) = 1;
        geom.offset(3) = 0.0;

%  % ****COMMENT/UNCOMMENT*************
%  % For TV_reg cut-down only. Uncomment if you wish to reconstruct using
%  % FDK reconstruction file but with a smaller number of z voxels.
%  % Do the same in tvreg_upn.m
%      else
%          geom.nVoxels(3) = 41;
%  % ****COMMENT/UNCOMMENT*************

    end

    dims = geom.nVoxels;
else
    disp('geom structures contains %d fields', length(fieldnames(geom)));
    return
end

image_vol = dims.*geom.voxel_size;
image_offset = -image_vol/2 + geom.offset;


b = double(b);
clear rho eta
prodDims = prod(dims); % nVoxels
lenDims  = length(dims);

if(numel(dims)==2)
    nDd = 8;
else
    nDd = 12;
end

% Intitial settings of bmu and bL
bL = opt.bL;
bmu = opt.bmu;


% Initialize vectors to hold tv and fidelity of the iterates
ghxl = 0;
if nargout > 5
    ghxl=1;
end

% Array to hold iterates in columns
xl=0;
if nargout > 7 & k_max*prodDims<1e7
    xl = 1;
end


if(constraint.type==2)
    [xkp1 fxkp1 hxkp1 gxkp1 fxkp1l k hxkp1l gxkp1l xlist numGrad numBack numFunc numRest Lklist muklist rklist] = tvreg_upn_c(geom.voxel_size, b,alpha,tau,dims,bL, ...
                                                  bmu,epsb_rel,k_max,x, constraint.type, constraint.d, constraint.c,ghxl,xl,verbose, ...
                                                  geom.source.x,geom.source.y,geom.source.z,geom.dets.x,geom.dets.y,geom.dets.z, geom.angles, image_offset);

else
    [xkp1 fxkp1 hxkp1 gxkp1 fxkp1l k hxkp1l gxkp1l xlist numGrad numBack numFunc numRest Lklist muklist rklist] = tvreg_upn_c(geom.voxel_size, b,alpha,tau,dims,bL, ...
                                                  bmu,epsb_rel,k_max,x, constraint.type, 0, 0,ghxl,xl, verbose, ...
                                                  geom.source.x,geom.source.y,geom.source.z,geom.dets.x,geom.dets.y,geom.dets.z, geom.angles, image_offset);
end


if(~ghxl)
    clear gxkp1l;
    clear hxkp1l;
end

if(~xl)
    clear xlist;
else
    xlist =reshape(xlist,prodDims,k_max+1);
end
    
    
% Truncate excess zeros in fxkp1l, xlist, gxkp1l and hxkp1l
fxkp1l = fxkp1l(1:k+1);

if nargout > 8 & k_max*prodDims<1e7
        xlist = xlist(:,1:k+1);
end

if nargout >5
    info.numFunc = numFunc;
    info.numGrad = numGrad;
    info.numBack = numBack;
    info.numRest = numRest;
    info.Lklist = Lklist(1:k+1);
    info.muklist = muklist(1:k+1);
    info.rklist = rklist;
end

if nargout > 6
    gxkp1l      = gxkp1l(1:k+1);
    hxkp1l      = hxkp1l(1:k+1);
end



% If the iteration counter reaches the k_max, the algorithm did not
% converge
if k == k_max 
    disp('Did not find a epsb_rel solution in k_max iterations.')
end