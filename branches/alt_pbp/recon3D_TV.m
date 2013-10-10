% Example script to reconstruct in 2D using Total Variation regularisation

% N. Wadeson

% 09/07/2012

% Solving the algebraic problem: Ax = b.
%
% Parameter Definitions:
% Pixels_nY:	The number of detector pixel columns (max 2000)
% Pixels_nZ:    The number of detector pixel rows (max 2000)
%		If only performing a sub-section reconstruction such that
%		not all the rays are required, cut them down here. 
%		The detector array centre remains unchanged. 
%
% voxels:	A 3-vector giving the number of voxels for reconstruction 
%		in x, y and z respectively (with z pointing upwards). 
%		If you are reading in a reconstruction.xtekct file, this
%		value will be overwritten by the number in the file.
%
% alpha: 	The regularisation parameter
% tau: 		Huber function parameter
% opt.bL:	Initial guess of L, the Lipschitz constant
% opt.bmu:	Initial guess of mu, the strong convexity parameter
% opt.k_max: 	The maximum number of iterations to perform
%
% output:
% TV_out: 	The solution vector
% 		For details of other output values see tvreg_upn.m
%
% ***NOTE*** 
% (1) Check which version of load_data you wish to use (below)
% (2) Do you wish to use the beam hardening correction?


%------------------Parameters--------------------------------------------
pixels_nY = 1000;
pixels_nZ = 1000;
voxels = [500 500 1]; 
alpha = 0.005;
tau = 1e-4;
opt.bL = 6930;
opt.bmu = 0.5;
opt.k_max = 10000;
%------------------------------------------------------------------------

% Specify nonnegativity constraints
constraint.type = 2;
constraint.c    = 0*ones(prod(voxels),1);
constraint.d    = 1*ones(prod(voxels),1);

% Options
opt.epsb_rel = 1e-6; % user specified tolerance
opt.K        = 2;
opt.verbose  = 1;
opt.beta     = 0.95;


% TV Reg begins with 5 iterations of cgls
addpath /home/nwadeson/MATERIALS/CGLS_and_TVReg/cgls_XTek/
addpath /home/nwadeson/MATERIALS/CGLS_and_TVReg/cgls_XTek/c
addpath /home/nwadeson/MATERIALS/CGLS_and_TVReg/cgls_XTek/tools

%[data geom] = load_data(pathname, filename);
[data geom] = load_data(pathname, filename, pathname2, filename2);
aY = (2000 - pixels_nY)/2;
aZ = (2000 - pixels_nZ)/2;
geom.dets.z = geom.dets.z(aZ+1:aZ+pixels_nZ);
geom.dets.nz = length(geom.dets.z);
geom.dets.y = geom.dets.y(aY+1:aY+pixels_nY);
geom.dets.ny = length(geom.dets.y);

% Option to overwrite the original data
%data = foam_minus_casing;

data = data(aY+1:aY+pixels_nY, aZ+1:aZ+pixels_nZ,:);

geom = centre_geom(data, geom);
b = data(:);

% Optional beam hardening correction
b = b.^2;

[TV_out fxk_UPN hxk_UPN gxk_UPN fxkl_UPN info_UPN] = tvreg_upn(geom,b,alpha,tau,voxels,constraint,opt);
