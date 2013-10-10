% Example script to reconstruct a 2 cube phantom using CGLS algorithm.
% The number of rays and voxels are unrealistically small for fast solution
% Use this script to test everything is set-up properly and the code is working

% N. Wadeson

% 11/07/2012

% Solving the algebraic problem: Ax = b.
%
% Parameter Definitions:
% voxels:	A 3-vector giving the number of voxels for reconstruction 
%		in x, y and z respectively (with z pointing upwards). 
%
% iterations:	The number of CGLS iterations to perform.  The solution 
%		vector x is output for each multiple of 10 iterations
%		and/or the maximum iteration.
% 		For example:
%		iterations = 5, x is output at 5i
%		iterations = 15, x is output at 10i and 15i
%		iterations = 26, x is output at 10i, 20i, 26i
%
% output:
% x: 		The solution vector for iterations as above.
% rho:
% eta:

%------------------Parameters--------------------------------------------
voxels = [10 10 5];
iterations = 10;
%------------------------------------------------------------------------

test = 1;
% create simple 3 cube phantom
[b geom] = create_phantom(test);

% reconstruct at very low resolution (to save time)
[x rho eta] = cgls_XTek_single(b, iterations, geom, voxels);

% for one solution output
scrollView(reshape(x,voxels),3,0)

% if there is more than one output of the solution (i.e. iterations > 10)
%  nSols = 2 % set nSols to the number of output solutions
%  for i = 1:nSols
%   scrollView(reshape(x((i-1)*prod(voxels)+1:i*prod(voxels)),voxels),3,0)
%   s = sprintf('solution %d', i);
%   title(s)
%  end
disp('Test complete.  Congratulations!  Everything is working perfectly.') 

