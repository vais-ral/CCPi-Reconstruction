% Example script to check everything is working

addpath tools
addpath mex

% create simple 3 cube phantom
[b geom] = create_PB_phantom;

% reconstruct at very low resolution (to save time)
[x rho eta] = cgls_PB_single(b, 10, geom, [60 60 40]);
%[x rho eta] = cgls_PB_single(b, 10, geom, [459 459 364]);
%[x rho eta] = cgls_PB_single(b, 10, geom, [459 459 2]);

%write_tiff(x, '', 'res', [459 459 2], 16);
%write_tiff(x, '', 'res', [459 459 364], 16);

scrollView(reshape(x,60,60,40),3,0)
