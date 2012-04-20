% Example script to check everything is working

% W. Thompson

% 03/04/2012

% create simple 3 cube phantom
[b geom] = create_phantom;

% reconstruct at very low resolution (to save time)
[x rho eta] = cgls_XTek_single(b, 10, geom, [60 60 40]);

scrollView(reshape(x,60,60,40),3,0)