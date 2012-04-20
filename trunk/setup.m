% This just builds the mex files

% Included files:

% backproject_singledata.c         : Performs back projection for single ray
% CBbackproject_single.m           : Essentially just calls CBbackproject_single_newgeom_c.c
% CBbackproject_single_newgeom_c.c : Back projection for entire XTek geometry
% CBproject_single.m               : Essentially just calls CBproject_single_newgeom_c.c
% CBproject_single_newgeom_c.c     : Projection for entire XTek geometry
% centre_geom.m                    : Apply centre of rotation correction (for real data)
% cgls_XTek_single.m               : Reconstruction by CGLS (the main program)
% convert2D.m                      : Converts geometry to 2D slice
% create_phantom.m                 : Create simple cube phantom for testing
% example.m                        : Quick example to check it's working
% find_centre.m                    : Calculates true centre of rotation (for real data)
% load_data.m                      : Load real data from files
% project_singledata.c             : Performs projection for single ray
% scrollView.m                     : Useful for viewing reconstructed volumes

% W. Thompson

% 03/04/2012

mex CBbackproject_single_c.c backproject_singledata.c -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
mex CBproject_single_c.c project_singledata.c -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

