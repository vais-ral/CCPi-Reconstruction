function b = CBproject_single(x, geom, voxels, voxel_size, image_offset)
% function to perform cone beam back projection for algebraic reconstruction
% from Nikon Xtek data

% interface to C routine to do the back projection

% 06/09/2011

% x = volume data
% b = ray data
% geom = geometry structure array
% voxels, voxel_size, image_offset as defined in Jacobs rays code

b = CBproject_single_newgeom_c(voxels,geom.source.x,geom.source.y,geom.source.z,geom.dets.x,geom.dets.y,geom.dets.z,voxel_size,image_offset,x,geom.angles);


% mex CBproject_single_c.c project_singledata.c -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
