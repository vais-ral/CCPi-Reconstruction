function x = PBbackproject_single(b, geom, voxels, voxel_size, image_offset)
% function to perform parallel beam back projection for algebraic reconstruction

% interface to C routine to do the back projection

% 06/09/2011

% x = volume data
% b = ray data
% geom = geometry structure array
% voxels, voxel_size, image_offset as defined in Jacobs rays code

x = PBbackproject_single_newgeom_c(voxels,geom.source.x,geom.source.y,geom.source.z,geom.dets.x,geom.dets.y,geom.dets.z,voxel_size,image_offset,b,geom.angles);


% mex PBbackproject_single_c.c backproject_singledata.c -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
