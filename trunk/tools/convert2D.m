function [b geom] = convert2D(data, geom)
% function to convert 3D XTek data into 2D, using just (one of) the central
% slice(s).

% W. Thompson

% 29/11/2011

%c_slice = geom.dets.nz/2 + 1;
c_slice = floor(geom.dets.nz/2 + 1);

b = squeeze(data(:,c_slice,:));

geom.dets.z = geom.dets.z(c_slice);

geom.dets.nz = 1;

b = b(:);