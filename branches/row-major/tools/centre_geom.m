function geom = centre_geom(data, geom)
% function to calculate the centre of rotation for given XTek geometry and
% data, and apply this to the geometry structure.

% W. Thompson

% 29/11/2011

% check to see if 2D data, if not then work with the central slice
if geom.dets.nz ~= 1
    [data tmp_geom] = convert2D(data, geom);
else
    tmp_geom = geom;
end

centre = find_centre(reshape(data, geom.dets.ny, length(geom.angles)), tmp_geom);

geom.source.y = geom.source.y + centre;
geom.dets.y = geom.dets.y + centre;
