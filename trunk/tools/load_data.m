function [data geom] = load_data(pathname, filename, pathname2, filename2)
% function to load XTek data and corresponding machine geometry

% input:
% pathname: name of path where files are stored
% filename: name of files to load (name of .xtekct file without extension)
%
% optional input of FDK reconstruction.xtekct file for comparison with FDK 
% reconstruction, using the same parameters. 
% pathname2: pathname of .xtekct reconstruction file
% filename2: filename of .xtekct reconstruction file
%
% output:
% data:     returned data in MATLAB single format
% geom:     returned geometry structure array

% W. Thompson (adapted by N. Wadeson to read in a 2nd .xtekct file)

% 14/10/2011

fid = fopen([pathname filename '.xtekct']);  % open .xtekct file for reading

% read parameters into cell array of strings
params = textscan(fid, '%s %s', 'Delimiter', '=', 'HeaderLines', 1, 'CommentStyle', '[');

fclose(fid);            % close file

% extract parameter values into variables
ind = strcmp('SrcToObject', params{1});
SrcToObject = str2double(params{2}(ind));

ind = strcmp('SrcToDetector', params{1});
SrcToDetector = str2double(params{2}(ind));

ind = strcmp('MaskRadius', params{1});
MaskRadius = str2double(params{2}(ind));

ind = strcmp('DetectorPixelsX', params{1});
DetectorPixelsX = str2double(params{2}(ind));

ind = strcmp('DetectorPixelsY', params{1});
DetectorPixelsY = str2double(params{2}(ind));

ind = strcmp('DetectorPixelSizeX', params{1});
DetectorPixelSizeX = str2double(params{2}(ind));

ind = strcmp('DetectorPixelSizeY', params{1});
DetectorPixelSizeY = str2double(params{2}(ind));

ind = strcmp('Projections', params{1});
nProjections = str2double(params{2}(ind));

ind = strcmp('InitialAngle', params{1});
initAngle = str2double(params{2}(ind));

% write source and detector coordinate values into geom structure array
geom.source.x = -SrcToObject;
geom.source.y = 0;
geom.source.z = 0;

geom.dets.x = SrcToDetector - SrcToObject;
geom.dets.y = ((-(DetectorPixelsX - 1)*DetectorPixelSizeX/2):DetectorPixelSizeX:((DetectorPixelsX - 1)*DetectorPixelSizeX/2))';
geom.dets.z = ((-(DetectorPixelsY - 1)*DetectorPixelSizeY/2):DetectorPixelSizeY:((DetectorPixelsY - 1)*DetectorPixelSizeY/2))';
geom.dets.ny = DetectorPixelsX;
geom.dets.nz = DetectorPixelsY;

geom.d_sd = SrcToDetector;

% load in angle data
temp = importdata([pathname '_ctdata.txt'], '\t', 3);
angles = temp.data(:,2);
angles = angles(1:nProjections) - 90;
%angles = angles(1:nProjections) + initAngle - 90;
angles = angles*pi/180;     % convert to radians
geom.angles = angles;
clear temp angles

geom.mask_radius = MaskRadius;

% Set additional parameters if copying an FDK reconstruction
if (nargin > 2)
    fid = fopen([pathname2 filename2 '.xtekct']);  % open .xtekct file for reading
    % read parameters into cell array of strings
    params = textscan(fid, '%s %s', 'Delimiter', '=', 'HeaderLines', 1, 'CommentStyle', '[');

    fclose(fid); % close file

    % extract parameter values into variables
    ind = strcmp('VoxelsX', params{1});
    VoxelsX = str2double(params{2}(ind));

    ind = strcmp('VoxelsY', params{1});
    VoxelsY = str2double(params{2}(ind));

    ind = strcmp('VoxelsZ', params{1});
    VoxelsZ = str2double(params{2}(ind));

    ind = strcmp('VoxelSizeX', params{1});
    VoxelSize = str2double(params{2}(ind));

    ind = strcmp('OffsetX', params{1});
    OffsetX = str2double(params{2}(ind));

    ind = strcmp('OffsetY', params{1});
    OffsetY = str2double(params{2}(ind));

    ind = strcmp('OffsetZ', params{1});
    OffsetZ = str2double(params{2}(ind));

    geom.voxel_size = VoxelSize*[1 1 1];
    geom.nVoxels = [VoxelsX VoxelsY VoxelsZ];
    geom.offset = [OffsetX OffsetY OffsetZ];

end


% load projection data
data = uint16(zeros(DetectorPixelsY, DetectorPixelsX, nProjections));
for i = 1:nProjections
    data(:,:,i) = imread([pathname filename '_' dec2base(i,10,4) '.tif']);
end

% reshape array
data = permute(data,[2 1 3]);

% convert to single precision
data = single(data);

% scale and take -log
data = data/65535; % max value for 16bit integer
data = -log(data);











