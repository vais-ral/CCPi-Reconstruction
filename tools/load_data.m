function [data geom] = load_data(pathname, filename)
% function to load XTek data and corresponding machine geometry

% pathname: name of path where files are stored
% filename: name of files to load (name of .xtekct file without extension)
% data:     returned data in MATLAB single format
% geom:     returned geometry structure array

% W. Thompson

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
angles = angles(1:nProjections) + initAngle - 90;
angles = angles*pi/180;     % convert to radians
geom.angles = angles;
clear temp angles

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
data = data/65535;
data = -log(data);











