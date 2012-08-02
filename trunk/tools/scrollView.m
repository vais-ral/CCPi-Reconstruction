function hh = scrollView(vol, dim, limits)
% handle = scrollView(vold, dim, limits)
% displays slices through a 3d array, with a 
% slider to scroll through the missing dimension
% 
% scrollView(vol) scrolls thruogh the last dimension
% scrollViwe(vol, n) scrolls through dimension n
% scrollViwe(vol, n, limits) sets the colorlimits on each
% slice. use limits = 0 for max/min (vol(:))
  
% by David Szotten
% The University of Manchester 2009
%
% $Name:  $
% $Revision: 1.5 $
% $Date: 2009/02/13 12:13:31 $

app=[];

app.size = size(vol);

if min(app.size) == 1
	fprintf('Wrong size: ');
	disp(app.size);
	return
end

if ~isreal(vol)
	error('Only real values may displayed!')
end

if nargin > 2
	if limits == 0
		app.limits = [min(vol(:)), max(vol(:))];
		if app.limits(1) == app.limits(2)
			app.limits(1) = app.limits(1) -eps;
			app.limits(2) = app.limits(1) +eps;
		end
	else
		app.limits = limits;
	end
else
	app.limits = [];
end

if nargin < 2
	app.dim = 3;
else
	app.dim = dim;
end

maxval = size(vol, app.dim);

if nargout > 0
	hh = figure;
else
	figure
end
%imagesc( reshape( vol(1000:1013:1013*1013*601), 1013, [])')
colormap jet(1024)

app.axis = gca;
app.figure = gcf;

% Generate constants for use in uicontrol initialization
pos=get(app.axis,'position');
Newpos=[pos(1) pos(2)-0.1 pos(3) 0.05];
% % This will create a slider which is just underneath the axis
% % but still leaves room for the axis labels above the slider

%S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
%S=['title(get(gcbo,''value'')); imagesc( reshape( m.data(floor(get(gcbo,''value''))+1:1013:1013*1013*601), 1013, [])'')'];
% % Setting up callback string to modify XLim of axis (gca)
% % based on the position of the slider (gcbo)

% Creating Uicontrol
h=uicontrol('style','slider',...
    'units','normalized','position',Newpos,...
    'callback',@sfScroll,'min',1, 'max', maxval, 'value', floor(maxval/2), 'sliderstep', [1/maxval 10/maxval] );

%imagesc(zeros(50));
sfScroll( h, [] );
%%%%%%%%%%%%%%%%%%%%%% subfunctions %%%%%%%%%%%%%
% 
function sfScroll(varargin)

handle = varargin{1};
%scrollVal = floor(get(gcbo,'value'));
scrollVal = floor(get(handle,'value'));

	switch app.dim
		case 1
			sfImagescWrapper( vol(scrollVal, :, :), app.limits )
		case 2
			sfImagescWrapper( vol(:, scrollVal, :), app.limits )
		case 3
			sfImagescWrapper( vol(:, :, scrollVal), app.limits )
	end

	title(scrollVal);
	
end


function sfImagescWrapper(image, limits)
	
	if ( length(limits) >= 2 && limits(1) ~= limits(2))
		imagesc( squeeze(image), limits)
% 		hIm = get(app.axis, 'children');
% 		if isempty(hIm)
% 			imagesc( squeeze(image), limits)
% 		end
% 		set(hIm, 'cdata', squeeze(image) );
% 		set(app.axis, 'clim', limits);
	else
		imagesc( squeeze(image))
% 		hIm = get(app.axis, 'children');
% 		if isempty(hIm)
% 			imagesc( squeeze(image))
% 		end
% 		set(hIm, 'cdata', squeeze(image) );

	end

	axis image
	colorbar
end



end
