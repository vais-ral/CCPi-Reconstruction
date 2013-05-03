function centre = find_centre(data, geom)
% function to find the centre of rotation of the central slice for XTek
% circle scan data

% centre:   found centre of rotation value to be added to y geometry
% data:     XTek circle scan data
% geom:     geometry structure array

% based on the method described in T. Liu - "Direct central ray
% determination in computed microtomography" April 2009

% W. Thompson

% 21/10/2011

% Update to only include a ray if its opposite is contained within the
% detector region.  Instead of the 2-norm of the image error, find the
% average of the 2-norm of the error for the ray pairs.
% N. Wadeson 15/05/2012
% ***NOTE*** could set largest precision to depend on data size - or quit
% if not enough entries in nonzero.

% set up coordinate grids for testing the fit to data
[alpha s] = meshgrid(geom.angles,geom.dets.y);

% wrap grids in the angular direction to avoid value out of range errors
alpha = [(alpha - 2*pi) alpha (alpha + 2*pi)];
s = repmat(s,1,3);
test_data = double(repmat(data,1,3));


midpoint = 0;   % start search using midpoint at zero
precision = [1 0.1 0.01 0.001 0.0001];   % vector of precision values to search at

figure

for i = 1:length(precision)
    COR = (midpoint - 10*precision(i)):precision(i):(midpoint + 10*precision(i)); % values for centre of rotation
    
    M = zeros(length(COR),1);
    
    for j = 1:length(COR)
        gamma = atan(geom.dets.y / geom.d_sd);   % angle of each ray relative to theoretical central ray

        gamma_c = atan(COR(j) / geom.d_sd); % angle of assumed centre of rotation to central ray

        gamma_i = gamma - gamma_c;

        beta = 2 * gamma_i + pi;

        s2 = geom.d_sd * tan(2 * gamma_c - gamma);

        s2 = repmat(s2, 1, length(geom.angles));

        angles = repmat(geom.angles', geom.dets.ny, 1) + repmat(beta, 1, length(geom.angles));

        test = interp2(alpha, s, test_data, angles, s2, 'linear', 0);
        
        nonzero = find(test > 0);
% BGS surely we want the number of non-zero values for the average,
% not the sum of their positions
        M(j) = sum((test(nonzero) - data(nonzero)).^2)*(1/length(nonzero));
        %M(j) = sum((test(:) - data(:)).^2);
    end
    
    plot(COR * geom.source.x / geom.d_sd, M, '.') % plot to see if results are sensible
    title(['Precision: ' num2str(precision(i))])
    drawnow
    
    [minM indM] = min(M);   % minimum value and index

    fprintf('Precision %1.4f: COR = %1.4f, M = %4.4f\n', precision(i), COR(indM) * geom.source.x / geom.d_sd, minM)
    
    midpoint = COR(indM);   % set midpoint for next search
end

% transform centre to required value
centre = midpoint * geom.source.x / geom.d_sd;
