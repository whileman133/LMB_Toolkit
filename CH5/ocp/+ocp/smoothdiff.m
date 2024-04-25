function [dxdy_y, yout, dxdy_x, xout, dytrue] = smoothdiff(x, y, dy)
%SMOOTHDIFF Numerically differentiate noisy data using the 
%  "histogram counting" method.
% 
% Given vectors x and y, computes dxdy assuming regular spacing in the x
% samples. Specify the resolution (bin width), which controls the amount 
% of smoothing, using dy.
%
% The ideal (i.e. the underlying noise free) y vector should be monotonic
% (i.e. either increasing or decreasing over its entrire domain).
%
% Use the outputs dxdy_y and yout for viewing dxdy over y.
% Use the outputs dxdy_x and xout for viewing dxdy over x.
% Due to some bins containing zero samples (see explaination below),
% the lengths of dxdy_y,yout and dxdy_x,xout may differ.

% Determine histogram bin edges using the inputs y and dy.
% We force edges(1)=min(y) and edges(end)=max(y) by adjusting dy as needed.
% (dytrue is the actual dy value.)
miny = min(y); maxy = max(y); deltay = maxy - miny;
nbins = ceil(deltay/dy); dytrue = deltay / nbins;
yEdges = linspace(miny, maxy, nbins+1);  % one more edge than bin

% Count occurences of y in each bin. Use the result to approximate dx/dy. 
ny = histcounts(y, yEdges);
dx = ny*(max(x)-min(x))/length(x);  % Change in x across each bin, vector.
dxdy_y = dx/dytrue;

% Some bins may contain zero samples (i.e. ny=0, zero slope.) This is okay 
% when viewing dxdy over y, but when viewing dxdy over x there will be 
% more than one value of dxdy assigned to each x (i.e. duplicate entries 
% in the xout vector due to some dx=0.) For this reason we remove sample 
% points for which dx=0 in the dxdy vs x output.
dxdy_x = dxdy_y(ny~=0);
dx = dx(ny~=0);

% Edges of the x-bins corresponding to the y-bins.
% Order depends on whether y increases or decreases with x.
if y(1) > y(end)
    xEdges = max(x) - [0 cumsum(dx)];
else
    xEdges = min(x) + [0 cumsum(dx)];
end

% Derivative values correspond the center of the bins.
yout = (yEdges(1:end-1)+yEdges(2:end))/2;
xout = (xEdges(1:end-1)+xEdges(2:end))/2;

end