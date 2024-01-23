function yq = linearinterp(x,y,xq)
%LINEARINTERP Linearly interpolate noisy data.
%
% -- Usage --
% yq = LINEARINTERP(x,y,xq)
%
% -- Input --
% x,y = noisy vectors describing a monotonic function: y = f(x)
% xq  = query points at which to evalulate the output
%
% -- Output --
% yq  = vector of function values at the query points: yq = f(xq) 
%
% -- Changelog --
% 2024.01.16 | Clean up, mean y for dup x | Wes Hileman <whileman@uccs.edu>

% Initialize
x = x(:);    % force column vectors
y = y(:);
xq = xq(:);
yq = zeros(size(xq));

% Sort x to ensure it is monotonic for interpolation.
% (noise may corrupt underlying monotonic function).
% Strictly speaking, this alters the function y = f(x),
% but the alterations are within the noise.
[x,IX] = sort(x);
y = y(IX);

% Indicies to three different types of x query points.
lo = xq <= x(1);    % logical indicies to points off the low end
hi = xq >= x(end);  % logical indicies to points off the high end
mid = find(x(1) < xq & xq < x(end));  % indicies to mid points

% 1. xq < x(1): extrapolate off low end of table
if any(lo)
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    yq(lo) = (xq(lo)-x(1))*dy/dx + y(1);
end

% 2. xq > x(end): extrapolate off high end of table
if any(hi)
    dx = x(end) - x(end-1);
    dy = y(end) - y(end-1);
    yq(hi) = (xq(hi)-x(end))*dy/dx + y(end);
end

% 3. x(1) < xq < x(end): linearly interpolate
for k = 1:length(mid)
    ind = mid(k);
    xind = xq(ind);
    if any(x==xind)
        % Exactly at an x-value, no need to interpolate
        % !!! may be more than one if duplicate x values present - use mean
        yq(mid(k)) = mean(y(x==xind));
    else
        x2 = x(x>xind); x2 = x2(1);
        y2 = y(x>xind); y2 = y2(1);
        x1 = x(x<xind); x1 = x1(end);
        y1 = y(x<xind); y1 = y1(end);
        m = (y2-y1)/(x2-x1);
        yq(ind) = y2 + m*(xind-x2);
    end
end

end