function insetAx = addInset(varargin)
%ADDINSET Place inset detail inside of the current axes.
%
% Call this function after thesisFormat to place a "zoomed in" detail view
% of the current axes inside of the current axes.
%
% -- Usage --
% ADDINSET(xspan,xyinset) creates an inset view of the current axes over
%   the x-data interval specified by the 2-vector XSPAN. The lower left
%   corner of the inset axes is placed at the data-point specified by the
%   2-vector XYINSET. This function creates the inset by first cloning the
%   current axes, then setting the x-limits to XSPAN, the y-limits such 
%   that all plotted data are visible, and shrinking the size of cloned 
%   axes by a size-reduction factor (default is 2.8).
%
% ADDINSET(xspan,xyinset,red) also specifies the size-reduction 
%   factor, RED, of the inset axes relative to the original axes. 
%   Default red=2.8 (that is, the inset axes is 2.8 times smaller than
%   the original axes).
%
% -- Changelog --
% 2023.06.25 | Created | Wesley Hileman <whileman@uccs.edu>

isLinks = ...
    @(x)isempty(x)||...
    isvector(x)&&all(1<=x&x<=4)||...
    ischar(x)&&strcmpi(x,'auto');
is2Vector = @(x)isvector(x)&&length(x)==2;
parser = inputParser;
parser.addRequired('XSpan',is2Vector);
parser.addRequired('InsetPosition',is2Vector);
parser.addOptional('ReductionFactor',2.8,@(x)isnumeric(x)&&isscalar(x)&&x>1);
parser.addParameter('YSpan',[],is2Vector);
parser.addParameter('Links','auto',isLinks);
parser.addParameter('Axes',gca(),@(x)isgraphics(x,'axes'));
parser.addParameter('ScaleFactor',1.75);
parser.addParameter('ProjectionLineColor',[[0 0 0]/255 1]);
parser.addParameter('ProjectionLineWidth',1/1.75);
parser.addParameter('ProjectionLineStyle','-.');
parser.addParameter('BoxEdgeColor',[[0 0 0]/255 1]);
parser.addParameter('BoxLineWidth',1/1.75);
parser.addParameter('BoxLineStyle','-');
parser.addParameter('InsetAxisRulerColor',[0 0 0]/255);
parser.parse(varargin{:});
arg = parser.Results; % structure of validated arguments

% Fetch thesisFormat metadata.
fig = arg.Axes.Parent;
metadata = fig.UserData;
if ~isstruct(metadata)||...
   ~isfield(metadata,'origin__')||...
   ~strcmp(metadata.origin__,'thesisFormat')
    error('Must call thesisFormat() before addInsetAxes().');
end
scale = metadata.scale;  % screen scale factor

% Scale line widths.
arg.ProjectionLineWidth = arg.ProjectionLineWidth*scale;
arg.BoxLineWidth = arg.BoxLineWidth*scale;

% Determine inset position in relative figure units.
% First, get position relative to axes.
xl = arg.Axes.XLim(1);
yl = arg.Axes.YLim(1);
dx = diff(arg.Axes.XLim);
dy = diff(arg.Axes.YLim);
x0 = (arg.InsetPosition(1)-xl)/dx;
y0 = (arg.InsetPosition(2)-yl)/dy;
% Next, get position relative to figure window.
xl = arg.Axes.Position(1);
yl = arg.Axes.Position(2);
dx = arg.Axes.Position(3);
dy = arg.Axes.Position(4);
x0 = x0*dx+xl;
y0 = y0*dy+yl;
w = arg.Axes.Position(3)/arg.ReductionFactor;
h = arg.Axes.Position(4)/arg.ReductionFactor;
InsetPosition = [x0 y0 w h];

xspan = arg.XSpan;
yspan = arg.YSpan;
if isempty(yspan)
    % Determine yspan from data.
    yspan = [+inf -inf];
    for line = findall(arg.Axes,'Type','line')'
        span = xspan(1)<=line.XData&line.XData<=xspan(2); % logical indicies to x-span
        ymin = min(line.YData(span));
        ymax = max(line.YData(span));
        if ymin<yspan(1)
            yspan(1) = ymin;
        end
        if ymax>yspan(2)
            yspan(2) = ymax;
        end
    end % for
end

% Create inset axes.
insetAx = copyobj(arg.Axes,fig);
insetAx.XLim = xspan;
insetAx.YLim = yspan;
insetAx.Position = InsetPosition;
insetAx.XLabel = [];
insetAx.YLabel = [];
insetAx.Title = [];
insetAx.XAxis.Color = arg.InsetAxisRulerColor;
insetAx.YAxis.Color = arg.InsetAxisRulerColor;
insetAx.FontSize = insetAx.FontSize/arg.ScaleFactor;
for line = findall(insetAx,'Type','line')'
    line.LineWidth = line.LineWidth/arg.ScaleFactor;
    line.MarkerSize = line.MarkerSize/arg.ScaleFactor;
end

% Calculate bounds of inset axes in data units.
xl = arg.Axes.XLim;
yl = arg.Axes.YLim;
po = arg.Axes.Position;
pi = insetAx.Position;
x1 = (pi(1)-po(1))/po(3)*diff(xl)+xl(1);
x2 = (pi(1)+pi(3)-po(1))/po(3)*diff(xl)+xl(1);
y1 = (pi(2)-po(2))/po(4)*diff(yl)+yl(1);
y2 = (pi(2)+pi(4)-po(2))/po(4)*diff(yl)+yl(1);

% Determine start/end points for projection lines.
pointsInset = [
    x1 y1; 
    x2 y1; 
    x2 y2; 
    x1 y2
];
pointsBox = [
    xspan(1) yspan(1); 
    xspan(2) yspan(1); 
    xspan(2) yspan(2); 
    xspan(1) yspan(2);
];
midpointInset = [mean([x1 x2]) mean([y1 y2])];
quad2boxlinks.N = [4 3];
quad2boxlinks.NE = [2 4];
quad2boxlinks.E = [2 3];
quad2boxlinks.SE = [1 3];
quad2boxlinks.S = [1 2];
quad2boxlinks.SW = [2 4];
quad2boxlinks.W = [1 4];
quad2boxlinks.NW = [1 3];
quad2boxlinks.C = [1 2 3 4];  % center quadrant
quad2insetlinks.N = [1 2];
quad2insetlinks.NE = [2 4];
quad2insetlinks.E = [1 4];
quad2insetlinks.SE = [1 3];
quad2insetlinks.S = [4 3];
quad2insetlinks.SW = [2 4];
quad2insetlinks.W = [2 3];
quad2insetlinks.NW = [1 3];
quad2insetlinks.C = [1 2 3 4];  % center quadrant
if ischar(arg.Links) && strcmpi(arg.Links,'auto')
    % Automatically select links.
    % First, determine quadrant of inset relative to box.
    if midpointInset(1)>xspan(2)
        if midpointInset(2)>yspan(2)
            quad = 'NE';
        elseif midpointInset(2)<yspan(1)
            quad = 'SE';
        else
            quad = 'E';
        end
    elseif midpointInset(1)<xspan(1)
        if midpointInset(2)>yspan(2)
            quad = 'NW';
        elseif midpointInset(2)<yspan(1)
            quad = 'SW';
        else
            quad = 'W';
        end
    else
        if midpointInset(2)>yspan(2)
            quad = 'N';
        elseif midpointInset(2)<yspan(1)
            quad = 'S';
        else
            quad = 'C';
        end
    end % if
    % Now, lookup links based on quadrant.
    linksBox = quad2boxlinks.(quad);
    linksInset = quad2insetlinks.(quad);
else
    linksBox = arg.links;
    linksInset = arg.links;
end % if
xedges = [pointsInset(linksInset,1) pointsBox(linksBox,1)];
yedges = [pointsInset(linksInset,2) pointsBox(linksBox,2)];

% Create box around data and projection lines between box and inset axes.
axes(arg.Axes);
arg.Axes.XLim = arg.Axes.XLim;  % prevent automatic rescaling of x-axes
arg.Axes.YLim = arg.Axes.YLim;  % prevent automatic rescaling of y-axes
rectangle( ...
    'Position',[xspan(1) yspan(1) diff(xspan) diff(yspan)], ...
    'EdgeColor',arg.BoxEdgeColor, ...
    'LineWidth',arg.BoxLineWidth, ...
    'LineStyle',arg.BoxLineStyle ...
); hold on;
for k = 1:size(xedges,1)
    plot(xedges(k,:),yedges(k,:), ...
        'Color',arg.ProjectionLineColor, ...
        'LineWidth',arg.ProjectionLineWidth, ...
        'LineStyle',arg.ProjectionLineStyle,...
        'HandleVisibility','off' ...
    );
end

% Bring inset axes to front, keep focus on original axes.
uistack(arg.Axes,'bottom');
uistack(insetAx,'top');

end