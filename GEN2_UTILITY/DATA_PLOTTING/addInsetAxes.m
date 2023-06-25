function addInsetAxes(xspan,varargin)
%ADDINSETAXES

isLinks = ...
    @(x)isempty(x)||...
    isvector(x)&&all(1<=x&x<=4)||...
    ischar(x)&&strcmpi(x,'auto');
is2Vector = @(x)isvector(x)&&length(x)==2;
parser = inputParser;
parser.addRequired('XSpan',is2Vector);
parser.addRequired('InsetPosition',is2Vector);
parser.addOptional('ReductionFactor',2.8,@(x)isnumeric(x)&&isscalar(x)&&x>1);
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
parser.parse(xspan,varargin{:});
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

% Determine yspan from data.
xspan = arg.XSpan;
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
if ischar(arg.Links) && strcmpi(arg.Links,'auto')
    % Automatically select two links.
    edgeLengths = Inf(size(pointsInset,1),size(pointsBox,1));
    for ki = 1:size(pointsInset,1)
        for kb = 1:size(pointsBox,2)
            pointInset = pointsInset(ki,:);
            pointBox = pointsBox(kb,:);
            edgeLengths(ki,kb) = norm(pointInset-pointBox,2);
        end % for box
    end % for inset
    [~,ind1] = min(edgeLengths(:));
    ind1Inset = mod(ind1-1,size(edgeLengths,1))+1;
    %ind1Box = fix((ind1-1)/size(edgeLengths,2))+1;
    ind1Box = ind1Inset;
    edgeLengths(ind1Inset,:) = Inf;
    edgeLengths(:,ind1Box) = Inf;
    [~,ind2] = min(edgeLengths(:));
    ind2Inset = mod(ind2-1,size(edgeLengths,1))+1;
    %ind2Box = fix((ind2-1)/size(edgeLengths,2))+1;
    ind2Box = ind2Inset;
    indInset = [ind1Inset ind2Inset];
    indBox = [ind1Box ind2Box];
    xedges = [pointsInset(indInset,1) pointsBox(indBox,1)];
    yedges = [pointsInset(indInset,2) pointsBox(indBox,2)];
else
    ind = arg.Links;
    xedges = [pointsInset(ind,1) pointsBox(ind,1)];
    yedges = [pointsInset(ind,2) pointsBox(ind,2)];
end % if

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