function param = thesisFormat(varargin)
%THESISFORMAT Format the current figure for publication.
%
% THESISFORMAT, by itself, makes the current figure suitable for
%   publication in a Thesis or article.
%
% THESISFORMAT(addedSpace[,addedShift]) does not affect the plot 
%   formatting, but is included for compatibility with the legacy 
%   thesisFormat.
%
% THESISFORMAT(Key,Value) sets the keyword argument KEY to VALUE. 
%   The available keyword arguments are:
%     -- Plot formatting --
%     FigSizeInches: [] default (determined automatically); specify to make
%       the figure size abolute
%     FigMarginInches: [.25 .25 .25 .25] default (only used when
%       FigSizeInches is set explicitly)
%     PlotBoxPaddingInches: [left bottom right top] additional space 
%       to put around each plot axes
%     AxesWidthInches
%     AxesLineWidth
%     AxesFontSize
%     AxesLabelFontSize
%     AxesTitleFontSize
%     AxesAspectRatio: width/height ratio of plot box
%     AxesLimits: 'default' (default) or 'Nyquist'
%     LineLineWidth
%     LineMarkerLineWidth: 6.5 default
%     LineMarkerFaceColor: 'auto' default
%     TiledLayoutTitleFontSize
%     TiledLayoutLabelFontSize
%     FontName
%     -- Plot Image Export --
%     SaveName: '' default; if supplied, export figure to the given file
%       (Creates directories if they do not exist!)
%     SaveFlag: true default; if false, figure will not be saved
%     SaveExtension: 'eps png' (default)
%  Note: some arguments are automatically scaled for better fit on the
%  screen (see the scaleFig function defined in this file).
%
% param = THESISFORMAT(...) returns a structure of the parameter values 
%   used to style the figure.
%
% -- Changelog --
% 2023.06.11 | Make compatible w/legacy thesisFormat | Wesley Hileman
% 2023.06.08 | Support subplots and tiledlayout | Wesley Hileman
% 2023.06.07 | Created | Wesley Hileman <whileman@uccs.edu>

% Remove legacy addedShift and addedSpace arguments if present!
args = varargin;
if nargin >= 1
    arg1 = varargin{1};
    if isnumeric(arg1)
        % Legacy addedSpace argument; ignore!
        args = args(2:end);
        % Check for addedShift argument.
        if nargin >= 2
            arg2 = varargin{2};
            if isnumeric(arg2)
                % Legacy addedShift argument; ignore!
                args = args(2:end);
            end
        end
    end
end

% Parse Key-Value parameters.
parser = inputParser;
parser.addParameter('FigSizeInches',[]); % default: automatically size
parser.addParameter('FigMarginInches',[0.25 0.25 0.25 0.25]);
parser.addParameter('AxesWidthInches',3.7);
parser.addParameter('AxesLineWidth',0.7);
parser.addParameter('AxesFontSize',11);
parser.addParameter('AxesLabelFontSize',12);
parser.addParameter('AxesTitleFontSize',14);
parser.addParameter('AxesAspectRatio',2/(sqrt(5)-1)); % width/height, golden ratio
parser.addParameter('AxesLimits','default',@(x)any(strcmpi(x,{'default','Nyquist'})));
parser.addParameter('LineLineWidth',1.8);
parser.addParameter('LineMarkerLineWidth',1);
parser.addParameter('LineMarkerSize',6.5);
parser.addParameter('LineMarkerFaceColor','auto');
parser.addParameter('SaveName','');
parser.addParameter('SaveExtension','-depsc -dpng');
parser.addParameter('SaveFlag',true);
parser.addParameter('TiledLayoutTitleFontSize',18);
parser.addParameter('TiledLayoutLabelFontSize',14);
parser.addParameter('FontName','Times');
parser.addParameter('PlotBoxPaddingInches',[0 0 0 0]);
parser.parse(args{:});
p = parser.Results; % structure of validated parameters
param = p;

% Fetch figure.
fig = gcf;
figUnits = fig.Units;
fig.Units = 'inches';

% Set default axes and line properties. (Applies only to new objects, not
% any existing objects!)
set(fig,'defaultAxesPlotBoxAspectRatioMode','manual');
set(fig,'defaultAxesPlotBoxAspectRatio',[p.AxesAspectRatio 1 1]);
set(fig,'defaultAxesLineWidth',p.AxesLineWidth);
set(fig,'defaultAxesFontName',p.FontName);
set(fig,'defaultAxesFontSize',p.AxesFontSize);
set(fig,'defaultAxesTitleFontSizeMultiplier',p.AxesTitleFontSize/p.AxesFontSize);
set(fig,'defaultAxesTitleFontWeight','normal');
set(fig,'defaultAxesLabelFontSizeMultiplier',p.AxesLabelFontSize/p.AxesFontSize);
set(fig,'defaultAxesXMinorTick','On');
set(fig,'defaultAxesYMinorTick','On');
set(fig,'defaultLineLineWidth',p.LineLineWidth);
set(fig,'defaultLineMarkerFaceColor',p.LineMarkerFaceColor);
set(fig,'defaultTiledLayoutTileSpacing','compact');
set(fig,'defaultTiledLayoutPadding','compact');

tiledLayout = findobj(fig,'Type','tiledlayout');
if length(tiledLayout)==1
    % Tiled layout.

    tUnits = tiledLayout.Units;
    tiledLayout.Units = 'inches';
    
    nrows = tiledLayout.GridSize(1);
    ncols = tiledLayout.GridSize(2);
    scale = getScreenScaleFactor(nrows,ncols);
    preprocessFig(fig,p);
    p = scaleFig(fig,p,scale);
    postprocessFig(fig,p);

    if ~isempty(p.FigSizeInches)
        % User-provided figure size/margins.
        figWidth = p.FigSizeInches(1);
        figHeight = p.FigSizeInches(2);
        figMargins = p.FigMarginInches;
    else
        % Infer figure size from number of tiles.
        figWidth = ncols*p.AxesWidthInches;
        figHeight = nrows*p.AxesWidthInches/p.AxesAspectRatio;
        figMargins = [0 0 0 0];
    end

    % Redefine figure and tiledlayout sizes.
    boxWidth = figWidth-sum(figMargins([1 3]));
    boxHeight = figHeight-sum(figMargins([2 4]));
    tiledLayout.OuterPosition = [figMargins(1) figMargins(2) boxWidth boxHeight];
    fig.Position(3:4) = [figWidth figHeight];
    fig.PaperSize = [figWidth figHeight];

    % Restore tiledlayout units.
    tiledLayout.Units = tUnits;
else
    % Single axis or subplots.

    AX = getFigureAxes(fig);  % matrix for subplots
    [nrows, ncols] = size(AX);
    axVect = AX(:)';
    axUnit = arrayfun(@(ax)ax.Units,AX(:),'UniformOutput',false);
    for ax = axVect
        ax.Units = 'inches';
    end
    scale = getScreenScaleFactor(nrows,ncols);
    preprocessFig(fig,p);
    p = scaleFig(fig,p,scale);
    postprocessFig(fig,p);

    % Redefine figure size to accomodate all subplots.
    if ~isempty(p.FigSizeInches)
        % User-provided figure size.
        figWidth = p.FigSizeInches(1);
        figHeight = p.FigSizeInches(2);
    else
        % Infer size from Axes size and number of subplots.
        [L, B, R, T] = getAxesMargins(AX,p);
        boxWidth = p.AxesWidthInches;
        boxHeight = p.AxesWidthInches/p.AxesAspectRatio;
        
        % Redefine position and size of each axes object.
        boundingHeight = zeros(size(AX));
        boundingWidth = zeros(size(AX));
        for row = 1:nrows
            for col = 1:ncols
                ax = AX(row,col);

                % Calculate width/height needed to contain each subplot. 
                left = max(L(:,col));
                right = max(R(:,col));
                top = max(T(row,:));
                bottom = max(B(row,:));
                boundingHeight(row,col) = boxHeight + top + bottom;
                boundingWidth(row,col) = boxWidth + left + right;

                % Redefine position/size of the present axes.
                posleft = sum(boundingWidth(row,1:(col-1))) + left;
                posbottom = sum(boundingHeight(1:(row-1),col)) + bottom;
                ax.InnerPosition = [posleft posbottom boxWidth boxHeight];
            end
        end
        
        figWidth = max(sum(boundingWidth,2));
        figHeight = max(sum(boundingHeight,1));
    end

    % Redefine figure size.
    fig.Position(3:4) = [figWidth figHeight];
    fig.PaperSize = [figWidth figHeight];

    % Restore axes units.
    for k = 1:length(axVect)
        ax = axVect(k);
        ax.Units = axUnit{k};
    end
end

% Restore figure units.
fig.Units = figUnits;

% Store metadata to figure.
metadata.scale = scale;
metadata.arg = p;
metadata.origin__ = 'thesisFormat';
fig.UserData = metadata;

% - Finish up -
% Handle axes-limit option.
if strcmpi(p.AxesLimits,'Nyquist')
    axObjects = findall(fig,'Type','axes');
    axObjects = axObjects(:)';
    for ax = axObjects
        setAxesNyquist('axes',ax);
    end
end
% Handle export-figure option.
if ~isempty(p.SaveName) && p.SaveFlag
    [filepath, filename] = fileparts(p.SaveName);
    extensions = split(p.SaveExtension);
    if ~isfolder(filepath)
        mkdir(filepath);
    end
    for k = 1:length(extensions)
        print(fullfile(filepath,filename),extensions{k});
    end
end
% Focus on current figure.
figure(fig);

end % thesisFinish()


% Utility functions -------------------------------------------------------

function AX = getFigureAxes(fig)
    % Empty tag to disqualify legend & colorbar
    ax = findobj(fig,'Type','Axes','tag','');

    if length(ax) <= 1
        % Single plot!
        ncols = 1; 
        nrows = 1;
        rows = 1; 
        cols = 1; 
    else
        % Get outer position of each subplot axes in a matrix P;
        % dim1=axes, dim2=[left bottom width height].
        P = getAxesPositionNormalized(ax);
        left = P(:,1);
        bottom = P(:,2);

        % Determine number of rows and columns.
        dx = diff(sort(left))/max(left); 
        dy = diff(sort(bottom))/max(bottom);
        ncols = sum(dx>0.1)+1;
        nrows = sum(dy>0.1)+1;

        % Determine row/column indicies of each axes.
        rows = ceil(bottom*nrows);
        cols = ceil(left*ncols);
    end

    % Construct matrix of axes objects.
    AX = gobjects(nrows,ncols);
    for k = 1:length(ax)
        row = rows(k);
        col = cols(k);
        AX(row,col) = ax(k);
    end
end

function P = getAxesPositionNormalized(ax)
    P = zeros(length(ax),4);
    for k = 1:length(ax)
        a = ax(k);
        aunits = a.Units;
        a.Units = 'normalized';
        P(k,:) = a.Position;
        a.Units = aunits;
    end
end

function h = getTextWidthInches(text)
    units = text.Units;
    text.Units = 'inches';
    h = text.Extent(3);
    text.Units = units;
end

function h = getTextHeightInches(text)
    units = text.Units;
    text.Units = 'inches';
    h = text.Extent(4);
    text.Units = units;
end

function xTickWidth = fudgeXTickHeight(ax,scale)
    % Hack to get the approximate height of the x-axis tick labels.
    lab = ax.XTickLabels(:)';
    if isempty(lab)
        xTickWidth = 0;
        return;
    end
    allTickLabels = strjoin(lab);
    phantom = text(0,0,allTickLabels, ...
        'FontSize',ax.FontSize,'FontName',ax.FontName,'Units','inches');
    xTickWidth = phantom.Extent(4) + scale*0.12; % fudge factor
    delete(phantom);
end

function yTickWidth = fudgeYTickWidth(ax,scale)
    % Hack to get the approximate width of the y-axis tick labels.
    lab = ax.YTickLabels(:)';
    if isempty(lab)
        yTickWidth = 0;
        return;
    end
    allTickLabels = strjoin(lab,'\n');
    phantom = text(0,0,allTickLabels, ...
        'FontSize',ax.FontSize,'FontName',ax.FontName,'Units','inches');
    yTickWidth = phantom.Extent(3) + scale*0.12; % fudge factor
    delete(phantom);
end

function scale = getScreenScaleFactor(nrows,ncols)
    if nrows > 3
        scale = 0.75;
    elseif max(nrows,ncols) > 2
        scale = 1;
    else
        scale = 1.75;
    end
end

function preprocessFig(fig,p)
    %FORMATFIG Apply figure formatting that does not scale.

    % Axes formatting.
    axObjects = findall(fig,'Type','axes');
    axObjects = axObjects(:)';
    for ax = axObjects
        ax.PlotBoxAspectRatioMode = 'manual';
        ax.PlotBoxAspectRatio = [p.AxesAspectRatio 1 1];
        ax.LineWidth = p.AxesLineWidth;
        ax.FontName = p.FontName;
        ax.TitleFontWeight = 'normal';
        ax.TitleFontSizeMultiplier = p.AxesTitleFontSize/p.AxesFontSize;
        ax.LabelFontSizeMultiplier = p.AxesLabelFontSize/p.AxesFontSize;
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        % Axes line formatting.
        lines = findobj(ax,'Type','Line');
        for line = lines(:)'
            line.MarkerFaceColor = p.LineMarkerFaceColor;
        end
    end

    % Tiled layout formatting.
    tiledLayout = findobj(fig,'Type','tiledlayout');
    if length(tiledLayout)==1
        tiledLayout.TileSpacing = 'compact';
        tiledLayout.Padding = 'compact';
        tiledLayout.Title.FontName = p.FontName;
        tiledLayout.XLabel.FontName = p.FontName;
        tiledLayout.YLabel.FontName = p.FontName;
    end
end

function p = scaleFig(fig,p,scale)
    %SCALEFIG Apply figure formatting that scales with the screen.

    % Redefine sizes for scale.
    p.AxesWidthInches = p.AxesWidthInches*scale;
    p.AxesLineWidth = p.AxesLineWidth*scale;
    p.AxesFontSize = p.AxesFontSize*scale;
    p.AxesLabelFontSize = p.AxesLabelFontSize*scale;
    p.AxesTitleFontSize = p.AxesTitleFontSize*scale;
    p.LineLineWidth = p.LineLineWidth*scale;
    p.LineMarkerLineWidth = p.LineMarkerLineWidth*scale;
    p.LineMarkerSize = p.LineMarkerSize*scale;
    p.TiledLayoutTitleFontSize = p.TiledLayoutTitleFontSize*scale;
    p.TiledLayoutLabelFontSize = p.TiledLayoutLabelFontSize*scale;
    
    % Apply scaled sizes to all axes/lines.
    axObjects = findall(fig,'Type','axes');
    axObjects = axObjects(:)';
    for ax = axObjects
        ax.LineWidth = p.AxesLineWidth;
        ax.FontSize = p.AxesFontSize;
        lines = findobj(ax,'Type','Line');
        for line = lines(:)'
            line.LineWidth = p.LineLineWidth;
            line.MarkerSize = p.LineMarkerSize;
        end
    end

    % Apply scaled sizes to tiledlayout.
    tiledLayout = findobj(fig,'Type','tiledlayout');
    if length(tiledLayout)==1
        tiledLayout.Title.FontSize = p.TiledLayoutTitleFontSize;
        tiledLayout.XLabel.FontSize = p.TiledLayoutLabelFontSize;
        tiledLayout.YLabel.FontSize = p.TiledLayoutLabelFontSize;
    end
end

function postprocessFig(fig,p)
    % Apply post-processing operations to a figure. i.e., those that must
    % be performed after all other properties have been updated.

    if p.LineLineWidth ~= p.LineMarkerLineWidth
        % Use different width for lines and markers.
        axObjects = findall(fig,'Type','axes');
        axObjects = axObjects(:)';
        for ax = axObjects
            lines = findall(ax,'Type','Line');
            for k = 1:length(lines)
                line = lines(k);
                if ~strcmpi(line.Marker,'none')
                    if ~strcmpi(line.LineStyle,'none')
                        % MATLAB doesn't support markers and lines with different
                        % widths, so we create two separate line objects with different
                        % line widths.
                        lineonly = copyobj(line,ax);
                        lineonly.Marker = 'none'; % line only!
                        lineonly.LineWidth = p.LineLineWidth;
                        lineonly.HandleVisibility = 'off'; % don't show in legend
                        markersonly = copyobj(line,ax);
                        markersonly.LineStyle = 'none'; % markers only!
                        markersonly.LineWidth = p.LineMarkerLineWidth;
                        markersonly.HandleVisibility = 'off'; % don't show in legend
                        % Hack to make original line invisible without 
                        % "graying out" the legend entry (otherwise, could set
                        % Visible property to 'off').
                        % (https://www.mathworks.com/matlabcentral/answers/353013)
                        line.XData = nan;
                        line.YData = nan;
                    else
                        line.LineWidth = p.LineMarkerLineWidth;
                    end % if
                end % if
            end % for line
        end % for axes
    end % if
end

function [L, B, R, T] = getAxesMargins(AX,p)
    L = zeros(size(AX));
    B = zeros(size(AX));
    R = zeros(size(AX));
    T = zeros(size(AX));
    [nrows, ncols] = size(AX);
    scale = getScreenScaleFactor(nrows,ncols);
    for row = 1:nrows
        for col = 1:ncols
            % Note: row, col indices are relative to the bottom left corner
            % of the figure window! e.g., (1,1) -> bottom left subplot.
            ax = AX(row,col);
    
            % Calculate margins needed to accomodate axis tick marks and labels.
            % Note: 'TightInset' margins are wrong when we mess with the font
            % size (see above), so we infer the top/bottom from the title and x-axis
            % text heights (plus font line height for axis tick marks); we 
            % fudge a bit to get the left margin.
            approxLineHeightInches = 1.5*p.AxesFontSize/72;
            approxXTickHeight = fudgeXTickHeight(ax,scale);
            approxYTickWidth = fudgeYTickWidth(ax,scale);
            leftMargin = getTextWidthInches(ax.YLabel) + ...
                approxYTickWidth + 0.12*scale; % fudge factor
            rightMargin = ax.TightInset(3);
            topMargin = getTextHeightInches(ax.Title) + ...
                0.1*scale; % title should have some margin above
            bottomMargin = getTextHeightInches(ax.XLabel)...
                + approxXTickHeight;
            if leftMargin < approxLineHeightInches
                % Add space for y-axis tick marks.
                leftMargin = approxLineHeightInches;
            end
            if rightMargin < approxLineHeightInches/2
                % Add space for x-axis tick mark.
                rightMargin = approxLineHeightInches/2;
            end
            if topMargin < approxLineHeightInches/2
                % Add space for y-axis tick mark.
                topMargin = approxLineHeightInches/2;
            end
            T(row,col) = topMargin + p.PlotBoxPaddingInches(4);
            B(row,col) = bottomMargin + p.PlotBoxPaddingInches(2);
            L(row,col) = leftMargin + p.PlotBoxPaddingInches(1);
            R(row,col) = rightMargin + p.PlotBoxPaddingInches(3);
        end % for col
    end % for row
end