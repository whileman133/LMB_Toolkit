function thesisFormat3d(addedSpace,addedShift)
% NOTESFORMAT Set up current figure with format for course notes,
% works for plots having multiple subplots
%
% Call this routine after plotting.  Optional input argument = 
% extra space on [left bottom right top] to leave around *each* frame
% Optional second input argument = space between horizontal and 
% vertical subplots
if nargin < 1,
  %addedSpace = [0 0 0 0];
  addedSpace = [0.2 0.1 0.3 0]; % Ryan Jobman
end
if nargin < 2,
  addedShift = [14/72 14/72];
end

% determine number of subplot rows and columns
ax = findobj(gcf,'type','axes','tag',''); % empty tag to disqualify
if length(ax) == 1,                       % legend & colorbar
%   pos = get(ax,'position');
  rows = 1; cols = 1; ncols = 1; nrows = 1;
else
  pos = cell2mat(get(ax,'position'));
  dx = diff(sort(pos(:,1)))/max(pos(:,1)); ncols = sum(dx>0.1)+1;
  dy = diff(sort(pos(:,2)))/max(pos(:,2)); nrows = sum(dy>0.1)+1;
  rows = ceil(pos(:,2)*nrows);
  cols = ceil(pos(:,1)*ncols);
end
% pos
% xlocs = round(pos(:,1)*20)/20;
% ylocs = round(pos(:,2)*20)/20;
% xlocs = unique(pos(:,1))
% ylocs = unique(pos(:,2))
% % ncols = numel(xlocs) % same X positions
% % nrows = numel(ylocs) % same Y positions
% stop

if max(nrows,ncols)>2,
  SCREENSCALE = 1;
else
  SCREENSCALE = 1.75;
end
if nrows > 3,
  SCREENSCALE = 0.75;
end
addedSpace = addedSpace*SCREENSCALE;

% fix axes limits so labels don't shift around
for k = 1:length(ax),
  xlim(ax(k),get(ax(k),'xlim')); 
  ylim(ax(k),get(ax(k),'ylim'));
  zlim(ax(k),get(ax(k),'zlim'));
end

% Note that the final figure is scaled down from this...
THINLINEWIDTH = 0.7*SCREENSCALE;
THICKLINEWIDTH = 1.8*SCREENSCALE;
AXISFONTSIZE   = 11*SCREENSCALE;
XLABELFONTSIZE = 12*SCREENSCALE;
YLABELFONTSIZE = 12*SCREENSCALE;
TITLEFONTSIZE  = 14*SCREENSCALE;
AWIDTH = 3.7*SCREENSCALE; % axis width, in inches
AHEIGHT = AWIDTH*(sqrt(5)-1)/2; % axis height, using golden ratio

% Set default drawing properties
set(gcf,'DefaultAxesLineWidth',THINLINEWIDTH);
set(gcf,'DefaultLineLineWidth',THICKLINEWIDTH);
for k = 1:length(ax),
  set(gcf,'CurrentAxes',ax(k));
  set(gca,'LineWidth',THINLINEWIDTH);
  set(gca,'FontSize',AXISFONTSIZE,'FontName','Times','FontWeight','normal');
  set(gca,'XMinorTick','On');
  set(gca,'YMinorTick','On');
  set(gca,'ZMinorTick','On');
  set(get(gca,'XLabel'),'FontSize',XLABELFONTSIZE,'FontName','Times','FontWeight','normal');
  set(get(gca,'YLabel'),'FontSize',YLABELFONTSIZE,'FontName','Times','FontWeight','normal');
  set(get(gca,'ZLabel'),'FontSize',YLABELFONTSIZE,'FontName','Times','FontWeight','normal');
  set(get(gca,'Title'),'FontSize',TITLEFONTSIZE,'FontName','Times','FontWeight','normal');
  set(findobj(gca,'Type','Line'),'LineWidth',THICKLINEWIDTH);
end

% Save present system of units
funits=get(gcf,'units');
punits=get(gcf,'paperunits');
set(gcf,'units','inches');
set(gcf,'paperunits','inches');
aunits = cell(size(ax));
for k = 1:length(ax),
  set(gcf,'CurrentAxes',ax(k));
  aunits{k}=get(gca,'units');
  set(gca,'units','inches')
end

% Fix position of axes
padBottom = 1.3*(AXISFONTSIZE + XLABELFONTSIZE)/72 + addedSpace(2);
padTop = 1.2*TITLEFONTSIZE/72 + addedSpace(4);
padRight = 0.08*SCREENSCALE + addedSpace(3);
padLeft = 1.1*YLABELFONTSIZE/72 + 0.25*SCREENSCALE + addedSpace(1);
% and, for multiple sub-plots horizontally or vertically...
padHoriz = addedShift(1)*SCREENSCALE;
padVert  = addedShift(2)*SCREENSCALE;

FWIDTH = AWIDTH + padLeft + padRight;
FHEIGHT = AHEIGHT + padTop + padBottom;
a = get(gcf,'position');
set(gcf,'position',[a(1) a(2) FWIDTH*ncols+padHoriz*(ncols-1),...
                              FHEIGHT*nrows+padVert*(nrows-1)]);
set(gcf,'paperposition',[0 0 FWIDTH*ncols+padHoriz*(ncols-1),...
                             FHEIGHT*nrows]+padVert*(nrows-1));
for k = 1:length(ax),
  set(gcf,'CurrentAxes',ax(k));
  thiscol = cols(k) - 1; % 0..ncols-1
  thisrow = rows(k) - 1; % 0..nrows-1
  posnRect = [FWIDTH*thiscol+padLeft+thiscol*padHoriz,...
                      FHEIGHT*thisrow+padBottom+thisrow*padVert,...
                      AWIDTH AHEIGHT];
  set(gca,'position',posnRect);
end

% Fix labels
for k = 1:length(ax),
  set(gcf,'CurrentAxes',ax(k));

  xL = get(gca,'xlim'); yL = get(gca,'ylim'); %zL = get(gca,'zlim');
  t1 = text(xL(1),yL(1),'X','verticalalignment','middle','horizontalalignment','center','fontsize',30);
  set(t1,'units','inches'); t1p =get(t1,'position');
  t2 = text(xL(2),yL(1),'Y','verticalalignment','middle','horizontalalignment','center','fontsize',30);
  set(t2,'units','inches'); t2p = get(t2,'position');
  t3 = text(xL(1),yL(2),'Z','verticalalignment','middle','horizontalalignment','center','fontsize',30);
  set(t3,'units','inches'); t3p = get(t3,'position');
  t4 = text(xL(2),yL(2),'W','verticalalignment','middle','horizontalalignment','center','fontsize',30);
  set(t4,'units','inches'); t4p = get(t4,'position');
  delete(t1); delete(t2); delete(t3); delete(t4);

%   [posnx,posny] = ds3nfu(gca,xl(1),yl(1),zl(1))
%   set(t1,'position',[posnx,posny]);
  
  [~,minY] = min([t1p(2), t2p(2), t3p(2), t4p(2)]);
  [~,minX] = min([t1p(1), t2p(1), t3p(1), t4p(1)]);
  leftAxis = 'X'; % default
  xPts = 'WZ'; if minY <=2, xPts = 'XY'; end
  yPts = 'WY'; if minY == 1 || minY == 3, yPts = 'XZ'; end

%   if minX == 2, leftAxis = 'Y'; end % true, but not helpful
  if minX == 1, leftAxis = 'Y'; end % false, but helpful for rotation
  if minX == 2, leftAxis = 'X'; end % false, but helpful for rotation
  if minX == 3, leftAxis = 'Y'; end
  if minX == 4, leftAxis = 'X'; end
%   leftAxis
  
  % fix xlabel position
  xLabel = get(gca,'xlabel');
  xunits = get(xLabel,'units'); set(xLabel,'units','inches');
  % The following works for test cases to date, but not exhaustively
  % tested
  if leftAxis == 'X',
    theta = atan2(t3p(2)-t4p(2),t3p(1)-t4p(1));
  else
    theta = atan2(t2p(2)-t1p(2),t2p(1)-t1p(1));
  end
  set(xLabel,'VerticalAlignment','middle');
  set(xLabel,'HorizontalAlignment','center');  
  set(xLabel,'Rotation',theta*180/pi);
%   a = get(xLabel,'position')
  if strcmp(xPts,'WZ'),
    % find center position along axis
    xpos = (t3p(1)+t4p(1))/2; ypos = (t3p(2)+t4p(2))/2;
    L2 = 2*XLABELFONTSIZE/72;
    d1 = L2*cos(pi/2-abs(theta)); d2 = L2*sin(pi/2-abs(theta));    
    set(xLabel,'position',[xpos-d1 ypos-d2]);    
  else
    xpos = (t2p(1)+t1p(1))/2; ypos = (t2p(2)+t1p(2))/2;
    L2 = 2*XLABELFONTSIZE/72;
    d1 = L2*cos(pi/2-abs(theta)); d2 = L2*sin(pi/2-abs(theta));    
    set(xLabel,'position',[xpos-d1 ypos-d2]);    
  end
%   set(xLabel,'position',[a(1) -1.25*(AXISFONTSIZE+XLABELFONTSIZE)/72-1*addedSpace(2)]);
  set(xLabel,'units',xunits);

  % fix ylabel position
  yLabel = get(gca,'ylabel');
  yunits = get(yLabel,'units'); set(yLabel,'units','inches');
  % The following works for test cases to date, but not exhaustively
  % tested
  if leftAxis == 'Y',
    theta = atan2(t1p(2)-t3p(2),t1p(1)-t3p(1));
    disp('case 1')
  else
    theta = atan2(t1p(2)-t3p(2),t1p(1)-t3p(1)); % yup
  end
  set(yLabel,'VerticalAlignment','middle');
  set(yLabel,'HorizontalAlignment','center');  
  set(yLabel,'Rotation',theta*180/pi);
  if strcmp(yPts,'WY'),
    xpos = (t4p(1)+t2p(1))/2; ypos = (t4p(2)+t2p(2))/2;
    L2 = 2*YLABELFONTSIZE/72;
    d1 = L2*cos(pi/2+theta); d2 = L2*sin(pi/2+theta);    
    set(yLabel,'position',[xpos-d1 ypos-d2]);
    warning('ylabel case not yet verified...');
  else % Next case works for cases tried
    % find center position along axis
    xpos = (t3p(1)+t1p(1))/2; ypos = (t3p(2)+t1p(2))/2;
    L2 = 2*YLABELFONTSIZE/72;
    d1 = L2*cos(pi/2+theta); d2 = L2*sin(pi/2+theta);    
    set(yLabel,'position',[xpos-d1 ypos-d2]);
  end  
  set(yLabel,'units',yunits);

  % fix zlabel position
  zLabel = get(gca,'zlabel');
  zunits = get(zLabel,'units'); set(zLabel,'units','inches');
  a = get(zLabel,'position');
  set(zLabel,'VerticalAlignment','top');
  set(zLabel,'position',[-padLeft a(2)]);
  set(zLabel,'units',zunits);

  % fix title position
  tLabel = get(gca,'title');
  tunits = get(tLabel,'units'); set(tLabel,'units','inches');
  a = get(tLabel,'position');
  set(tLabel,'VerticalAlignment','bottom');
  set(tLabel,'position',[a(1) AHEIGHT + 0.2*TITLEFONTSIZE/72]);
  set(tLabel,'units',tunits);
end

set(gcf,'units',funits);
set(gcf,'paperunits',punits);
for k = 1:length(ax),
  set(gcf,'CurrentAxes',ax(k));
  set(gca,'units',aunits{k})
end
% the following is supposed to help with R2014b graphics, but doesn't seem
% to do much
D = get(gcf,'PaperPosition'); % Returns 1x4 vector [left right width height]
set(gcf,'PaperSize',[D(3) D(4)]); %default PaperSize is [8.5 11]

% Ryan Jobman code for positioning labels
hx = get(gca,'xlabel');
x_lim=get(gca,'xlim');
x_lim_scale=x_lim(2)-x_lim(1);
x_lim_mid=x_lim_scale/2+x_lim(1);
hy = get(gca,'ylabel');
y_lim=get(gca,'ylim');
y_lim_scale=y_lim(2)-y_lim(1);
y_lim_mid=y_lim_scale/2+y_lim(1);
hz = get(gca,'zlabel');
z_lim=get(gca,'zlim');
z_lim_scale=z_lim(2)-z_lim(1);
z_lim_mid=z_lim_scale/2+z_lim(1);
set(hx,'position',[x_lim_mid-.1*x_lim_scale,y_lim(1)-.4*y_lim_scale,z_lim(1)]);
set(hy,'position',[x_lim(1)-.4*x_lim_scale,y_lim_mid-.1*y_lim_scale,z_lim(1)]); 
set(hz,'position',[x_lim(1),y_lim(2)+.4*y_lim_scale,z_lim_mid-.1*z_lim_scale]);
ht = get(gca,'title');
set(ht,'position',[x_lim_mid+.35*x_lim_scale,y_lim(2),z_lim(2)+.05*z_lim_scale]);
drawnow; pause(0.01);