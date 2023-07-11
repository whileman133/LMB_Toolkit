% demoThesisFormat.m
%
% Demonstrate usage of the THESISFORMAT function for generating
% publication-quality figures in MATLAB. Hip hip hoorah!
%
% -- Changelog --
% 2023.06.11 | Make compatible with legacy thesisFormat | Wesley Hileman
% 2023.06.08 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;

% The following determines whether the generated figures will be saved 
% to disk. Note: 
% - thesisFormat defaults to saving figures in both eps and png formats.
% - if any directory does not yet exist, thesisfig creates it automatically.
saveflag = false;
savedir = 'demothesisfig';

% Generate data to plot...
freq = logspace(-3,3,50);
s = 1j*2*pi*freq;
Z1 = 2e-3*s + 10./s./0.01./(10 + 1./s./0.01);
Z2 = 1./10./s + 10./s./0.01./(10 + 1./s./0.01);

% Plotting ----------------------------------------------------------------

% Basic nyquist plot.
figure;
plot(real(Z1),imag(Z1),'ro-');
xlabel('Z_1'' [\Omega]');
ylabel('Z_1'''' [\Omega]');
title('Nyquist: Z_1(f)');
thesisFormat('AxesLimits','Nyquist', ...
    'SaveFlag',saveflag,'SaveName',fullfile(savedir,'1-nyquist'));

% If you need to fudge the label margins (hopefully a rare event), use 
% the 'PlotBoxPaddingInches' option, which specifies additional padding 
% around each plot axes in the form [left bottom right top].
figure;
plot(real(Z1),imag(Z1),'ro-');
xlabel('Z_1'' [\Omega]');
ylabel('Z_1'''' [\Omega]');
title('Nyquist: Z_1(f) (Fudged Plot Margins)');
thesisFormat('AxesLimits','Nyquist','PlotBoxPaddingInches',[.3 .3 .3 .3], ...
    'SaveFlag',saveflag,'SaveName',fullfile(savedir,'2-fudgedmargins'));

% And yes, it can even handle line breaks within titles/labels!
figure;
plot(real(Z1),imag(Z1),'ro-');
xlabel(['Z_1'' [\Omega]' sprintf('\n(Something else here)')]);
ylabel(['Z_1'''' [\Omega]' sprintf('\n(Something more here)')]);
title(sprintf('Nyquist: Z_1(f)\n(A special note here)'));
thesisFormat('AxesLimits','Nyquist', ...
    'SaveFlag',saveflag,'SaveName',fullfile(savedir,'3-labellinebreaks'));

% Two sublots
figure;
subplot(1,2,1);
plot(real(Z1),imag(Z1),'ro-');
xlabel('Z_1'' [\Omega]');
ylabel('Z_1'''' [\Omega]');
title('Nyquist: Z_1(f)');
subplot(1,2,2);
plot(real(Z2),imag(Z2),'bd-');
xlabel('Z_2'' [\Omega]');
ylabel('Z_2'''' [\Omega]');
title('Nyquist: Z_2(f)');
thesisFormat('AxesLimits','Nyquist', ...
    'SaveFlag',saveflag,'SaveName',fullfile(savedir,'4-twosubplots'));

% Six subplots.
% Axes limits are not be set using the 'AxesLimits','Nyquist' option 
% since four of the subplots are bode plots; use setAxesNyquist on the
% Nqyuist subplots instead.
figure;
subplot(3,2,1);
plot(real(Z1),imag(Z1),'ro-');
xlabel('Z_1'' [\Omega]');
ylabel('Z_1'''' [\Omega]');
title('Nyquist: Z_1(f)');
setAxesNyquist;
subplot(3,2,3);
loglog(freq,abs(Z1),'ro-');
xlabel('Cyclic Frequency, f [Hz]');
ylabel('||Z_1(f)|| [\Omega]');
title('Bode Magnitude: Z_1(f)');
subplot(3,2,5);
semilogx(freq,angle(Z1)*180/pi,'ro-');
xlabel('Cyclic Frequency, f [Hz]');
ylabel('\angleZ_1(f) [\circ]');
title('Bode Phase: Z_1(f)');
subplot(3,2,2);
plot(real(Z2),imag(Z2),'bd-');
xlabel('Z_2'' [\Omega]');
ylabel('Z_2'''' [\Omega]');
title('Nyquist: Z_2(f)');
setAxesNyquist;
subplot(3,2,4);
loglog(freq,abs(Z2),'bd-');
xlabel('Cyclic Frequency, f [Hz]');
ylabel('||Z_2(f)|| [\Omega]');
title('Bode Magnitude: Z_2(f)');
subplot(3,2,6);
semilogx(freq,angle(Z2)*180/pi,'bd-');
xlabel('Cyclic Frequency, f [Hz]');
ylabel('\angleZ_2(f) [\circ]');
title('Bode Phase: Z_2(f)');
thesisFormat( ...
    'SaveFlag',saveflag,'SaveName',fullfile(savedir,'5-sixsubplots'));

% Tiled layout: fill entire US letter with plots. Useful when you need to
% show a lot of figures on the same page.
figure;
l = tiledlayout(6,3);
l.Title.String = 'Shared Title';
l.XLabel.String = 'Shared X-Label';
l.YLabel.String = 'Shared Y-Label';
for k = 1:floor(prod(l.GridSize)/2)
    nexttile; plot(1:10,1:10,'-bo');
    nexttile; plot(1:10,10:-1:1,'-md');
end
thesisFormat( ...
    'FigSizeInches',[8.5 11], ...  % US letter size
    'FigMarginInches',[0.25 0.25 0.25 0.25], ... % the default anyway
    'SaveFlag',saveflag,'SaveName',fullfile(savedir,'6-tiledlayout'));