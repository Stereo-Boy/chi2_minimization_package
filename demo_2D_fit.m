try

close all, clear all;
thisPath=fileparts(mfilename('fullpath')); %path to folder where code is located
shared_path = fullfile(thisPath,'shared_functions'); %path to common functions
addpath(genpath(shared_path));  % add that path to use these functions
fun_path = fullfile(thisPath,'fitting_functions'); %path to fitting functions
check_folder(fun_path,1,'verboseON'); addpath(genpath(fun_path)); % add that path to use these functions

x = -10:0.1:10; % data range for a square grid
[xx,yy] = meshgrid(x,x); % data square grid
params = [300,0,0,1,3,deg2rad(20),20]; % parameters of a simulated 2D Gaussian
zz = reshape(elliptical_2D_gaussian_withd([xx(:),yy(:)],params),size(xx)); % actual simulated 2D Gaussian for that square grid

% plot the ideal
subplot(1,3,1)
contourf(xx,yy,zz); title('Simulated data'); colorbar

% create and plot the data
fraction2keep = 0.5; % we keep only 50% of the data
noiseLevel = 0.1;   % noise in %
sampling = rand(size(xx(:)))<fraction2keep;
xxData = xx(sampling); yyData = yy(sampling); zzData = zz(sampling)+rand(size(zz(sampling))).*noiseLevel-noiseLevel/2;
subplot(1,3,2)
zz2plot = zz; zz2plot(~sampling) = nan; zz2plot(sampling) = zzData; zz2plot = zz2plot(end:-1:1,:); % invert x
imagesc(zz2plot); title('Noisy data'); colorbar

% normalization is needed
mn = 1; mx = 1300; amp = mx - mn;
zzData = (zzData - mn)./amp; % now this varies between 0 and 1

% fit the data
opts = optimset('display','off');
pMin = [0.005,min(x),min(x),0.1,0.1,deg2rad(0),0];
pMax = [0.70,max(x),max(x),10,10,deg2rad(180),0.02];
[bestParams,bestChiSq,zzModel,SE,stdFits] = chi2minimFit(@elliptical_2D_gaussian_withd, pMin, pMax, [xxData(:),yyData(:)], zzData(:), 20*ones(size(zzData)), 50, 'verboseON', [xx(:),yy(:)],opts);

% retransform the data
zzModel = zzModel.*amp+mn;
bestParams([1,7]) = bestParams([1,7]).*amp+mn;
dispi('Parameters simulated: ', params)
dispi('Retransformed parameters: ', bestParams)

% plot the fit
subplot(1,3,3)
contourf(xx,yy,reshape(zzModel,size(xx))); title('Fitted model'); colorbar

catch err
    keyboard
end