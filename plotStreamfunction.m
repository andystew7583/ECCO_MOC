%%% 
%%% plotStreamfunction.m
%%%
%%% Creates various plots of various streamfunctions for testing purposes.
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Start by clearing memory
clear all;

%%% Required paths
addpath colormaps;
addpath CDT/cdt;

%%% Load global variables
isopDefinitions;

%%% Add required paths
p = genpath('gcmfaces/'); addpath(p);
p = genpath('gcmfaces/'); addpath(p);
p = genpath('m_map/'); addpath(p);
p = genpath('rw_namelist/'); addpath(p);

%%% Load all grid variables from nctiles_grid/ into mygrid
grid_load([ECCO_grid_dir filesep],5,'nctiles',0,1);

%%% Make mygrid accessible in current workspace:
gcmfaces_global;

%%% Directory from which to read ECCO output data
myenv.nctilesdir = fullfile(ECCO_data_dir,filesep);







%%% Plotting options
fontsize = 14;
scrsz = get(0,'ScreenSize');
framepos = [0.25*scrsz(3) 0.25*scrsz(4) 1000 500];
fignum = 0;

%%% Isopycnal depths
load([products_dir filesep 'Zisop.mat']);
Nlats = length(lat);

%%% Total (i.e. residual) streamfunction
load([products_dir filesep 'PSItot.mat']);
% load([products_dir filesep 'PSIeul.mat']);
% PSIbol = PSItot - PSIeul;
% PSIbol_mean = PSItot_mean - PSIeul_mean;
% PSIbol_zonmean = PSItot_zonmean - PSIeul_zonmean;

%%% Select streamfunction to plot
psitype = 'tot';
eval(['PSI = PSI',psitype,';']);
eval(['PSI_mean = PSI',psitype,'_mean;']);
eval(['PSI_zonmean = PSI',psitype,'_zonmean;']);

%%% Adjusted streamfunction (enforces psi=0 at sea floor)
PSI_adj = PSI-repmat(PSI(:,end,:),[1 Nd+1 1]);

%%% Time-averaged streamfunction
PSI_tavg = mean(PSI,3);

%%% De-seasonalized streamfunction
PSI_anom = PSI - PSI_tavg;
PSI_seasononly = repmat(squeeze(mean(reshape(PSI_anom,[Nlats Nd+1 12 Nt/12]),4)),[1 1 Nt/12]);
PSI_noseason = PSI - PSI_seasononly;

%%% De-seasoned and de-trended streamfunction
PSI_trend = zeros(Nlats,Nd+1);
PSI_trend_rval = zeros(Nlats,Nd+1);
PSI_trend_pval = zeros(Nlats,Nd+1);
PSI_trendonly = zeros(Nlats,Nd+1,Nt);
PSI_anom = PSI_noseason - PSI_tavg;
for j=1:Nlats
  j
  for k=1:Nd+1
    coeff = polyfit(tt',squeeze(PSI_anom(j,k,:)),1);
    PSI_trend(j,k) = coeff(1);
    [r,p] = corr(coeff(2)+coeff(1)*tt',squeeze(PSI_anom(j,k,:)));
    PSI_trend_rval(j,k) = r;
    PSI_trend_pval(j,k) = p;
    PSI_trendonly(j,k,:) = reshape(coeff(2) + coeff(1)*tt,[1 1 Nt]);
  end
end
PSI_notrend = PSI_noseason - PSI_trendonly;

%%% Calculate EOFs of de-seasoned, de-trended, de-meaned streamfunction
PSI_anom = PSI_notrend - PSI_tavg;
% PSI_anom = PSI_noseason - PSI_tavg;
eof_msk = (PSI_tavg~=0);
eof_idx = find(lat>-40);
eof_msk(eof_idx,:) = 0;
[eof_maps,pc,expvar] = eof(PSI_anom,'mask',eof_msk);

%%% Interpolate zonal-mean streamfunction onto mean isopycnal depths
PSI_zonmean_isop = zeros(Nlats,Nd+1);
for j = 1:Nlats
  PSI_zonmean_isop(j,:) = interp1(mygrid.RF,PSI_zonmean(j,:),Zisop_mean(j,:),'linear','extrap');
end

%%% Residual streamfunction in density space
[DD,LL] = meshgrid(dens_bnds-1000,lat);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,DD,mean(PSI,3)/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-40 40]);
axis([-80 80 29 38]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('\sigma_2 (kg/m^3)');
title('Residual overturning streamfunction in density coordinates');

%%% Adjusted residual streamfunction in density space
[DD,LL] = meshgrid(dens_bnds-1000,lat);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,DD,mean(PSI_adj,3)/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-40 40]);
axis([-80 80 29 38]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('\sigma_2 (kg/m^3)');
title('Adjusted residual overturning streamfunction in density coordinates');

%%% Residual streamfunction in density space, time mean component
[DD,LL] = meshgrid(dens_bnds-1000,lat);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,DD,PSI_mean/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-40 40]);
axis([-80 80 29 38]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('\sigma_2 (kg/m^3)');
title('Residual overturning streamfunction in density coordinates, time-mean component');

%%% Residual streamfunction in density space, time mean component
[DD,LL] = meshgrid(dens_bnds-1000,lat);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,DD,(mean(PSI,3)-PSI_mean)/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-10 10]);
axis([-80 80 29 38]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('\sigma_2 (kg/m^3)');
title('Residual overturning streamfunction in density coordinates, fluctuating component');

%%% Residual overturning in depth space
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,mean(PSI_noseason,3)/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-20 20]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Residual overturning streamfunction in depth coordinates');

%%% Eulerian-mean overturning
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,PSI_zonmean_isop/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-20 20]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Eulerian-mean overturning streamfunction');

%%% Standing eddy overturning
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,(mean(PSI,3)-PSI_zonmean_isop)/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-20 20]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Standing eddy overturning streamfunction');

%%% Standard deviation of residual overturning
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,std(PSI,0,3)/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-20 20]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Standard deviation of residual overturning streamfunction');

%%% Standard deviation of de-seasonalized residual overturning
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,std(PSI_noseason,0,3)/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-20 20]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Standard deviation of de-seasonalized residual overturning streamfunction');

%%% Standard deviation of de-deasoned and de-trended residual overturning
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,std(PSI_notrend,0,3)/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-20 20]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Standard deviation of de-seasoned and de-trended residual overturning streamfunction');

%%% Linear trend of residual overturning
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,PSI_trend*Nt/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-10 10]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Linear trend in the residual overturning streamfunction');

%%% r^2 for linear trend 
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,PSI_trend_rval.^2);
shading interp;
colormap(hot(100));
colorbar;
% caxis([- 10]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Linear trend r^2');

%%% p-value for linear trend 
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,PSI_trend_pval);
shading interp;
colormap(hot(100));
colorbar;
% caxis([- 10]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Linear trend p-value');

%%% EOF1
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,eof_maps(:,:,end));
shading interp;
colormap(redblue(100));
colorbar;
caxis([-.05 .05]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title(['EOF1, expvar=',num2str(expvar(end))]);

%%% EOF2
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,eof_maps(:,:,end-1));
shading interp;
colormap(redblue(100));
colorbar;
caxis([-.05 .05]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title(['EOF2, expvar=',num2str(expvar(end-1))]);

%%% EOF1 in density space
[DD,LL] = meshgrid(dens_bnds-1000,lat);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,DD,eof_maps(:,:,end));
shading interp;
colormap(redblue(100));
colorbar;
caxis([-.05 .05]);
axis([-80 -40 34 37.5]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title(['EOF1, expvar=',num2str(expvar(end))]);

%%% EOF2 in density space
[DD,LL] = meshgrid(dens_bnds-1000,lat);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,DD,eof_maps(:,:,end-1));
shading interp;
colormap(redblue(100));
colorbar;
caxis([-.05 .05]);
axis([-80 -40 34 37.5]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title(['EOF2, expvar=',num2str(expvar(end-1))]);



%%% Final residual overturning in depth space
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,(PSI_tavg+PSI_trendonly(:,:,end))/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-20 20]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Final residual overturning streamfunction in depth coordinates');



%%% Initial residual overturning in depth space
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,(PSI_tavg+PSI_trendonly(:,:,11))/1e6);
shading interp;
colormap(cmocean('balance',40));
colorbar;
caxis([-20 20]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
xlabel('Latitude');
ylabel('Depth (m)');
title('Initial residual overturning streamfunction in depth coordinates');





%%% Pretty version of the adjusted residual streamfunction in density space
LL = repmat([lat]',[1 Nd+1]);
fignum = fignum + 1;
handle = figure(fignum);
set(handle,'Position',framepos);
pcolor(LL,-Zisop_mean,mean(PSI_adj,3)/1e6);
shading interp;
colormap(cmocean('balance',40));
% colormap(redblue(40));
colorbar;
caxis([-20 20]);
axis([-80 80 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
set(gca,'Color',[.8 .8 .8]);
set(gca,'Position',[0.07 0.1 0.84 0.85]);
xlabel('Latitude');
ylabel('Depth (m)');
% title('Residual overturning streamfunction in depth coordinates');

Zisop_idx = [53 83 108 116 125];
Zisop_plot = Zisop_mean;
for j=1:Nlats
  idx = find(Zisop_plot(j,:) == 0);
  Zisop_plot(j,idx) = NaN;
  idx = find(Zisop_plot(j,:) == Zisop_plot(j,end));
  Zisop_plot(j,idx) = NaN;
end
Zisop_plot(lat>0,Zisop_idx(end)) = NaN; %%% Get rid of distracting deep density in NH

hold on;
for k=Zisop_idx
  plot([lat],-Zisop_plot(:,k),'Color',[.5 .5 .5],'LineWidth',1);  
%   plot([lat],-Zisop_plot(:,k),'Color',[1 1 1],'LineWidth',1);  
end
hold off;

figure1 = gcf;
annotation(figure1,'arrow',[0.802 0.7],[0.676 0.648]);
annotation(figure1,'arrow',[0.262 0.203],[0.684 0.804000000000004]);
annotation(figure1,'arrow',[0.261000000000001 0.32],...
  [0.878000000000001 0.85]);
annotation(figure1,'arrow',[0.687 0.777],[0.878 0.916]);
annotation(figure1,'arrow',[0.17 0.217],[0.61 0.452]);
annotation(figure1,'arrow',[0.221 0.178],[0.604 0.776]);
annotation(figure1,'arrow',[0.618000000000001 0.504],...
  [0.482000000000001 0.492]);
annotation(figure1,'arrow',[0.491000000000001 0.599],...
  [0.248000000000001 0.264]);