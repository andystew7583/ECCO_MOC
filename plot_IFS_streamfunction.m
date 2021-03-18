%%%
%%% plot_IFS_streamfunction.m
%%%
%%% Creates some plots of the IFS and overturning streamfunction derived
%%% from the ECCOv4r4 daily output.
%%%

%%% Required paths
addpath colormaps;

%%% Static definitions
isopDefinitions;

%%% Load precomputed streamfunction and IFS time series
load(fullfile(products_dir,'PSItot.mat'));
load(fullfile(products_dir,'IFS.mat'));
load(fullfile(products_dir,'WS.mat'));

%%% Coriolis parameter matrix on IFS points
Omega = 366/365*2*pi/86400;
ff = 2*Omega*sind(secLats);
FF = repmat(ff',[1 length(dens_bnds)]);

%%% Streamfunction in density space
figure(1);
[DD,LL] = meshgrid(dens_bnds,lat);
pcolor(LL,DD,mean(PSI,3)/1e6);
shading interp;
caxis([-20 20]);
set(gca,'YDir','reverse');
colorbar;
colormap redblue;

%%% IFS-reconstructed streamfunction in density space
figure(2);
[DD,LL] = meshgrid(dens_bnds,secLats);
IFSmean = mean(IFS,3);
IFSmean = IFSmean - repmat(IFSmean(:,end),[1 length(dens_bnds)]);
pcolor(LL,DD,-IFSmean./rhoConst./FF/1e6);
shading interp;
caxis([-20 20]);
set(gca,'YDir','reverse');
colorbar;
colormap redblue;



%%% Latitudinally-averaged time series metrics of overturning and IFS
ymin = -60;
ymax = -50;
dens_psimax = 1037.1;
% dens_psimax = 1037;
didx_psimax = find(dens_bnds>dens_psimax,1,'first');
yidx_psi = find((lat>ymin) & (lat<ymax));
yidx_ifs = find((secLats>ymin) & (secLats<ymax));
PSImean = squeeze(mean(PSI(yidx_psi,didx_psimax,:),1));
IFSmean = squeeze(mean(IFS(yidx_ifs,didx_psimax,:),1));
TFSmean = squeeze(mean(TFS(yidx_ifs,:),1))';
WSmean = squeeze(mean(WS(yidx_ifs,:)))';
tt = startdate + (0:Nt-1);
smoothlen = 365;

%%% Time series of streamfunction and form stresses
figure(3);
plot(tt,PSImean);
hold on;
plot(tt,-IFSmean/rhoConst/mean(ff(yidx_ifs)));
plot(tt,-TFSmean/rhoConst/mean(ff(yidx_ifs)));
hold off;
datetick('x')

%%% Time series with running monthly filters
figure(4);
plot(tt,smooth(PSImean,smoothlen));
hold on;
plot(tt,-smooth(IFSmean/rhoConst/mean(ff(yidx_ifs)),smoothlen));
plot(tt,-smooth(TFSmean/rhoConst/mean(ff(yidx_ifs)),smoothlen));
hold off;
datetick('x')

%%% Actual vs predicted AABW transports with running monthly filters
figure(5);
plot(tt,smooth(PSImean,smoothlen));
hold on;
plot(tt,-smooth((IFSmean-TFSmean)/rhoConst/mean(ff(yidx_ifs)),smoothlen));
hold off;
datetick('x')


%%% Actual vs predicted AABW transports with running monthly filters
figure(6);
plot(tt,-smooth(PSImean-mean(PSImean),smoothlen));
hold on;
plot(tt,smooth((WSmean-mean(WSmean))/rhoConst/mean(ff(yidx_ifs)),smoothlen));
hold off;
datetick('x')

figure(7);
plot(tt,-smooth(PSImean-mean(PSImean),smoothlen));
hold on;
plot(tt,-smooth((TFSmean-mean(TFSmean))/rhoConst/mean(ff(yidx_ifs)),smoothlen));
hold off;
datetick('x')

figure(8);
scatter(smooth(PSImean,smoothlen),-smooth((WSmean)/rhoConst/mean(ff(yidx_ifs)),smoothlen));

%%% Time series correlations
corr(smooth(IFSmean,smoothlen),smooth(PSImean,smoothlen))
corr(-smooth(TFSmean,smoothlen),smooth(PSImean,smoothlen))
corr(smooth(IFSmean-TFSmean,smoothlen),smooth(PSImean,smoothlen))
corr(smooth(WSmean,smoothlen),smooth(PSImean,smoothlen))
corr(smooth(WSmean,smoothlen),smooth(TFSmean,smoothlen))


lag = 1;
corr(WSmean(1:end-lag),PSImean(1+lag:end))

%%% Barotropic momentum balance
figure(9);
plot(secLats,mean(WS,2));
hold on;
plot(secLats,mean(TFS,2));
plot(secLats,mean(IFS(:,didx_psimax,:),3));
hold off;
legend('Wind Stress','Topographic Form Stress');

%%% Cross-correlation between wind and TFS
[r,lags]=xcorr(WSmean,TFSmean,30);
figure(10);
plot(lags,r);

figure(11)
plot(dens_bnds-1037,std(squeeze(mean(IFS(yidx_ifs,:,:),1)),[],2));
xlabel('\sigma_2 - 1037 kg/m^3');
ylabel('std(IFS)')

figure(12)
plot(dens_bnds-1037,mean(squeeze(mean(IFS(yidx_ifs,:,:),1)),2));
xlabel('\sigma_2 - 1037 kg/m^3');
ylabel('mean(IFS)')




