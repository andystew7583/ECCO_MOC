%%%
%%% plot_IFS_streamfunction.m
%%%
%%% Creates some plots of the IFS and overturning streamfunction derived
%%% from the ECCOv4r4 daily output.
%%%

%%% Required paths
addpath colormaps;
addpath CDT/cdt;

%%% Static definitions
isopDefinitions;

%%% Load precomputed streamfunction and IFS time series
load(fullfile(products_dir,'PSIbol.mat'))
PSIbol = PSI;
load(fullfile(products_dir,'PSIeul.mat'))
PSIeul = PSI;
load(fullfile(products_dir,'PSItot.mat'))
load(fullfile(products_dir,'IFS.mat'));
load(fullfile(products_dir,'WS.mat'));

%%% Clean up streamfunctions
PSIbol_mean = mean(PSIbol,3);
for j = 1:length(lat)
  for m=Nd:-1:1
    if (PSIbol_mean(j,m)==PSIbol_mean(j,m+1))
      continue;
    end
    PSIbol_mean(j,m+2:end) = 0;
    break;
  end
end

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

%%% Streamfunction in density space
figure(20);
[DD,LL] = meshgrid(dens_bnds,lat);
pcolor(LL,DD,mean(PSIeul,3)/1e6);
shading interp;
caxis([-20 20]);
set(gca,'YDir','reverse');
colorbar;
colormap redblue;

%%% Streamfunction in density space
figure(21);
[DD,LL] = meshgrid(dens_bnds,lat);
pcolor(LL,DD,mean(PSIbol_mean,3)/1e6);
shading interp;
caxis([-20 20]);
set(gca,'YDir','reverse');
colorbar;
colormap redblue;

EFSmean = zeros(length(secLats),Nd+1);
for m=1:Nd+1
  EFSmean(:,m) = interp1(lat,mean(PSIbol(:,m),3),secLats,'linear');
end
EFSmean = -EFSmean.*rhoConst.*FF;

%%% IFS-reconstructed streamfunction in density space
figure(2);
[DD,LL] = meshgrid(dens_bnds,secLats);
IFSmean = mean(IFS,3);
IFSmean = IFSmean - repmat(IFSmean(:,end),[1 length(dens_bnds)]);
pcolor(LL,DD,-(IFSmean+EFSmean)./rhoConst./FF/1e6);
% pcolor(LL,DD,-(EFSmean)./rhoConst./FF/1e6);
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
PSIeul_mean = squeeze(mean(PSIeul(yidx_psi,didx_psimax,:),1));
PSIbol_mean = squeeze(mean(PSIbol(yidx_psi,didx_psimax,:),1));
IFSmean = squeeze(mean(IFS(yidx_ifs,didx_psimax,:),1));
EFSmean = -PSIbol_mean*rhoConst*mean(ff(yidx_ifs));
TFSmean = squeeze(mean(TFS(yidx_ifs,:),1))';
WSmean = squeeze(mean(WS(yidx_ifs,:)))';
tt = startdate + (0:Nt-1);
smoothlen = 30;

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
plot(tt,-smooth((EFSmean+IFSmean-TFSmean)/rhoConst/mean(ff(yidx_ifs)),smoothlen));
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
corr(smooth(IFSmean+EFSmean-TFSmean,smoothlen),smooth(PSImean,smoothlen))
corr(smooth(WSmean,smoothlen),smooth(PSImean,smoothlen))
corr(smooth(WSmean,smoothlen),smooth(TFSmean,smoothlen))


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


eof_msk = (mean(PSI,3)~=0);
eof_idx = find(lat>-40);
eof_msk(eof_idx,:) = 0;
[eof_maps,pc,expvar] = eof(PSI,10,'mask',eof_msk);

figure(13);
[DD,LL] = meshgrid(dens_bnds,lat);
pcolor(LL,DD,eof_maps(:,:,1));
shading interp;
colormap redblue;
set(gca,'YDir','reverse');
set(gca,'XLim',[-70 -40]);
set(gca,'YLim',[1032 1038]);
caxis([-0.04 0.04])

figure(14);
[DD,LL] = meshgrid(dens_bnds,lat);
pcolor(LL,DD,eof_maps(:,:,2));
shading interp;
colormap redblue;
set(gca,'YDir','reverse');
set(gca,'XLim',[-70 -40]);
set(gca,'YLim',[1032 1038]);
caxis([-0.04 0.04])

figure(15);
[DD,LL] = meshgrid(dens_bnds,lat);
pcolor(LL,DD,eof_maps(:,:,3));
shading interp;
colormap redblue;
set(gca,'YDir','reverse');
set(gca,'XLim',[-70 -40]);
set(gca,'YLim',[1032 1038]);
caxis([-0.04 0.04])

figure(16);
plot(pc(1,:))

corr(smooth(pc(1,:)',smoothlen),smooth(IFSmean-TFSmean,smoothlen))
corr(smooth(pc(2,:)',smoothlen),smooth(IFSmean-TFSmean,smoothlen))



%%% RF calculation
Nrf = 10*365;
PSI_WS_RF = calcRF(WSmean,PSImean*rhoConst*mean(ff(yidx_ifs)),Nrf);
[TFS_WS_RF,TFS_WS_eps] = calcRF(WSmean,TFSmean,Nrf);
IFS_WS_RF = calcRF(WSmean,IFSmean,Nrf);
EFS_WS_RF = calcRF(WSmean,EFSmean,Nrf);
RFS_WS_RF = calcRF(WSmean,IFSmean+EFSmean-TFSmean,Nrf);
tt_rf = 0:Nrf-1;

figure(17);
plot(tt_rf,cumsum(PSI_WS_RF));

figure(18);
plot(tt_rf,cumsum(TFS_WS_RF));

figure(19);
plot(tt_rf,cumsum(IFS_WS_RF));

figure(20);
plot(tt_rf,cumsum(EFS_WS_RF));


figure(21);
plot(tt_rf,cumsum(RFS_WS_RF));


Lwin = 21*365;
Dwin = 30;
Nrf = 10*365;
Nwin = floor((Nt-Lwin)/Dwin) + 1;
TFS_WS_RF = zeros(Nrf,Nwin);
for n = 1:Nwin
  n
  idx = (n-1)*Dwin + (1:Lwin);
  [rf,eps] = calcRF(WSmean(idx),TFSmean(idx),Nrf);
  TFS_WS_RF(:,n) = rf;
end

figure(22)
plot(tt_rf,mean(cumsum(TFS_WS_RF,1),2));
hold on
plot(tt_rf,mean(cumsum(TFS_WS_RF,1),2)+std(cumsum(TFS_WS_RF,1),[],2),'k--');
plot(tt_rf,mean(cumsum(TFS_WS_RF,1),2)-std(cumsum(TFS_WS_RF,1),[],2),'k--');
hold off

figure(23)
plot(tt_rf,smooth(cumsum(TFS_WS_RF(:,1),1),30))
hold on
for n = 2:Nwin
  plot(tt_rf,smooth(cumsum(TFS_WS_RF(:,n),1),30))
end
hold off