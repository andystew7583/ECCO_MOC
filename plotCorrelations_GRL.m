%%%
%%% plotCorrelations_GRL.m
%%%
%%% Plots relationships between wind, form stresses and overturning for our
%%% GRL paper.
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

%%% Coriolis parameter matrix on IFS points
Omega = 366/365*2*pi/86400;
ff = 2*Omega*sind(secLats);
FF = repmat(ff',[1 length(dens_bnds)]);

%%% Eddy component of IFS
EFSmean = zeros(length(secLats),Nd+1);
for m=1:Nd+1
  EFSmean(:,m) = interp1(lat,mean(PSIbol(:,m),3),secLats,'linear');
end
EFSmean = -EFSmean.*rhoConst.*FF;

%%% Latitudinally-averaged time series metrics of overturning and IFS
ymin = -60;
ymax = -50;
dens_psimax = 1037.1;
didx_psimax = find(dens_bnds>dens_psimax,1,'first');
yidx_psi = find((lat>ymin) & (lat<ymax));
yidx_ifs = find((secLats>ymin) & (secLats<ymax));
PSImean = squeeze(mean(PSI(yidx_psi,didx_psimax,:),1));
PSIeul_mean = squeeze(mean(PSIeul(yidx_psi,didx_psimax,:),1));
PSIbol_mean = squeeze(mean(PSIbol(yidx_psi,didx_psimax,:),1));
RFSmean = squeeze(mean(IFS(yidx_ifs,didx_psimax,:),1));
EFSmean = -PSIbol_mean*rhoConst*mean(ff(yidx_ifs));
TFSmean = squeeze(mean(TFS(yidx_ifs,:),1))';
WSmean = squeeze(mean(WS(yidx_ifs,:)))';
IFSmean = RFSmean + EFSmean;
tt = startdate + (0:Nt-1);

%%% Compute correlations and NSE between wind, form stresses, and AABW
%%% export
% smoothlens = [1:2:5*365];
smoothlens = [1:2:10 15:5:100 150:50:1000 1100:100:1500];
[r_WS_PSI,p_WS_PSI,rmin_WS_PSI,NSE_WS_PSI] = calcRelationships(WSmean,-PSImean*rhoConst*mean(ff(yidx_ifs)),smoothlens);
[r_WS_TFS,p_WS_TFS,rmin_WS_TFS,NSE_WS_TFS] = calcRelationships(WSmean,-TFSmean,smoothlens);
[r_WS_IFS,p_WS_IFS,rmin_WS_IFS,NSE_WS_IFS] = calcRelationships(WSmean,-IFSmean,smoothlens);
[r_PSI_TFS,p_PSI_TFS,rmin_PSI_TFS,NSE_PSI_TFS] = calcRelationships(PSImean*rhoConst*mean(ff(yidx_ifs)),TFSmean,smoothlens);
[r_PSI_IFS,p_PSI_IFS,rmin_PSI_IFS,NSE_PSI_IFS] = calcRelationships(-PSImean*rhoConst*mean(ff(yidx_ifs)),IFSmean,smoothlens);
[r_PSI_TFSIFS,p_PSI_TFSIFS,rmin_PSI_TFSIFS,NSE_PSI_TFSIFS] = calcRelationships(-PSImean*rhoConst*mean(ff(yidx_ifs)),IFSmean-TFSmean,smoothlens);








%%% Plotting options
fontsize = 14;
framepos = [417    526   791   800];
labelspacing = 200;
axpos = zeros(4,4);
axpos(1,:) = [0.08 0.57 0.4 0.4];
axpos(2,:) = [0.58 0.57 0.4 0.4];
axpos(3,:) = [0.08 0.3 0.9 0.19];
axpos(4,:) = [0.08 0.06 0.9 0.19];
cbpos = [0.95 0.37 0.015 0.6];
axlabels = {'(a)','(b)','(c)','(c)'};
linewidth = 1.5;

%%% Set up figure window
handle = figure(203);
clf;
set(handle,'Position',framepos);



subplot('Position',axpos(1,:));
colororder = get(gca,'ColorOrder');
semilogx(smoothlens,r_WS_PSI,'Color',colororder(1,:),'LineWidth',linewidth);
hold on;
semilogx(smoothlens,r_WS_TFS,'Color',colororder(4,:),'LineWidth',linewidth);
semilogx(smoothlens,r_WS_IFS,'Color',colororder(3,:),'LineWidth',linewidth);
ahandle = area(smoothlens,rmin_WS_PSI);
ahandle.FaceColor = [.8 .8 .8];
ahandle.FaceAlpha = 0.5;
ahandle.EdgeColor = 'None';
hold off;
xlabel('Smoothing window (days)','interpreter','latex');
ylabel('Correlation with Wind Stress (WS)','interpreter','latex');
set(gca,'XTick',([1 3 10 30 100 300 1000]));
set(gca,'FontSize',fontsize);
legend('$-T_\mathrm{AABW}$','TFS','IFS$_\mathrm{AABW}$','interpreter','latex','Location','SouthWest');
grid on;

subplot('Position',axpos(2,:));
semilogx(smoothlens,r_PSI_TFS,'Color',colororder(5,:),'LineWidth',linewidth);
hold on;
semilogx(smoothlens,r_PSI_IFS,'Color',colororder(6,:),'LineWidth',linewidth);
semilogx(smoothlens,r_PSI_TFSIFS,'Color',colororder(7,:),'LineWidth',linewidth);
ahandle = area(smoothlens,rmin_PSI_TFS);
ahandle.FaceColor = [.8 .8 .8];
ahandle.FaceAlpha = 0.5;
ahandle.EdgeColor = 'None';
hold off;
xlabel('Smoothing window (days)','interpreter','latex');
ylabel('Correlation with AABW flux ($T_\mathrm{AABW}$)','interpreter','latex');
set(gca,'XTick',([1 3 10 30 100 300 1000]));
set(gca,'FontSize',fontsize);
legend('$-$TFS','IFS$_\mathrm{AABW}$','IFS$_\mathrm{AABW}-$TFS','interpreter','latex','Location','West');
grid on;

%%% Actual vs predicted AABW transports with running monthly filters
subplot('Position',axpos(3,:));
plot(tt,-smooth(PSImean-mean(PSImean),30)/1e6,'Color',colororder(1,:),'LineWidth',linewidth);
hold on;
plot(tt,smooth((WSmean-mean(WSmean))/rhoConst/mean(ff(yidx_ifs)),30)/1e6,'Color',colororder(2,:),'LineWidth',linewidth);
hold off;
yhandle = ylabel('AABW flux anomaly $T^\prime_\mathrm{AABW}$ (Sv)','interpreter','latex');
yhandle.Position = [727355.5635     -20.88813834               -1];
set(gca,'FontSize',fontsize);
datetick('x')
set(gca,'XLim',[tt(15) tt(end/2)]);
set(gca,'YLim',[-19 22]);
leghandle = legend('Diagnosed','Reconstructed from wind stress','interpreter','latex','Location','NorthEast');
set(leghandle,'Orientation','horizontal');

%%% Actual vs predicted AABW transports with running monthly filters
subplot('Position',axpos(4,:));
plot(tt,-smooth(PSImean-mean(PSImean),30)/1e6,'Color',colororder(1,:),'LineWidth',linewidth);
hold on;
plot(tt,smooth((WSmean-mean(WSmean))/rhoConst/mean(ff(yidx_ifs)),30)/1e6,'Color',colororder(2,:),'LineWidth',linewidth);
hold off;
xlabel('Time','interpreter','latex');
% ylabel('AABW flux (T_A_A_B_W)','interpreter','latex');
set(gca,'FontSize',fontsize);
datetick('x')
set(gca,'XLim',[tt(end/2+1) tt(end-14)]);
set(gca,'YLim',[-19 22]);

%%% Add panel labels
for n=[1 2 4]
  annotation('TextBox',[axpos(n,1)-0.05 axpos(n,2)-0.05 0.03 0.03],'String',axlabels{n},'EdgeColor','None','FontSize',fontsize);
end








%%% Map correlations in latitude/density space
smoothlen = 30;
Wsmooth = 5;
r_map = NaN*ones(length(secLats),length(dens_bnds));
NSE_map = NaN*ones(length(secLats),length(dens_bnds));
for j=1:length(secLats)
  if (isnan(secLats(j)) || (secLats(j)<-70) || (secLats(j)>60))
    continue;
  end  
  secLats(j)
  ymin = secLats(j)-Wsmooth/2;
  ymax = secLats(j)+Wsmooth/2;
  yidx_psi = find((lat>ymin) & (lat<ymax));
  yidx_ifs = find((secLats>ymin) & (secLats<ymax));
  for k=1:length(dens_bnds)        
%     [r_yd,p_yd,rmin_yd,NSE_yd] = calcRelationships(squeeze(mean(WS(yidx_ifs,:),1))',-squeeze(mean(PSI(yidx_psi,k,:),1)*rhoConst*mean(ff(yidx_ifs))),smoothlens);
    WS_smooth = smooth(squeeze(mean(WS(yidx_ifs,:),1))',smoothlen);
    PSI_smooth = smooth(-squeeze(mean(PSI(yidx_psi,k,:),1)*rhoConst*mean(ff(yidx_ifs))),smoothlen);
    r_yd = corr(WS_smooth,PSI_smooth);
    r_map(j,k) = r_yd;
    
    %%% Compute NSE
    WS_smooth = WS_smooth - mean(WS_smooth);
    PSI_smooth = PSI_smooth - mean(PSI_smooth);
    NSE_map(j,k) = 1 - sum((WS_smooth-PSI_smooth).^2)/sum((PSI_smooth-mean(PSI_smooth)).^2);
  end
end

%%% Blank out areas with no MOC
NSE_map(isnan(r_map)) = NaN;






%%% Plotting options
fontsize = 14;
framepos = [417    526   791   800];
labelspacing = 200;
axpos = zeros(2,4);
axpos(1,:) = [0.08 0.57 0.9 0.4];
axpos(2,:) = [0.08 0.08 0.9 0.4];
axlabels = {'(a)','(b)'};
linewidth = 1.5;

%%% Set up figure window
handle = figure(205);
clf;
set(handle,'Position',framepos);


dticks = [11 19 27 36 43 53 63 75 87 99 109 119 129];
dticklabels = {'28','30','32','34','35','35.5','36','36.3','36.6','36.9','37','37.1','37.2'};
dval = find(dens_levs==1037.1); %%% Index at which the density is equal to our AABW density threshold


subplot('Position',axpos(1,:));
[DD,LL] = meshgrid(1:Nd+1,secLats);
pcolor(LL,DD,r_map);
shading interp;
colormap(gca,cmocean('balance',20));
colorbar;
% set(gca,'XLim',[ymin ymax_SO]);
set(gca,'YLim',[dticks(1) dticks(end)]);
set(gca,'YTick',dticks);
set(gca,'YTickLabel',dticklabels);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
set(gca,'Color',[.8 .8 .8]);
set(gca,'XLim',[-70 60]);
%xlabel('Latitude \phi');
ylabel('Density \sigma_2 (kg/m^3)');
caxis([-1 1]);
text(-68,dticks(2),'r(-\rho_0f\psi,WS)','FontSize',fontsize);

subplot('Position',axpos(2,:));
[DD,LL] = meshgrid(1:Nd+1,secLats);
% pcolor(LL,DD,1./(2-NSE_map));
pcolor(LL,DD,NSE_map);
shading interp;
colormap(gca,cmocean('amp',10));
colorbar;
% set(gca,'XLim',[ymin ymax_SO]);
set(gca,'YLim',[dticks(1) dticks(end)]);
set(gca,'YTick',dticks);
set(gca,'YTickLabel',dticklabels);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
set(gca,'Color',[.8 .8 .8]);
set(gca,'XLim',[-70 60]);
xlabel('Latitude \phi');
ylabel('Density \sigma_2 (kg/m^3)');
caxis([0 1]);
text(-68,dticks(2),'NSE(-\rho_0f\psi,WS)','FontSize',fontsize);

%%% Add panel labels
for n=[1 2]
  annotation('TextBox',[axpos(n,1)-0.05 axpos(n,2)-0.05 0.03 0.03],'String',axlabels{n},'EdgeColor','None','FontSize',fontsize);
end







%%%
%%% Helper function to calculate correlations and Nash-Sutcliffe efficiencies
%%% for a range of different temporal smoothing time scales.
%%%
function [r,p,rmin,NSE] = calcRelationships (X,Y,smoothlens)
  
  r = 0*smoothlens;
  p = 0*smoothlens;
  rmin = 0*smoothlens;
  NSE = 0*smoothlens;
  for m=1:length(smoothlens)
    
    %%% Smooth time series
    smoothlen = smoothlens(m);
    X_smooth = smooth(X,smoothlen);
    Y_smooth = smooth(Y,smoothlen);
    
    %%% Remove end points where smoothing produces artifices
    idx = ceil(smoothlen/2):length(X)-ceil(smoothlen/2);
    N_smooth = length(idx);
    X_smooth = X_smooth(idx);
    Y_smooth = Y_smooth(idx);
    
    %%% Compute correlation
    [r_tmp,p_tmp] = corr(X_smooth,Y_smooth);
    r(m) = r_tmp;
    p(m) = p_tmp;
    
    %%% Compute effective DOF for X and Y
    lag_max = 3*365;
    Xautocov = zeros(1,lag_max+1);
    Yautocov = zeros(1,lag_max+1);
    Lautocov = zeros(1,lag_max+1);
    for lag = 1:length(Lautocov)      
      Lautocov(lag) = lag-1;
      Xautocov(lag) = corr(X_smooth(1:end-lag+1),X_smooth(lag:end));
      Yautocov(lag) = corr(Y_smooth(1:end-lag+1),Y_smooth(lag:end));
    end
    Te_X = Lautocov(find(Xautocov<1/exp(1),1,'first'));
    Te_Y = Lautocov(find(Xautocov<1/exp(1),1,'first'));    
    if (isempty(Te_X))
      Te_X = lag_max;
      disp(['Covariance limit reached, smoothing = ',num2str(smoothlen)]);
    end
    if (isempty(Te_Y))
      Te_Y = lag_max;
      disp(['Covariance limit reached, smoothing = ',num2str(smoothlen)]);
    end
    Te = max(Te_X,Te_Y);
    
    %%% Compute min. correlation coefficient for statistical significance
    r0 = 0; %%% Initial guess for iteration
    p0 = 0.01; %%% Statistical significance criterion            
    Neff = N_smooth/Te;
    opts1=  optimset('display','off');
    rmin(m) = lsqnonlin(@(r) calcPval(r,Neff)-p0,r0,0,1,opts1);   
    
    %%% Compute NSE
    X_smooth = X_smooth - mean(X_smooth);
    Y_smooth = Y_smooth - mean(Y_smooth);
    NSE(m) = 1 - sum((X_smooth-Y_smooth).^2)/sum((Y_smooth-mean(Y_smooth)).^2);
    
  end
end

%%%
%%% Helper function to perform Gaussian smooting with a specified grid
%%% width w.
%%%
function s = gauss_smooth (r,w)
  
  x = [1:length(r)]';
  s = 0*r;
  for m = 1:length(r)
    filt = exp(-((x-m)/(2*w)).^2);    
    s(m) = sum(r.*filt) / sum(filt);
  end
  
end

%%%
%%% Helper function to compute p-value from correlation coefficient (r) and
%%% number of degrees of freedom (N)
%%%
function p = calcPval (r,N)

  %%% Compute t-statistic
  t=r*sqrt((N-2)/(1-r^2));
  
  %%% Compute p-value  
  p=2*(1-tcdf(t,N-2));
  
end

