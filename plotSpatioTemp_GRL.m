%%%
%%% plotSpatioTemp_GRL.m
%%% 
%%% Plots spatio-temporal relationships between

%%% Required paths
addpath shadedErrorBar;

%%% Load constants
isopDefinitions;

%%% Time index at which to switch from short timescale to long timescale
%%% plot
endIdx = 18;

%%% Set true to smooth long timescale plot
do_smooth = true;

%%% Length of box filter to apply to long timescale plot
smoothLen = 30;

%%% Load pre-computed RFs
load('RFs.mat');

%%% Load wind stress time series
load(fullfile(products_dir,'WS.mat'));
ymin = -60;
ymax = -50;
yidx_ifs = find((secLats>ymin) & (secLats<ymax));
WSmean = squeeze(mean(WS(yidx_ifs,:)))';

%%% Cumulative responses
PSI_WS_CRF = cumsum(PSI_WS_RF,1);
TFS_WS_CRF = cumsum(TFS_WS_RF,1);
IFS_WS_CRF = cumsum(IFS_WS_RF,1);
EFS_WS_CRF = cumsum(EFS_WS_RF,1);
PSI_WS_err = sqrt(std(PSI_WS_CRF,[],2).^2 + (std(PSI_WS_eps(:))/std(WSmean)).^2);
TFS_WS_err = sqrt(std(TFS_WS_CRF,[],2).^2 + (std(TFS_WS_eps(:))/std(WSmean)).^2);
IFS_WS_err = sqrt(std(IFS_WS_CRF,[],2).^2 + (std(IFS_WS_eps(:))/std(WSmean)).^2);
EFS_WS_err = sqrt(std(EFS_WS_CRF,[],2).^2 + (std(EFS_WS_eps(:))/std(WSmean)).^2);

%%% Smoothed responses
PSI_WS_CRF_smooth = PSI_WS_CRF;
TFS_WS_CRF_smooth = TFS_WS_CRF;
IFS_WS_CRF_smooth = IFS_WS_CRF;
EFS_WS_CRF_smooth = EFS_WS_CRF;
if (do_smooth)
  for m = 1:Nwin
    PSI_WS_CRF_smooth(:,m) = smooth(PSI_WS_CRF(:,m),smoothLen);
    TFS_WS_CRF_smooth(:,m) = smooth(TFS_WS_CRF(:,m),smoothLen);
    IFS_WS_CRF_smooth(:,m) = smooth(IFS_WS_CRF(:,m),smoothLen);
    EFS_WS_CRF_smooth(:,m) = smooth(EFS_WS_CRF(:,m),smoothLen);
  end
end
PSI_WS_err_smooth = sqrt(std(PSI_WS_CRF_smooth,[],2).^2 + (std(PSI_WS_eps(:))/std(WSmean)).^2);
TFS_WS_err_smooth = sqrt(std(TFS_WS_CRF_smooth,[],2).^2 + (std(TFS_WS_eps(:))/std(WSmean)).^2);
IFS_WS_err_smooth = sqrt(std(IFS_WS_CRF_smooth,[],2).^2 + (std(IFS_WS_eps(:))/std(WSmean)).^2);
EFS_WS_err_smooth = sqrt(std(EFS_WS_CRF_smooth,[],2).^2 + (std(EFS_WS_eps(:))/std(WSmean)).^2);
  

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE THE PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
fontsize = 14;
framepos = [417    526   791   600];
labelspacing = 200;
axpos = zeros(6,4);
axpos(1,:) = [0.08 0.1 0.29 0.85];
axpos(2,:) = [0.38 0.1 0.6 0.85];
ylim = [-1.4 1.4];


figure(204);                
clf;
set(gcf,'Position',framepos);   

%%% Short timescale response
subplot('Position',axpos(1,:));
colororder = get(gca,'ColorOrder');
tt_plot = tt_rf(1:endIdx);
% plot(tt_plot,mean(PSI_WS_CRF(1:endIdx,:),2));
% shadedErrorBar(tt_plot,mean(PSI_WS_CRF(1:endIdx,:),2),std(PSI_WS_CRF(1:endIdx,:),[],2),'lineProps',{'Color',colororder(1,:)});
shadedErrorBar(tt_plot,mean(PSI_WS_CRF(1:endIdx,:),2),PSI_WS_err(1:endIdx,:),'lineProps',{'Color',colororder(1,:)});
hold on
% plot(tt_plot,mean(TFS_WS_CRF(1:endIdx,:),2));
% plot(tt_plot,mean(IFS_WS_CRF(1:endIdx,:),2));
% plot(tt_plot,mean(EFS_WS_CRF(1:endIdx,:),2));
% shadedErrorBar(tt_plot,mean(TFS_WS_CRF(1:endIdx,:),2),std(TFS_WS_CRF(1:endIdx,:),[],2),'lineProps',{'Color',colororder(2,:)});
% shadedErrorBar(tt_plot,mean(IFS_WS_CRF(1:endIdx,:),2),std(IFS_WS_CRF(1:endIdx,:),[],2),'lineProps',{'Color',colororder(3,:)});
% shadedErrorBar(tt_plot,mean(EFS_WS_CRF(1:endIdx,:),2),std(EFS_WS_CRF(1:endIdx,:),[],2),'lineProps',{'Color',colororder(4,:)});
shadedErrorBar(tt_plot,-mean(TFS_WS_CRF(1:endIdx,:),2),TFS_WS_err(1:endIdx,:),'lineProps',{'Color',colororder(2,:)});
shadedErrorBar(tt_plot,-mean(IFS_WS_CRF(1:endIdx,:),2),IFS_WS_err(1:endIdx,:),'lineProps',{'Color',colororder(3,:)});
shadedErrorBar(tt_plot,-mean(EFS_WS_CRF(1:endIdx,:),2),EFS_WS_err(1:endIdx,:),'lineProps',{'Color',colororder(4,:)});
plot([-1 tt_plot],0*[-1 tt_plot],'k--');
plot([0 0],ylim,'k:');
hold off;
set(gca,'XLim',[-1 endIdx-1]);
set(gca,'YLim',ylim);
xlabel('Lag time \tau (days)');
ylabel('Normalized response to a step change in wind stress');
set(gca,'FontSize',fontsize);
box on;

%%% Long timescale response
subplot('Position',axpos(2,:));
tt_plot = tt_rf(endIdx+1:end)/365;
% plot(tt_plot,mean(PSI_WS_CRF_smooth(endIdx+1:end,:),2));
% shadedErrorBar(tt_plot,mean(PSI_WS_CRF_smooth(endIdx+1:end,:),2),std(PSI_WS_CRF_smooth(endIdx+1:end,:),[],2),'lineProps',{'Color',colororder(1,:)});
shadedErrorBar(tt_plot,mean(PSI_WS_CRF_smooth(endIdx+1:end,:),2),PSI_WS_err_smooth(endIdx+1:end,:),'lineProps',{'Color',colororder(1,:)});
hold on;
% plot(tt_plot,mean(TFS_WS_CRF_smooth(endIdx+1:end,:),2));
% plot(tt_plot,mean(IFS_WS_CRF_smooth(endIdx+1:end,:),2));
% plot(tt_plot,mean(EFS_WS_CRF_smooth(endIdx+1:end,:),2));
% shadedErrorBar(tt_plot,mean(TFS_WS_CRF_smooth(endIdx+1:end,:),2),std(TFS_WS_CRF_smooth(endIdx+1:end,:),[],2),'lineProps',{'Color',colororder(2,:)});
% shadedErrorBar(tt_plot,mean(IFS_WS_CRF_smooth(endIdx+1:end,:),2),std(IFS_WS_CRF_smooth(endIdx+1:end,:),[],2),'lineProps',{'Color',colororder(3,:)});
% shadedErrorBar(tt_plot,mean(EFS_WS_CRF_smooth(endIdx+1:end,:),2),std(EFS_WS_CRF_smooth(endIdx+1:end,:),[],2),'lineProps',{'Color',colororder(4,:)});
shadedErrorBar(tt_plot,-mean(TFS_WS_CRF_smooth(endIdx+1:end,:),2),TFS_WS_err_smooth(endIdx+1:end,:),'lineProps',{'Color',colororder(2,:)});
shadedErrorBar(tt_plot,-mean(IFS_WS_CRF_smooth(endIdx+1:end,:),2),IFS_WS_err_smooth(endIdx+1:end,:),'lineProps',{'Color',colororder(3,:)});
shadedErrorBar(tt_plot,-mean(EFS_WS_CRF_smooth(endIdx+1:end,:),2),EFS_WS_err_smooth(endIdx+1:end,:),'lineProps',{'Color',colororder(4,:)});
plot(tt_plot,0*tt_plot,'k--');
hold off;
set(gca,'YLim',ylim);
set(gca,'XLim',[tt_plot(1) tt_plot(end)]);
set(gca,'YTickLabel','');
xlabel('Lag time \tau (years)');
set(gca,'FontSize',fontsize);
box on;
handle = legend('Northward AABW transport ($T_{\mathrm{AABW}}$)','Topographic form stress (TFS)','Resolved isopycnal form stress (RIFS$_\mathrm{AABW}$)','Eddy isopycnal form stress (EIFS$_\mathrm{AABW}$)','Location','SouthEast');
set(handle,'interpreter','latex');
set(handle,'FontSize',fontsize+2);

