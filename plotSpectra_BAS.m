%%%
%%% plotSpectra_BAS.m
%%%
%%% Plots ECCO overturning circulation and variability for our BAS seminar.
%%%

%%% Required paths
addpath colormaps;
addpath CDT/cdt;

%%% Load constants
isopDefinitions;

%%% Load streamfunction and mean isopycnal depths
load(fullfile(products_dir,'PSItot.mat'));
load(fullfile(products_dir,'Zisop_mean.mat'));
Nlats = length(lat);

%%% Compute AABW transport
ymin = -60;
ymax = -50;
dens_psimax = 1037.1;
didx_psimax = find(dens_bnds>dens_psimax,1,'first');
yidx_psi = find((lat>ymin) & (lat<ymax));
PSImean = squeeze(mean(PSI(yidx_psi,didx_psimax,:),1));

%%% Compute EOFs
eof_msk = (mean(PSI,3)~=0);
eof_idx = find(lat>-40);
eof_msk(eof_idx,:) = 0;
[eof_maps,pc,expvar] = eof(PSI,10,'mask',eof_msk);
eof_maps(isnan(eof_maps)) = 0;

%%% Compute FFTs
pc1fft = fft(pc(1,:)) / std(pc(1,:)) / Nt;
pc2fft = fft(pc(2,:)) / std(pc(2,:)) / Nt;
pc3fft = fft(pc(3,:)) / std(pc(3,:)) / Nt;
PSIfft = fft(PSImean) / Nt;
pc1fft(1) = 0;
pc2fft(1) = 0;
pc3fft(1) = 0;
PSIfft(1) = 0;
TT = Nt ./ (0:1:Nt/2-1);

%%% Wind time series
load(fullfile(products_dir,'WS.mat'));
ymin = -60;
ymax = -50;
yidx_ifs = find((secLats>ymin) & (secLats<ymax));
WSmean = squeeze(mean(WS(yidx_ifs,:)))';
WSfft = fft(WSmean)/Nt;
WSfft(1) = 0;

%%% Plotting options
fontsize = 14;
framepos = [417    526   1200   600];
labelspacing = 200;
axpos = zeros(5,4);
axpos(1,:) = [0.08 0.7 0.85 0.27];
axpos(2,:) = [0.08 0.37 0.24 0.27];
axpos(3,:) = [0.38 0.37 0.24 0.27];
axpos(4,:) = [0.68 0.37 0.24 0.27];
axpos(5,:) = [0.08 0.07 0.85 0.23];
cbpos = [0.95 0.37 0.015 0.6];
drange = [34 37.5];
ymin = -78;
ymax_SO = -40;
ymax_GO = 78;
psimax = 20;
dticks = [36 43 53 63 75 87 99 109 119 129];
dticklabels = {'34','35','35.5','36','36.3','36.6','36.9','37','37.1','37.2'};
dval = find(dens_levs==1037.1); %%% Index at which the density is equal to our AABW density threshold
axlabels = {'(a)','(b)','(c)','(d)','(e)'};

%%% Set up figure window
handle = figure(220);
clf;
set(handle,'Position',framepos);


%%% Spectra
axes('Position',[0.1 0.1 0.8 0.8]);
colororder = get(gca,'ColorOrder');
p1 = semilogx(TT,1-2*cumsum(abs(PSIfft(1:Nt/2).^2))/var(PSImean),'LineWidth',1.5);
hold on
p2 = semilogx(TT,(1-2*cumsum(abs(pc1fft(1:Nt/2).^2)))*expvar(1)/100,'LineWidth',1.5); % loglog(TT,abs(pc1fft(1:Nt/2)));
p3 = semilogx(TT,(1-2*cumsum(abs(pc2fft(1:Nt/2).^2)))*expvar(2)/100,'Color',colororder(4,:),'LineWidth',1.5); % loglog(TT,abs(pc2fft(1:Nt/2)),'Color',colororder(4,:));
p4 = semilogx(TT,(1-2*cumsum(abs(pc3fft(1:Nt/2).^2)))*expvar(3)/100,'Color',colororder(5,:),'LineWidth',1.5); % loglog(TT,abs(pc2fft(1:Nt/2)),'Color',colororder(4,:));
p5 = semilogx(TT,(1-2*cumsum(abs(WSfft(1:Nt/2).^2))/var(WSmean)),'LineStyle','--','Color','k','LineWidth',1.5); % loglog(TT,abs(pc2fft(1:Nt/2)),'Color',colororder(4,:));
hold off
axis tight
xlabel('Period (days)');
ylabel('Cumulative variance');
set(gca,'XTick',([3 10 30 100 300 1000 3000]));
set(gca,'FontSize',fontsize);
% p2.Color(4) = 0.25;
p3.Color(4) = 0.5;
handle = legend('Abyssal overturning $T_{\mathrm{AABW}}$ (Sv)','1st principal component','2nd principal component','3rd principal component','Wind stress','Location','NorthWest');
set(handle,'interpreter','latex');
grid on;
% text(1.5,0.9,['\psi(\phi,z)'],'FontSize',fontsize);






