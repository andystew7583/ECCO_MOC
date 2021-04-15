%%%
%%% plotOverturning_GRL.m
%%%
%%% Plots ECCO overturning circulation and variability for our GRL paper.
%%%

%%% Required paths
addpath colormaps;
addpath CDT/cdt;

%%% Load constants
isopDefinitions;

%%% Load streamfunction and mean isopycnal depths
% load(fullfile(products_dir,'PSItot.mat'));
% load(fullfile(products_dir,'Zisop_mean.mat'));
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
pc2fft = fft(pc(2,:)) / std(pc(1,:)) / Nt;
PSIfft = fft(PSImean) / Nt;
pc1fft(1) = 0;
pc2fft(1) = 0;
PSIfft(1) = 0;
TT = Nt ./ (0:1:Nt/2-1);

%%% Plotting options
fontsize = 14;
framepos = [417    526   791   800];
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
handle = figure(201);
clf;
set(handle,'Position',framepos);





%%% Pretty version of the adjusted residual streamfunction in density space
subplot('Position',axpos(1,:));
LL = repmat([lat]',[1 Nd+1]);
pcolor(LL,-Zisop_mean,mean(PSI,3)/1e6);
shading interp;
colormap(gca,cmocean('balance',40));
% colormap(redblue(40));
cbhandle = colorbar;
set(cbhandle,'Position',cbpos);
caxis([-psimax psimax]);
axis([ymin ymax_GO 0 6000]);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
set(gca,'Color',[.8 .8 .8]);
% xlabel('Latitude');
ylabel('Depth (m)');
% title('Residual overturning streamfunction in depth coordinates');

%%% Add isopycnals
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
end
jmin = find(lat>=-60,1,'first');
jmax = find(lat<=-50,1,'last');
plot(lat(jmin:jmax),-Zisop_plot(jmin:jmax,dval),'w--','LineWidth',1.5);
plot(lat([jmin jmax]),-Zisop_plot([jmin jmax],dval),'wo');
hold off;

text(22,430,num2str(dens_levs(Zisop_idx(1))-1000),'FontSize',fontsize,'Color',[.5 .5 .5]);
text(47,870,num2str(dens_levs(Zisop_idx(2))-1000),'FontSize',fontsize,'Color',[.5 .5 .5]);
text(50,2870,num2str(dens_levs(Zisop_idx(3))-1000),'FontSize',fontsize,'Color',[.5 .5 .5]);
text(44,5100,num2str(dens_levs(Zisop_idx(4))-1000),'FontSize',fontsize,'Color',[.5 .5 .5]);
text(-40,4900,num2str(dens_levs(Zisop_idx(5))-1000),'FontSize',fontsize,'Color',[.5 .5 .5]);

anArrow = annotation(gcf,'arrow',[0 0],[0 0]);
anArrow.Parent = gca;
anArrow.Position = [-20 5000 20 0];
anArrow = annotation(gcf,'arrow',[0 0],[0 0]);
anArrow.Parent = gca;
anArrow.Position = [0 3000 -20 0];
anArrow = annotation(gcf,'arrow',[0 0],[0 0]);
anArrow.Parent = gca;
anArrow.Position = [0 2200 -20 0];
anArrow = annotation(gcf,'arrow',[0 0],[0 0]);
anArrow.Parent = gca;
anArrow.Position = [-20 500 20 0];
% anArrow = annotation(gcf,'arrow',[0 0],[0 0]);
% anArrow.Parent = gca;
% anArrow.Position = [-48 50 10 100];
% anArrow = annotation(gcf,'arrow',[0 0],[0 0]);
% anArrow.Parent = gca;
% anArrow.Position = [-48 2000 -10 -1300];
% anArrow = annotation(gcf,'arrow',[0 0],[0 0]);
% anArrow.Parent = gca;
% anArrow.Position = [-48 2750 -10 2000];
% anArrow = annotation(gcf,'arrow',[0 0],[0 0]);
% anArrow.Parent = gca;
% anArrow.Position = [-66 2800 10 2000];

%%% Add arrows indicating directions of circulation
%%% TODO FIX
% figure1 = gcf;
% annotation(figure1,'arrow',[0.802 0.7],[0.676 0.648]);
% annotation(figure1,'arrow',[0.262 0.203],[0.684 0.804000000000004]);
% annotation(figure1,'arrow',[0.261000000000001 0.32],...
%   [0.878000000000001 0.85]);
% annotation(figure1,'arrow',[0.687 0.777],[0.878 0.916]);
% annotation(figure1,'arrow',[0.17 0.217],[0.61 0.452]);
% annotation(figure1,'arrow',[0.221 0.178],[0.604 0.776]);
% annotation(figure1,'arrow',[0.618000000000001 0.504],...
%   [0.482000000000001 0.492]);
% annotation(figure1,'arrow',[0.491000000000001 0.599],...
%   [0.248000000000001 0.264]);








%%% Southern Ocean overturning in density space
subplot('Position',axpos(2,:));
% [DD,LL] = meshgrid(dens_bnds,lat);
[DD,LL] = meshgrid(1:Nd+1,lat);
% pcolor(LL,DD-1000,mean(PSI,3)/1e6);
pcolor(LL,DD,mean(PSI,3)/1e6);
shading interp;
hold on;
plot([-60 -50],[dval dval],'w--o','LineWidth',1.5);
hold off;
colormap(gca,cmocean('balance',40));
% colorbar;
caxis([-psimax psimax]);
set(gca,'XLim',[ymin ymax_SO]);
set(gca,'YLim',[dticks(1) dticks(end)]);
set(gca,'YTick',dticks);
set(gca,'YTickLabel',dticklabels);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
set(gca,'Color',[.8 .8 .8]);
xlabel('Latitude \phi');
ylabel('Density \sigma_2 (kg/m^3)');
text(-78,dticks(2),['\psi(\phi,\sigma_2)'],'FontSize',fontsize);

%%% EOF1
subplot('Position',axpos(3,:));
[DD,LL] = meshgrid(1:Nd+1,lat);
pcolor(LL,DD,eof_maps(:,:,1)*std(pc(1,:))/1e6);
shading interp;
hold on;
plot([-60 -50],[dval dval],'w--o','LineWidth',1.5);
hold off;
colormap(gca,cmocean('balance',40));
% colorbar;
caxis([-psimax psimax]);
set(gca,'XLim',[ymin ymax_SO]);
set(gca,'YLim',[dticks(1) dticks(end)]);
set(gca,'YTick',dticks);
set(gca,'YTickLabel',dticklabels);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
set(gca,'Color',[.8 .8 .8]);
xlabel('Latitude \phi');
text(-77,dticks(2),['EOF1'],'FontSize',fontsize);
text(-77,dticks(3),['Exp. var.: ',num2str(expvar(1),'%.0f'),'%'],'FontSize',fontsize);
% ylabel('Density \sigma_2 (kg/m^3)');

%%% EOF2
subplot('Position',axpos(4,:));
[DD,LL] = meshgrid(1:Nd+1,lat);
pcolor(LL,DD,-eof_maps(:,:,2)*std(pc(1,:))/1e6);
shading interp;
hold on;
plot([-60 -50],[dval dval],'w--o','LineWidth',1.5);
hold off;
colormap(gca,cmocean('balance',40));
caxis([-psimax psimax]);
set(gca,'XLim',[ymin ymax_SO]);
set(gca,'YLim',[dticks(1) dticks(end)]);
set(gca,'YTick',dticks);
set(gca,'YTickLabel',dticklabels);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
set(gca,'Color',[.8 .8 .8]);
xlabel('Latitude \phi');
text(-77,dticks(2),['EOF2'],'FontSize',fontsize);
text(-77,dticks(3),['Exp. var.: ',num2str(expvar(2),'%.0f'),'%'],'FontSize',fontsize);
% ylabel('Density \sigma_2 (kg/m^3)');

%%% Spectra
subplot('Position',axpos(5,:));
colororder = get(gca,'ColorOrder');
p1 = loglog(TT,abs(PSIfft(1:Nt/2))/1e6);
hold on
p2 = loglog(TT,abs(pc1fft(1:Nt/2)));
p3 = loglog(TT,abs(pc2fft(1:Nt/2)),'Color',colororder(4,:));
hold off
axis tight
xlabel('Period (days)');
ylabel('Spectral power');
set(gca,'XTick',([3 10 30 100 300 1000 3000]));
set(gca,'FontSize',fontsize);
% p2.Color(4) = 0.25;
p3.Color(4) = 0.5;
handle = legend('Abyssal overturning $T_{\mathrm{AABW}}$ (Sv)','1st principal component','2nd principal component','Location','SouthEast');;

set(handle,'interpreter','latex');

for n=1:size(axpos,1)
  annotation('TextBox',[axpos(n,1)-0.05 axpos(n,2)-0.05 0.03 0.03],'String',axlabels{n},'EdgeColor','None','FontSize',fontsize);
end