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
load(fullfile(products_dir,'PSItot.mat'));
load(fullfile(products_dir,'Zisop_mean.mat'));
Nlats = length(lat);

%%% Compute EOFs
eof_msk = (mean(PSI,3)~=0);
eof_idx = find(lat>-40);
eof_msk(eof_idx,:) = 0;
[eof_maps,pc,expvar] = eof(PSI,10,'mask',eof_msk);
eof_maps(isnan(eof_maps)) = 0;

%%% Plotting options
fontsize = 14;
framepos = [417    526   791   800];
labelspacing = 200;
axpos = zeros(4,4);
axpos(1,:) = [0.08 0.55 0.85 0.4];
axpos(2,:) = [0.08 0.05 0.25 0.4];
axpos(3,:) = [0.38 0.05 0.25 0.4];
axpos(4,:) = [0.68 0.05 0.25 0.4];
cbpos = [0.95 0.05 0.015 0.9];
drange = [34 37.5];
ymin = -78;
ymax_SO = -40;
ymax_GO = 78;
psimax = 20;

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
hold off;

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







dticks = [36 43 53 63 75 87 99 109 119 129];
dticklabels = {'34','35','35.5','36','36.3','36.6','36.9','37','37.1','37.2'};

%%% Southern Ocean overturning in density space
subplot('Position',axpos(2,:));
% [DD,LL] = meshgrid(dens_bnds,lat);
[DD,LL] = meshgrid(1:Nd+1,lat);
% pcolor(LL,DD-1000,mean(PSI,3)/1e6);
pcolor(LL,DD,mean(PSI,3)/1e6);
shading interp;
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
xlabel('Latitude');
ylabel('Density \sigma_2 (kg/m^3)');

%%% EOF1
subplot('Position',axpos(3,:));
[DD,LL] = meshgrid(1:Nd+1,lat);
pcolor(LL,DD,eof_maps(:,:,1)*std(pc(1,:))/1e6);
shading interp;
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
xlabel('Latitude');
% ylabel('Density \sigma_2 (kg/m^3)');

%%% EOF2
subplot('Position',axpos(4,:));
[DD,LL] = meshgrid(1:Nd+1,lat);
pcolor(LL,DD,eof_maps(:,:,2)*std(pc(1,:))/1e6);
shading interp;
colormap(gca,cmocean('balance',40));
caxis([-psimax psimax]);
set(gca,'XLim',[ymin ymax_SO]);
set(gca,'YLim',[dticks(1) dticks(end)]);
set(gca,'YTick',dticks);
set(gca,'YTickLabel',dticklabels);
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
set(gca,'Color',[.8 .8 .8]);
xlabel('Latitude');
% ylabel('Density \sigma_2 (kg/m^3)');
