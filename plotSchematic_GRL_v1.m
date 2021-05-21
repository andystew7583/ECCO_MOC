%%%
%%% plotSchematic_GRL.m6
%%%
%%% Plots a schematic of the isopycnal momentum balance for our GRL paper.
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load global variables
isopDefinitions;

%%% Add required paths
p = genpath('gcmfaces/'); addpath(p);
p = genpath('gcmfaces/'); addpath(p);
p = genpath('m_map/'); addpath(p);
p = genpath('rw_namelist/'); addpath(p);

%%% Load all grid variables from nctiles_grid/ into mygrid
grid_load([ECCO_grid_dir filesep],5,'nctiles',0,1);

%%% Make mygrid accessible in current workspace
gcmfaces_global;

%%% Directory from which to read ECCO output data
myenv.nctilesdir = fullfile(ECCO_data_dir,filesep);

%%% Latitude along which to plot the section
lat_plot = -58;

%%% Density surface to plot for IFS
dens_plot = 37.1;



%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%

%%% Load mean density
load(fullfile(products_dir,'DENS_mean.mat'));
DENS = convert2array(DENS_mean);


%%% Grids
XC = convert2array(mygrid.XC);
YC = convert2array(mygrid.YC);
RC = mygrid.RC;
Depth = convert2array(mygrid.Depth);
hFacC = convert2array(mygrid.hFacC);

%%% Find target latitude index
secLats = mean(YC,1);
latidx = find(secLats>=lat_plot,1);

%%% Remove land points and extract density dection
DENS(hFacC==0) = 1037.2;
DENSxz = squeeze(DENS(:,latidx,:));
secDepth = Depth(:,latidx);

%%% Rearrange matrices so that x-grid is monotonically increasing
xx = XC(:,latidx)';
idlidx = find(xx>0,1,'last')+1;
sortidx = [idlidx:length(xx) 1:idlidx-1];
xx = xx(sortidx);
DENSxz = DENSxz(sortidx,:);
secDepth = secDepth(sortidx)';








%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAKE THE PLOT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set true to show perturbed state, false to show mean state
show_pert = true;

%%% Plotting options
fontsize = 16;
framepos = [417    526   1000   700];
labelspacing = 200;
axpos = zeros(1,4);
axpos(1,:) = [0.08 0.1 0.9 0.8];
cmap = [ ...
  14*16+4 15*16+6 15*16+8 ; ...
  8*16+7 12*16+14 15*16+10 ; ...
  ] / 255;
darkblue = [0 0 8*16+11]/255;
if (show_pert)
  titlestr = 'Response to high-frequency wind fluctuations';
  WSlabel = '$\mathrm{WS}^\prime$';
  IFSlabel = '$\mathrm{IFS}^\prime \approx 0$';  
  TFSlabel = ['$\mathrm{TFS}^\prime \approx \mathrm{WS}^\prime$'];
  AABWlabel = ['$|f|\rho_0 T_{\mathrm{AABW}}^\prime \approx - \mathrm{WS}^\prime$'];  
  WSlabel_x = -117;
  IFSlabel_x = -60;
  TFSlabel_x = 30;
  WSlabel_y = -500;
  IFSlabel_y = 800;
  TFSlabel_y = 3600;
  WSlinewidth = 4;
  IFSlinewidth = 0.5;
  pfontsize = fontsize + 8;
  WSheadwidth = 15;
  figlabel = '(b)';
  IFSlinestyle = '--';
else
  titlestr = 'Multi-decadal mean state';
  WSlabel = 'Wind stress ($\overline{\mathrm{WS}}$)';
  IFSlabel = ['Isopycnal',char(10),'form stress ($\overline{\mathrm{IFS}}$)'];  
  TFSlabel = ['Topographic',char(10),'form stress ($\overline{\mathrm{TFS}}$)'];
  AABWlabel = ['$|f|\rho_0 {T_{\mathrm{AABW}}} \approx {\mathrm{IFS}} - {\mathrm{TFS}}$'];
  WSlabel_x = -130;
  IFSlabel_x = -80;
  TFSlabel_x = 12;
  WSlabel_y = -500;
  IFSlabel_y = 700;
  TFSlabel_y = 3500;
  WSlinewidth = 2;
  IFSlinewidth = 2;
  pfontsize = fontsize+2;
  WSheadwidth = 10;
  figlabel = '(a)';
  IFSlinestyle = '-';
end

%%% cnoidal waveform for the surface
M = 0.99;
[S,C,D] = ellipj(xx/2,M);

%%% Set up figure window
handle = figure(202);
clf;
set(handle,'Position',framepos);
[ZZ,XX] = meshgrid(RC,xx);
contourf(XX,-ZZ,squeeze(DENSxz)-1000,[0 37.1 100]);
hold on;
% [C,h] = contour(XX,-ZZ,DENSxz-1000,[36 36.5 37 37.05 37.15 37.2],'EdgeColor',[.85 .85 .85]);
% clabel(C,h);
ahandle = area(xx',[secDepth; -RC(end)-secDepth]');
ahandle(1).FaceAlpha = 0;
ahandle(2).FaceColor = [12*16+4 10*16+4 8*16+4]/255;
ahandle = area(xx',-150*D');
ahandle.EdgeColor = cmap(1,:);
% ahandle(1).FaceAlpha = 0;
ahandle.FaceColor = cmap(1,:);
plot(xx,0*xx,'Color',cmap(1,:));
hold off;
set(gca,'YDir','reverse');
colormap(cmap);
box off;
set(gca,'FontSize',fontsize);
xlabel('Longitude');
ylabel('Depth (m)')
set(gca,'Position',axpos);

%%% Labels
text(-15,800,'$p^+$','FontSize',pfontsize,'interpreter','latex','Color',darkblue);
text(30,800,'$p^-$','FontSize',pfontsize,'interpreter','latex','Color',darkblue);
text(62,3000,'$p^+$','FontSize',pfontsize,'interpreter','latex','Color',darkblue);
text(90,3000,'$p^-$','FontSize',pfontsize,'interpreter','latex','Color',darkblue);
text(IFSlabel_x,IFSlabel_y,IFSlabel,'FontSize',fontsize+2,'interpreter','latex');
text(TFSlabel_x,TFSlabel_y,TFSlabel,'FontSize',fontsize+2,'interpreter','latex');
text(WSlabel_x,WSlabel_y,WSlabel,'FontSize',fontsize+2,'interpreter','latex');
text(-23,2750,AABWlabel,'FontSize',fontsize+2,'interpreter','latex');

if (~show_pert)
  text(-150,5000,'Sea floor','FontSize',fontsize+2,'interpreter','latex');
  text(100,1500,'$\sigma_2 = 37.1$ kg/m$^3$','FontSize',fontsize+2,'interpreter','latex');
  line([103.5 113.5],[1914.15 1634.18],'Color','k');
  line([-135.5 -123.5],[4800 4264.25],'Color','k');
end

%%% Form stress arows
line([-110 -110],[-350 200],'Color',darkblue,'LineWidth',WSlinewidth);
anArrow = annotation('arrow',[0.5 0.5],[0.5 0.5],'Color',darkblue,'LineWidth',WSlinewidth,'HeadStyle','plain','HeadWidth',WSheadwidth);
anArrow.Parent = gca;
anArrow.Position = [-110 200 25 0];
line([-44 -44],[1000 1550],'Color',darkblue,'LineWidth',IFSlinewidth,'LineStyle',IFSlinestyle);
anArrow = annotation('arrow',[0.5 0.5],[0.5 0.5],'Color',darkblue,'LineWidth',IFSlinewidth,'LineStyle',IFSlinestyle,'HeadStyle','plain');
anArrow.Parent = gca;
anArrow.Position = [-44 1550 25 0];
line([50 50],[3800 4350],'Color',darkblue,'LineWidth',WSlinewidth);
anArrow = annotation('arrow',[0.5 0.5],[0.5 0.5],'Color',darkblue,'LineWidth',WSlinewidth,'HeadStyle','plain','HeadWidth',WSheadwidth);
anArrow.Parent = gca;
anArrow.Position = [50 4350 25 0];

%%% Add AABW arrow into the page
if (show_pert)
  annotation(gcf,'line',[0.563 0.53],[0.553 0.599],'LineWidth',1.5);
  annotation(gcf,'line',[0.53 0.563],[0.553 0.599],'LineWidth',1.5);
  annotation(gcf,'ellipse',[0.5235000000000001 0.544285714285714 0.0459999999999988 0.0642857142857139],'LineWidth',1.5);
else
  annotation(gcf,'line',[0.568 0.525],[0.545 0.604],'LineWidth',1.5);
  annotation(gcf,'line',[0.525 0.568],[0.545 0.604],'LineWidth',1.5);
  annotation(gcf,'ellipse',[0.516000000000001 0.534285714285714 0.0609999999999988 0.0842857142857139],'LineWidth',1.5);
end

%%% Add title
title(titlestr);

%%% Add figure label
text(-200,6350,figlabel,'FontSize',fontsize+2,'interpreter','latex');




