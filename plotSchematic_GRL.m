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
fig_mode = 1;

%%% Plotting options
framepos = [417    526   780   730];
labelspacing = 200;
axpos = zeros(1,4);
axpos(1,:) = [0.12 0.1 0.86 0.8];
cmap = [ ...
  14*16+4 15*16+6 15*16+8 ; ...
  8*16+7 12*16+14 15*16+10 ; ...
  ] / 255;
darkblue = [0 0 8*16+11]/255;

%%% Mode-dependent options
switch (fig_mode)
  
  case 1
    titlestr = 'Southern Ocean momentum balance';
    WSlabel = 'Wind stress ($\mathrm{WS}$)';
    IFSlabel = ['Isopycnal',char(10),'form stress ($\mathrm{IFS}$)'];  
    TFSlabel = ['Topographic',char(10),'form stress ($\mathrm{TFS}$)'];
    AABWlabel = ['$|f|\rho_0 {T_{\mathrm{AABW}}} \approx {\mathrm{IFS}} - {\mathrm{TFS}}$'];
    WSlabel_x = -110;
    IFSlabel_x = -80;
    TFSlabel_x = 12;
    WSlabel_y = -500;
    IFSlabel_y = 700;
    TFSlabel_y = 3500;
    AABWlabel_x = -23;
    AABWlabel_y = 2750;
    WSlinewidth = 2;
    IFSlinewidth = 2;
    pfontsize = fontsize+2;
    WSheadwidth = 10;
    figlabel = '(a)';
    IFSlinestyle = '-';
    fontsize = 16;
    
  case 2
    titlestr = 'Multi-annual mean state';
    WSlabel = '$\overline{\mathrm{WS}}$';
    IFSlabel = ['$\overline{\mathrm{IFS}} \approx |f|\rho_0 \overline{T_{\mathrm{AABW}}} + \overline{\mathrm{TFS}}$'];
    TFSlabel = ['$\overline{\mathrm{TFS}} \approx \overline{\mathrm{WS}}$'];
    AABWlabel = ['$\overline{T_{\mathrm{AABW}}}$'];
    WSlabel_x = -97;
    IFSlabel_x = -100;
    TFSlabel_x = 25;
    WSlabel_y = -500;
    IFSlabel_y = 800;
    TFSlabel_y = 3600;
    AABWlabel_x = -10;
    AABWlabel_y = 2800;
    WSlinewidth = 2;
    IFSlinewidth = 2;
    pfontsize = fontsize + 8;
    WSheadwidth = 15;
    figlabel = '(b)';
    IFSlinestyle = '-';
    fontsize = 24;
    
  case 3
    titlestr = 'Response to high-frequency winds';
    WSlabel = '$\mathrm{WS}^\prime$';
    IFSlabel = '$\mathrm{IFS}^\prime \approx 0$';  
    TFSlabel = ['$\mathrm{TFS}^\prime \approx \mathrm{WS}^\prime$'];
    AABWlabel = ['$|f|\rho_0 T_{\mathrm{AABW}}^\prime \approx - \mathrm{WS}^\prime$'];  
    WSlabel_x = -97;
    IFSlabel_x = -60;
    TFSlabel_x = 25;
    WSlabel_y = -500;
    IFSlabel_y = 800;
    TFSlabel_y = 3600;
    AABWlabel_x = -23;
    AABWlabel_y = 2750;
    WSlinewidth = 4;
    IFSlinewidth = 0.5;
    pfontsize = fontsize + 8;
    WSheadwidth = 15;
    figlabel = '(c)';
    IFSlinestyle = '--';
    fontsize = 24;
    
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
set(gca,'XLim',[-110 110]);

%%% Labels
text(IFSlabel_x,IFSlabel_y,IFSlabel,'FontSize',fontsize+2,'interpreter','latex');
text(TFSlabel_x,TFSlabel_y,TFSlabel,'FontSize',fontsize+2,'interpreter','latex');
text(WSlabel_x,WSlabel_y,WSlabel,'FontSize',fontsize+2,'interpreter','latex');
text(AABWlabel_x,AABWlabel_y,AABWlabel,'FontSize',fontsize+2,'interpreter','latex');

if (fig_mode == 1)
  text(-15,800,'$p^+$','FontSize',pfontsize,'interpreter','latex','Color',darkblue);
  text(30,800,'$p^-$','FontSize',pfontsize,'interpreter','latex','Color',darkblue);
  text(62,3000,'$p^+$','FontSize',pfontsize,'interpreter','latex','Color',darkblue);
  text(90,3000,'$p^-$','FontSize',pfontsize,'interpreter','latex','Color',darkblue);
  text(-70,4700,'Sea floor','FontSize',fontsize+2,'interpreter','latex');
  text(-100,1900,'$\sigma_2 = 37.1$ kg/m$^3$','FontSize',fontsize+2,'interpreter','latex');
  line([-51 -67.5],[2276 2050],'Color','k');
  line([-73.5 -60],[4264 4550],'Color','k');
end

%%% Form stress arows
line([-90 -90],[-350 200],'Color',darkblue,'LineWidth',WSlinewidth);
anArrow = annotation('arrow',[0.5 0.5],[0.5 0.5],'Color',darkblue,'LineWidth',WSlinewidth,'HeadStyle','plain','HeadWidth',WSheadwidth);
anArrow.Parent = gca;
anArrow.Position = [-90 200 25 0];
line([-44 -44],[1000 1550],'Color',darkblue,'LineWidth',IFSlinewidth,'LineStyle',IFSlinestyle);
anArrow = annotation('arrow',[0.5 0.5],[0.5 0.5],'Color',darkblue,'LineWidth',IFSlinewidth,'LineStyle',IFSlinestyle,'HeadStyle','plain');
anArrow.Parent = gca;
anArrow.Position = [-44 1550 25 0];
line([50 50],[3800 4350],'Color',darkblue,'LineWidth',WSlinewidth);
anArrow = annotation('arrow',[0.5 0.5],[0.5 0.5],'Color',darkblue,'LineWidth',WSlinewidth,'HeadStyle','plain','HeadWidth',WSheadwidth);
anArrow.Parent = gca;
anArrow.Position = [50 4350 25 0];

%%% Add AABW arrow into the page
if (fig_mode == 3)
  annotation(gcf,'line',[0.57 0.523],[0.553 0.599],'LineWidth',1.5);
  annotation(gcf,'line',[0.523 0.57],[0.553 0.599],'LineWidth',1.5);
  annotation(gcf,'ellipse',[0.5135000000000001 0.544285714285714 0.0659999999999988 0.0642857142857139],'LineWidth',1.5);
else
  annotation(gcf,'line',[0.576 0.517],[0.545 0.604],'LineWidth',1.5);
  annotation(gcf,'line',[0.517 0.576],[0.545 0.604],'LineWidth',1.5);
  annotation(gcf,'ellipse',[0.506000000000001 0.534285714285714 0.0809999999999988 0.0842857142857139],'LineWidth',1.5);
end

%%% Add title
title(titlestr);

%%% Add figure label
text(-130,6350,figlabel,'FontSize',fontsize+2,'interpreter','latex');




