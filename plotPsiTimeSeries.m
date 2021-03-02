%%% 
%%% plotPsiTimeSeries.m
%%%
%%% Plots a time series of the MOC calculated end-to-end across the Southern Ocean.
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD STREAMFUNCTION AND CALCULATE TIME SERIES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Isopycnal depths
load([products_dir filesep 'Zisop.mat']);
Nlats = length(lat);

%%% Total (i.e. residual) streamfunction
load([products_dir filesep 'PSItot.mat']);

%%% Select streamfunction to plot
psitype = 'tot';
eval(['PSI = PSI',psitype,';']);
eval(['PSI_mean = PSI',psitype,'_mean;']);
eval(['PSI_zonmean = PSI',psitype,'_zonmean;']);

%%% Adjusted streamfunction (enforces psi=0 at sea floor)
% PSI_adj = PSI-repmat(PSI(:,end,:),[1 Nd+1 1]);
PSI_adj = PSI;

%%% Calculate strength of AABW cell across the ACC
ACC_idx = (lat>-60) & (lat<-40);
PSI_aabw = zeros(1,Nt);
for n=1:Nt 
  PSI_aabw(n) = -max(min(PSI_adj(ACC_idx,:,n),[],2),[],1);
end

Nmonths = 12;
Nyears = Nt/Nmonths;
PSI_aabw_yearly = zeros(1,Nyears);
PSI_yearly = reshape(PSI_adj,[size(PSI_adj,1) size(PSI_adj,2) Nmonths Nyears]);
PSI_yearly = squeeze(mean(PSI_yearly,3));
for n=1:Nyears
  PSI_aabw_yearly(n) = -max(min(PSI_yearly(ACC_idx,:,n),[],2),[],1);
end





%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CREATE PLOTS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plotting options
fontsize = 14;
scrsz = get(0,'ScreenSize');
framepos = [0.25*scrsz(3) 0.25*scrsz(4) 1000 500];
fignum = 0;

%%% Time series of AABW cell strength
fignum = fignum + 1;
handle = figure(fignum);
clf;
set(handle,'Position',framepos);
plot(1:Nt,PSI_aabw);
set(gca,'FontSize',fontsize);
xlabel('Time (months)');
ylabel('\Psi_A_A_B_W, monthly (Sv)');

%%% Time series of AABW cell strength
fignum = fignum + 1;
handle = figure(fignum);
clf;
set(handle,'Position',framepos);
plot(1:Nyears,PSI_aabw_yearly);
set(gca,'FontSize',fontsize);
xlabel('Time (years)');
ylabel('\Psi_A_A_B_W, annual-mean (Sv)');