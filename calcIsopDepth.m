%%%
%%% calcIsopDepth.m
%%%
%%% Calculates depths of density surfaces from each ECCO output snapshot.
%%% Must be called after calcDensity.
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Start by clearing memory
clear all;

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

%%% Tools for calculating properties in latitude bands
if ~isfield(mygrid,'LATS_MASKS')  
  gcmfaces_lines_zonal;   
end

%%% Grids
Nlats = length([mygrid.LATS_MASKS.lat]);
Nz = length(mygrid.RC);
drf = mk3D(mygrid.DRF,mygrid.hFacC);
dxg = mk3D(mygrid.DXG,mygrid.hFacC);
dyg = mk3D(mygrid.DYG,mygrid.hFacC);
rf = mygrid.RF;

%%% This is a hack to make the calc_MeridionalTransport function calculate
%%% areas properly
for j = 1:Nlats
  mygrid.LATS_MASKS(j).mskWedge = abs(mygrid.LATS_MASKS(j).mskWedge);
  mygrid.LATS_MASKS(j).mskSedge = abs(mygrid.LATS_MASKS(j).mskSedge);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATE AREA ABOVE GEOPOTENTIALS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% To store computed depths
zisop = zeros(Nlats,Nd,Nt);
zisop_mean = zeros(Nlats,Nd);

%%% Calculate cell face areas 
AREAW_z = mygrid.mskW .* mygrid.hFacW .* dyg .* drf;
AREAS_z = mygrid.mskS .* mygrid.hFacS .* dxg .* drf;

%%% To store ocean area above each cell lower face for each latitude band
%%% Note we add an additional row of zeros along the top because the area
%%% above zero is, obviously, zero
Aocean = zeros(Nlats,Nz+1);

%%% Calculate area of the ocean along each latitude band at each z-level
for k = 1:Nz
  Aocean(:,k+1) = calc_MeridionalTransport(AREAW_z(:,:,k),AREAS_z(:,:,k));
end

%%% Integrate vertically to get total ocean area above cell bottom faces,
%%% for each latitude band
Aocean = cumsum(Aocean,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATE AREA ABOVE ISOPYCNALS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loop through snapshots to calculate areas above isopycnals
Aisop = zeros(Nlats,Nd,Nt);
for n = 1:Nt
% for n = 1
    
  n
  tstart = tic();
  
  %%% Load density field
  DENS = read_nctiles([density_dir 'DENS' num2str(n)],'DENS'); 
  
  %%% Compute areas in density bins
  [AREAW,AREAS] = layers_remap({AREAW_z,AREAS_z},'extensive',DENS,dens_levs,nDblRes);  
    
  %%% Compute vertically integrated areas above each isopycnal
  AREAW = cumsum(AREAW,3,'omitnan');
  AREAS = cumsum(AREAS,3,'omitnan');
  
  %%% Integrate along latitude bands to get the total area above each isopycnal    
  for k = 1:Nd
    Aisop(:,k,n) = calc_MeridionalTransport(AREAW(:,:,k),AREAS(:,:,k));
  end
  
  toc(tstart)
  
end

%%% Calculate time-mean area above isopycnals
Aisop_mean = mean(Aisop,3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATE ISOPYCNAL DEPTHS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% loop through snapshots of areas above isopycnals and interpolate to
%%% find depths
Zisop = zeros(Nlats,Nd+1,Nt);
Zisop_mean = zeros(Nlats,Nd+1);
for j = 1:Nlats
  
  %%% No ocean here so can't interpolate
  if (Aocean(j,2)==0)
    continue;
  end
  
  %%% Truncate Aocean and rf vectors to remove repeated entries at the ends
  %%% of the Aocean vector
  [Aocean_trunc,idx] = unique(Aocean(j,:));
  rf_trunc = rf(idx); 
  
  %%% Interpolate vertically. Note that Aisop(k) is the total area over
  %%% density bins 1 through k, i.e. it corresponds to the density level
  %%% k+1/2.
  for n = 1:Nt 
%   for n = 1
    Zisop(j,2:Nd+1,n) = interp1(Aocean_trunc,rf_trunc,Aisop(j,:,n),'linear','extrap');
  end
  Zisop_mean(j,2:Nd+1) = interp1(Aocean_trunc,rf_trunc,Aisop_mean(j,:),'linear','extrap');
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% POST-PROCESSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Write to a .mat file
lat = [mygrid.LATS_MASKS.lat];
save([products_dir 'Zisop.mat'],'dens_levs','dens_bnds','lat','Aisop','Aocean','Zisop','Aisop_mean','Zisop_mean');

%%% Close Matlab
exit;