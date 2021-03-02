%%%
%%% calcStreamfunction.m
%%%
%%% Computes the isopycnal overturning streamfunction from pre-computed
%%% isopycnal volume fluxes.
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
Nlats = length([mygrid.LATS_MASKS.lat]);
Nz = length([mygrid.RC]);
drf = mk3D(mygrid.DRF,mygrid.hFacC);
dxg = mk3D(mygrid.DXG,mygrid.hFacC);
dyg = mk3D(mygrid.DYG,mygrid.hFacC);



%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% To store computed streamfunctions
PSI = zeros(Nlats,Nd+1,Nt);
PSI_mean = zeros(Nlats,Nd+1);
PSI_zonmean = zeros(Nlats,Nz+1);

%%% Loop through snapshots and compute streamfunction
for n = 1:Nt

  disp(n);
  tstart = tic();
  
  %%% Current date as a datenum (for v4r4 daily data)
  thedate = startdate + n - 1;
  fileidstr = ['_',datestr(thedate,'yyyy_mm_dd')];
  disp(fileidstr);
  
  %%% Load volume fluxes in density bins
%   UFLUX = read_nctiles([ufluxes_dir 'UFLUX' num2str(n)],'UFLUX',1:Nd);       
%   VFLUX = read_nctiles([vfluxes_dir 'VFLUX' num2str(n)],'VFLUX',1:Nd); 
  UFLUX = read_nctiles([ufluxes_dir 'UFLUX' fileidstr],'UFLUX',1:Nd);       
  VFLUX = read_nctiles([vfluxes_dir 'VFLUX' fileidstr],'VFLUX',1:Nd); 
    
  %%% Compute vertically integrated fluxes
  PSIU = cumsum(UFLUX,3,'omitnan');
  PSIV = cumsum(VFLUX,3,'omitnan');
  
  %%% Integrate along latitude bands to get the streamfunction  
  for k = 1:Nd
    PSI(:,k+1,n) = calc_MeridionalTransport(PSIU(:,:,k),PSIV(:,:,k));
  end
  
  toc(tstart)
  
end

%%% Do it all again for the time-mean velocities
UFLUX_mean = read_nctiles([products_dir 'UFLUX_mean'],'UFLUX_mean',1:Nd);       
VFLUX_mean = read_nctiles([products_dir 'VFLUX_mean'],'VFLUX_mean',1:Nd); 
PSIU_mean = cumsum(UFLUX_mean,3,'omitnan');
PSIV_mean = cumsum(VFLUX_mean,3,'omitnan');
for k = 1:Nd
  PSI_mean(:,k+1) = calc_MeridionalTransport(PSIU_mean(:,:,k),PSIV_mean(:,:,k));
end

%%% Compute the zonal-mean overturning streamfunction too
eval(['load ',products_dir,'UVEL',psitype,'_mean.mat;']); 
eval(['load ',products_dir,'VVEL',psitype,'_mean.mat;']); 
eval(['UFLUX_zonmean = UVEL',psitype,'_mean;']);
eval(['VFLUX_zonmean = VVEL',psitype,'_mean;']);
UFLUX_zonmean = UFLUX_zonmean .* mygrid.mskW .* dyg .* drf;
VFLUX_zonmean = VFLUX_zonmean .* mygrid.mskS .* dxg .* drf;
PSIU_zonmean = cumsum(UFLUX_zonmean,3,'omitnan');
PSIV_zonmean = cumsum(VFLUX_zonmean,3,'omitnan');
for k = 1:Nz
  PSI_zonmean(:,k+1) = calc_MeridionalTransport(PSIU_zonmean(:,:,k),PSIV_zonmean(:,:,k));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% POST-PROCESSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Assign to a streamfunction variable according to the selected
%%% streamfunction type
eval(['PSI',psitype,' = PSI;']);
eval(['PSI',psitype,'_mean = PSI_mean;']);
eval(['PSI',psitype,'_zonmean = PSI_zonmean;']);

%%% Write to a .mat file
lat = [mygrid.LATS_MASKS.lat];
save([products_dir 'PSI',psitype,'.mat'],'dens_levs','dens_bnds','lat',['PSI',psitype],['PSI',psitype,'_mean'],['PSI',psitype,'_zonmean']);

%%% Close matlab
exit;
