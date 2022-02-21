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

%%% Modify latitude masks based on selected region type
switch (regionType)

  case 'AtlOnly'
    [atlMskC,atlMskW,atlMskS] = v4_basin('atlExt');
    for j=1:length(mygrid.LATS_MASKS)
      mygrid.LATS_MASKS(j).mskCint = mygrid.LATS_MASKS(j).mskCint .* atlMskC;
      mygrid.LATS_MASKS(j).mskCedge = mygrid.LATS_MASKS(j).mskCedge .* atlMskC;
      mygrid.LATS_MASKS(j).mskWedge = mygrid.LATS_MASKS(j).mskWedge .* atlMskW;
      mygrid.LATS_MASKS(j).mskSedge = mygrid.LATS_MASKS(j).mskSedge .* atlMskS;
    end

  case 'PacOnly'
    [pacMskC,pacMskW,pacMskS] = v4_basin('pacExt');
    [indMskC,indMskW,indMskS] = v4_basin('indExt');
    for j=1:length(mygrid.LATS_MASKS)
      mygrid.LATS_MASKS(j).mskCint = mygrid.LATS_MASKS(j).mskCint .* (pacMskC+indMskC);
      mygrid.LATS_MASKS(j).mskCedge = mygrid.LATS_MASKS(j).mskCedge .* (pacMskC+indMskC);
      mygrid.LATS_MASKS(j).mskWedge = mygrid.LATS_MASKS(j).mskWedge .* (pacMskW+indMskW);
      mygrid.LATS_MASKS(j).mskSedge = mygrid.LATS_MASKS(j).mskSedge .* (pacMskS+indMskS);
    end
    
  case 'AtlPac'
    %%% Do nothing
    
  otherwise
    error('Unrecognized region type specified in isopDefinitions');
    
end
  



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
  
  %%% File naming system depends on ECCO version
  if (isV4R4)
    thedate = startdate + n - 1;
    fileidstr = ['_',datestr(thedate,'yyyy_mm_dd')];        
  else
    fileidstr = num2str(n,'%.4d');    
  end
  
  %%% Load volume fluxes in density bins
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
save([products_dir 'PSI',psitype,'_',regionType,'.mat'],'dens_levs','dens_bnds','lat',['PSI',psitype],['PSI',psitype,'_mean'],['PSI',psitype,'_zonmean']);

%%% Close matlab
exit;

