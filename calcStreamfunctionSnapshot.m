%%%
%%% calcStreamfunctionSnapshot.m
%%%
%%% Computes the isopycnal overturning streamfunction from pre-computed
%%% isopycnal volume fluxes.
%%%
function calcStreamfunctionSnapshot (n)


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

  %%% Tools for calculating properties in latitude bands
  if ~isfield(mygrid,'LATS_MASKS')  
    gcmfaces_lines_zonal;   
  end
  Nlats = length([mygrid.LATS_MASKS.lat]);




  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CALCULATION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%

  %%% To store computed streamfunctions
  PSI = zeros(Nlats,Nd+1);

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

  %%% Write to a .mat file
  lat = [mygrid.LATS_MASKS.lat];
  save([psi_dir 'PSI',psitype,fileidstr,'.mat'],'dens_levs','dens_bnds','lat','PSI');

end