%%%
%%% calcDensity.m
%%%
%%% Calculates density at each ECCO output snapshot
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

%%% Make mygrid accessible in current workspace:
gcmfaces_global;

%%% Directory from which to read ECCO outputdata
myenv.nctilesdir = fullfile(ECCO_data_dir,filesep);

%%% Output directory
mkdir(density_dir);



%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Vector of reference pressures
pp = -rhoConst*gravity*mygrid.RC;

%%% Create gmfaces object of pressures at cell centers
PP = mk3D(pp,mygrid.hFacC);

%%% To calculate time-mean density while we're at it
DENS_mean = zeros(mygrid.hFacC);

%%% Loop through snapshots and compute density
for n = 1:Nt
  
  disp(n);
  
  if (isV4R4)
    %%% Current date as a datenum (for v4r4 daily data)
    thedate = startdate + n - 1;
    fileidstr = ['_',datestr(thedate,'yyyy_mm_dd')];
    disp(fileidstr);
  else
    fileidstr = num2str(n,'%.4d');
  end
    
  %%% Load potential temperature and salinity
  listVars={'THETA','SALT'};
  for vvv=1:length(listVars)
    vv=listVars{vvv};
    if (isV4R4)
      tmp1 = read_nctiles_daily([myenv.nctilesdir vv filesep vv fileidstr '.nc'],vv);
    else
      tmp1 = read_nctiles([myenv.nctilesdir vv filesep vv],vv,n);        
    end    
    eval([vv '=tmp1;']);
  end
  
  %%% Calculate density
  [RHOP,RHOIS,RHOR] = density(THETA,SALT,PP,P_ref);
  DENS = RHOR;
  
  %%% Add to time-averages
  DENS_mean = DENS_mean + DENS;
    
  %%% Save to NetCDF files in gcmfaces format    
  write2nctiles([density_dir 'DENS' fileidstr],DENS,1, ...
    {'fldName','DENS'}, ...
    {'longName',['Potential density referenced to ',num2str(P_ref),' db']}, ...
    {'units','kg/m^3'}, ...
    {'coord','lon lat dep tim'} ...
  );  
  
end

%%% Calculate time means
DENS_mean = DENS_mean / Nt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% POST-PROCESSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Save time-mean quantities to .mat files
save([products_dir 'DENS_mean.mat'],'DENS_mean');