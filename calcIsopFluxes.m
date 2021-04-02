%%%
%%% calcIsopFluxes.m
%%%
%%% Calculates volume fluxes in density bins at each ECCO output snapshot.
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

%%% Make mygrid accessible in current workspace:
gcmfaces_global;

%%% Directory from which to read ECCO output data
myenv.nctilesdir = fullfile(ECCO_data_dir,filesep);
drf = mk3D(mygrid.DRF,mygrid.hFacC);
dxg = mk3D(mygrid.DXG,mygrid.hFacC);
dyg = mk3D(mygrid.DYG,mygrid.hFacC);

%%% Output directories
mkdir(ufluxes_dir);
mkdir(vfluxes_dir);



%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Loop through snapshots and compute fluxes in density classes. 
for n = 1:Nt
    
  tstart = tic();  
  
  %%% Load density field
  if (isV4R4)
    thedate = startdate + n - 1;
    fileidstr = ['_',datestr(thedate,'yyyy_mm_dd')];
    disp(fileidstr);
    DENS = read_nctiles([density_dir 'DENS' fileidstr],'DENS');       
  else
    disp(n);  
    DENS = read_nctiles([density_dir 'DENS' num2str(n)],'DENS');       
  end
    
  %%% Load Eulerian and bolus horizontal velocities
%   listVars={'UVEL','VVEL','UVELSTAR','VVELSTAR'};
  listVars={'UVELMASS','VVELMASS','GM_PsiX','GM_PsiY'};
  for vvv=1:length(listVars)    
    vv=listVars{vvv};
%     tmp1=read_nctiles([myenv.nctilesdir vv '/' vv],vv,n);
    tmp1 = read_nctiles_daily([myenv.nctilesdir vv filesep vv fileidstr '.nc'],vv);
    tmp1=mean(tmp1,4);
    tmp1(mygrid.mskC==0)=NaN;
    eval([vv '=tmp1;']);
  end
  
  %%% Eulerian-mean velocity, weighted by grid cells hFacs
  UVELeul = UVELMASS;
  VVELeul = VVELMASS;
  
  %%% Calculate bolus velocity, weighted by grid cell hFacs
  [UVELbol,VVELbol,fldWbolus]=calc_bolus(GM_PsiX,GM_PsiY);  
  UVELbol = UVELbol.* mygrid.hFacW;
  VVELbol = VVELbol.* mygrid.hFacS;
  
  %%% Calculate total transport velocity
  UVELtot = UVELeul + UVELbol;
  VVELtot = VVELeul + VVELbol;
  
  %%% Select the velocity field on which to perform calculations
  eval(['UVELcalc = UVEL',psitype,';']);
  eval(['VVELcalc = VVEL',psitype,';']);

  %%% Apply NaN-masks and weight velocities by fractional cell thicknesses    
  UVELcalc = UVELcalc .* mygrid.mskW .* dyg .* drf;
  VVELcalc = VVELcalc .* mygrid.mskS .* dxg .* drf;
  
  %%% Compute fluxes in density bins
  [UFLUX,VFLUX] = layers_remap({UVELcalc,VVELcalc},'extensive',DENS,dens_levs,nDblRes);  
    
  %%% Save to NetCDF files in gcmfaces format  
%   write2nctiles([ufluxes_dir 'UFLUX' num2str(n)],UFLUX,1, ...
  write2nctiles([ufluxes_dir 'UFLUX' fileidstr],UFLUX,1, ...
    {'fldName','UFLUX'}, ...
    {'longName',['u-volume fluxes in density bins at time level ' num2str(n)]}, ...
    {'units','m^3/s'}, ...
    {'coord','lon lat dep tim'} ...
  );  
%   write2nctiles([vfluxes_dir 'VFLUX' num2str(n)],VFLUX,1, ...
  write2nctiles([vfluxes_dir 'VFLUX' fileidstr],VFLUX,1, ...
    {'fldName','VFLUX'}, ...
    {'longName',['v-volume fluxes in density bins at time level ' num2str(n)]}, ...
    {'units','m^3/s'}, ...
    {'coord','lon lat dep tim'} ...
  );
  
  toc(tstart)  

end
