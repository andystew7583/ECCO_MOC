%%%
%%% calcMeanIsopFluxes.m
%%%
%%% Calculates mean volume fluxes in mean density bins in ECCO output.
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




%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Loop through snapshots and compute time-mean velocity fields
UVELeul_mean = zeros(mygrid.hFacC);
VVELeul_mean = zeros(mygrid.hFacC);
UVELbol_mean = zeros(mygrid.hFacC);
VVELbol_mean = zeros(mygrid.hFacC);
for n = 1:Nt
  
  n
  tstart = tic();
  
  %%% Load Eulerian and bolus horizontal velocities
%   listVars={'UVEL','VVEL','UVELSTAR','VVELSTAR'};
  listVars={'UVELMASS','VVELMASS','GM_PsiX','GM_PsiY'};
  for vvv=1:length(listVars)    
    vv=listVars{vvv};
    tmp1=read_nctiles([myenv.nctilesdir vv '/' vv],vv,n);
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
  
  %%% Add to time-averages  
  UVELeul_mean = UVELeul_mean + UVELMASS;
  VVELeul_mean = VVELeul_mean + VVELMASS;
  UVELbol_mean = UVELbol_mean + UVELbol;
  VVELbol_mean = VVELbol_mean + VVELbol;
  
  toc(tstart)  

end

%%% Calculate time means
UVELeul_mean = UVELeul_mean / Nt;
VVELeul_mean = VVELeul_mean / Nt;
UVELbol_mean = UVELbol_mean / Nt;
VVELbol_mean = VVELbol_mean / Nt;

%%% Load time-mean density field
load([products_dir 'DENS_mean.mat'],'DENS_mean');

%%% Do it all again for the time-mean velocities/density total transport velocity
UVELtot_mean = UVELeul_mean + UVELbol_mean;
VVELtot_mean = VVELeul_mean + VVELbol_mean;
eval(['UVELcalc_mean = UVEL',psitype,'_mean;']);
eval(['VVELcalc_mean = VVEL',psitype,'_mean;']);
UVELcalc_mean = UVELcalc_mean .* mygrid.mskW .* dyg .* drf;
VVELcalc_mean = VVELcalc_mean .* mygrid.mskS .* dxg .* drf;
[UFLUX_mean,VFLUX_mean] = layers_remap({UVELcalc_mean,VVELcalc_mean},'extensive',DENS_mean,dens_levs,nDblRes);  
write2nctiles([products_dir 'UFLUX_mean'],UFLUX_mean,1, ...
  {'fldName','UFLUX_mean'}, ...
  {'longName',['Mean u-volume fluxes in mean density bins']}, ...
  {'units','m^3/s'}, ...
  {'coord','lon lat dep tim'} ...
);
write2nctiles([products_dir 'VFLUX_mean'],VFLUX_mean,1, ...
  {'fldName','VFLUX_mean'}, ...
  {'longName',['Mean v-volume fluxes in mean density bins']}, ...
  {'units','m^3/s'}, ...
  {'coord','lon lat dep tim'} ...
);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% POST-PROCESSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Save time-mean quantities to .mat files
save([products_dir 'UVELeul_mean.mat'],'UVELeul_mean');
save([products_dir 'VVELeul_mean.mat'],'VVELeul_mean');
save([products_dir 'UVELbol_mean.mat'],'UVELbol_mean');
save([products_dir 'VVELbol_mean.mat'],'VVELbol_mean');
save([products_dir 'UVELtot_mean.mat'],'UVELtot_mean');
save([products_dir 'VVELtot_mean.mat'],'VVELtot_mean');




