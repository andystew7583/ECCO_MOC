%%%
%%% testUVELMASS.m
%%%
%%% Tests volume fluxes reported by different ECCO diagnostics.



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
grid_load;

%%% Make mygrid accessible in current workspace:
gcmfaces_global;

%%% Directory from which to read ECCO output data
myenv.nctilesdir = fullfile(ECCO_data_dir,filesep);



%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Select time index
n = 1;
   
%%% Load Eulerian and bolus horizontal velocities
listVars={'UVEL','VVEL','UVELSTAR','VVELSTAR','UVELMASS','VVELMASS','GM_PsiX','GM_PsiY'};
for vvv=1:length(listVars)    
  vv=listVars{vvv};
  tmp1=read_nctiles([myenv.nctilesdir vv '/' vv],vv,n);
  tmp1=mean(tmp1,4);
  tmp1(mygrid.mskC==0)=NaN;
  eval([vv '=tmp1;']);
end

%%% Calculate bolus velocity
[UVELbol,VVELbol,fldWbolus]=calc_bolus(GM_PsiX,GM_PsiY);
UVELbol=UVELbol.*mygrid.mskW; VVELbol=VVELbol.*mygrid.mskS;

%%% Check the difference
max(max(max(abs(UVEL.*mygrid.hFacW-UVELMASS))))
max(max(max(abs(VVEL.*mygrid.hFacS-VVELMASS))))
max(max(max(abs(UVELSTAR-UVELbol))))
max(max(max(abs(VVELSTAR-VVELbol))))

