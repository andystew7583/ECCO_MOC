%%% Start by clearing memory
clear all;

%add gcmfaces and MITprof directories to Matlab path:
p = genpath('gcmfaces/'); addpath(p);
p = genpath('gcmfaces/'); addpath(p);
p = genpath('m_map/'); addpath(p);

%load all grid variables from nctiles_grid/ into mygrid:
grid_load;

%make mygrid accessible in current workspace:
gcmfaces_global;

%%% Directory from which to read ECCO outputdata
myenv.nctilesdir=fullfile('nctiles_monthly',filesep);


%%%%%%%%%%%%%%%%%
%do computations:
%%%%%%%%%%%%%%%%%

listDiags={};

% [listTimes]=diags_list_times;
diags.listTimes=1;

%part 1:

% listVars={'UVEL','VVEL','UVELSTAR','UVELSTAR'};
listVars={'UVEL'};

for vvv=1:length(listVars);
  vv=listVars{vvv};
  tmp1=read_nctiles([myenv.nctilesdir vv '/' vv],vv,1);    
  eval([vv '=tmp1;']);
end

%apply NaN-masks:
UVELMASS=UVELMASS.*mygrid.mskW;
VVELMASS=VVELMASS.*mygrid.mskS;

if myenv.verbose>0; gcmfaces_msg('* call calc_bolus :  compute bolus component from GM streamfunction');end;

%add bolus component from GM scheme:
[UVELbol,VVELbol,fldWbolus]=calc_bolus(GM_PsiX,GM_PsiY);
UVELbol=UVELbol.*mygrid.mskW; VVELbol=VVELbol.*mygrid.mskS;
UVELtot=UVELMASS+UVELbol; VVELtot=VVELMASS+VVELbol;

%compute transports:
listDiags={listDiags{:},'fldBAR','gloOV','fldTRANSPORTS','gloMT_FW'};
if myenv.verbose>0; gcmfaces_msg('* call calc_barostream : comp. barotropic stream function');end;
[fldBAR]=calc_barostream(UVELMASS,VVELMASS);
if myenv.verbose>0; gcmfaces_msg('* call calc_overturn : comp. residual overturning stream function');end;
[gloOV]=calc_overturn(UVELtot,VVELtot);
if myenv.verbose>0; gcmfaces_msg('* call calc_transports : comp. transects transports');end;
[fldTRANSPORTS]=1e-6*calc_transports(UVELMASS,VVELMASS,mygrid.LINES_MASKS,{'dh','dz'});
if myenv.verbose>0; gcmfaces_msg('* call calc_MeridionalTransport : comp. meridional seawater transport');end;
[gloMT_FW]=1e-6*calc_MeridionalTransport(UVELMASS,VVELMASS,1);

end;%if ~isempty(missingVars);

%part 2:

listVars={'THETA','SALT','ADVx_TH','ADVy_TH','ADVx_SLT','ADVy_SLT'};
listVars={listVars{:},'DFxE_TH','DFyE_TH','DFxE_SLT','DFyE_SLT'};

missingVars={};
for vv=1:length(listVars);
  tmp1=[myenv.nctilesdir listVars{vv} filesep listVars{vv} '*nc'];
  if isempty(dir(tmp1)); missingVars={missingVars{:},tmp1}; end;
end;

if ~isempty(missingVars);
    fprintf('\n example_transports could not find the following files ---> skipping related computation!\n');
    disp(missingVars');
else;

if myenv.verbose>0; gcmfaces_msg('* load tracer and transports fields');end;
listDiags={listDiags{:},'fldTzonmean','fldSzonmean','gloMT_H','gloMT_SLT'};
for vvv=1:length(listVars);
    vv=listVars{vvv};
    tmp1=read_nctiles([myenv.nctilesdir vv '/' vv],vv);
    tmp1=mean(tmp1,4);
    tmp1(mygrid.mskC==0)=NaN;
    eval([vv '=tmp1;']);
end;

if myenv.verbose>0; gcmfaces_msg('* call calc_zonmean_T : comp. zonal mean temperature');end;
[fldTzonmean]=calc_zonmean_T(THETA);
if myenv.verbose>0; gcmfaces_msg('* call calc_zonmean_T : comp. zonal mean salinity');end;
[fldSzonmean]=calc_zonmean_T(SALT);

if myenv.verbose>0; gcmfaces_msg('* call calc_MeridionalTransport : comp. meridional heat transport');end;
tmpU=(ADVx_TH+DFxE_TH); tmpV=(ADVy_TH+DFyE_TH);
[gloMT_H]=1e-15*4e6*calc_MeridionalTransport(tmpU,tmpV,0);
if myenv.verbose>0; gcmfaces_msg('* call calc_MeridionalTransport : comp. meridional salt transport');end;
tmpU=(ADVx_SLT+DFxE_SLT); tmpV=(ADVy_SLT+DFyE_SLT);
[gloMT_SLT]=1e-6*calc_MeridionalTransport(tmpU,tmpV,0);

end;%if ~isempty(missingVars);

%part 3: format output

for ddd=1:length(listDiags);
    dd=listDiags{ddd};
    eval(['diags.' dd '=' dd ';']);
end;

if myenv.verbose>0;
    gcmfaces_msg('*** leaving example_transports');
    gcmfaces_msg('===============================================','');
end;

