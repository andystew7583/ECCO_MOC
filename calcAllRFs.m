%%%
%%% calcAllRFs.m
%%%
%%% Creates some plots of the IFS and overturning streamfunction derived
%%% from the ECCOv4r4 daily output.
%%%

%%% Required paths
addpath colormaps;
addpath CDT/cdt;

%%% Static definitions
isopDefinitions;

%%% Load precomputed streamfunction and IFS time series
load(fullfile(products_dir,'PSIbol.mat'))
PSIbol = PSI;
load(fullfile(products_dir,'PSIeul.mat'))
PSIeul = PSI;
load(fullfile(products_dir,'PSItot.mat'))
load(fullfile(products_dir,'IFS.mat'));
load(fullfile(products_dir,'WS.mat'));

%%% Coriolis parameter matrix on IFS points
Omega = 366/365*2*pi/86400;
ff = 2*Omega*sind(secLats);
FF = repmat(ff',[1 length(dens_bnds)]);

%%% Estimate eddy form stress
EFS = zeros(length(secLats),Nd+1);
for m=1:Nd+1
  EFS(:,m) = interp1(lat,mean(PSIbol(:,m),3),secLats,'linear');
end
EFS = -EFS.*rhoConst.*FF;

%%% Latitudinally-averaged time series metrics of overturning and IFS
ymin = -60;
ymax = -50;
dens_psimax = 1037.1;
didx_psimax = find(dens_bnds>dens_psimax,1,'first');
yidx_psi = find((lat>ymin) & (lat<ymax));
yidx_ifs = find((secLats>ymin) & (secLats<ymax));
PSImean = squeeze(mean(PSI(yidx_psi,didx_psimax,:),1));
PSIeul_mean = squeeze(mean(PSIeul(yidx_psi,didx_psimax,:),1));
PSIbol_mean = squeeze(mean(PSIbol(yidx_psi,didx_psimax,:),1));
IFSmean = squeeze(mean(IFS(yidx_ifs,didx_psimax,:),1));
EFSmean = -PSIbol_mean*rhoConst*mean(ff(yidx_ifs));
TFSmean = squeeze(mean(TFS(yidx_ifs,:),1))';
WSmean = squeeze(mean(WS(yidx_ifs,:)))';
tt = startdate + (0:Nt-1);

%%% Define window parameters for RFs
Lwin = 21*365;
Dwin = 30;
Nrf = 10*365;
Nwin = floor((Nt-Lwin)/Dwin) + 1;
tt_rf = 0:Nrf-1;

%%% Compute batches of RFs with shifted windows
[PSI_WS_RF,PSI_WS_eps] = calcRFbatch(WSmean,PSImean*rhoConst*mean(ff(yidx_ifs)),Nrf,Lwin,Dwin,Nwin);
[TFS_WS_RF,TFS_WS_eps] = calcRFbatch(WSmean,TFSmean,Nrf,Lwin,Dwin,Nwin);
[IFS_WS_RF,IFS_WS_eps] = calcRFbatch(WSmean,IFSmean,Nrf,Lwin,Dwin,Nwin);
[EFS_WS_RF,EFS_WS_eps] = calcRFbatch(WSmean,EFSmean,Nrf,Lwin,Dwin,Nwin);

%%% Store in a .mat file
save('RFs.mat','Lwin','Dwin','Nrf','Nwin','tt_rf',...
  'PSI_WS_RF','PSI_WS_eps',...
  'TFS_WS_RF','TFS_WS_eps',...
  'IFS_WS_RF','IFS_WS_eps',...
  'EFS_WS_RF','EFS_WS_eps');



%%%
%%% calcRFbatch
%%%
%%% Helper function to compute RFs of Y with respect to X,
%%% where each Rf is Nrf points long, and the full time series is split
%%% into Nwin windows of length Lwin, separated by Dwin points.
%%%
function [rf,eps] = calcRFbatch(X,Y,Nrf,Lwin,Dwin,Nwin)

  rf = zeros(Nrf,Nwin);
  eps = zeros(Lwin-Nrf+1,Nwin);
  for n = 1:Nwin
    n    
    idx = (n-1)*Dwin + (1:Lwin);
    [rf_tmp,eps_tmp] = calcRF(X(idx),Y(idx),Nrf);
    rf(:,n) = rf_tmp;
    eps(:,n) = eps_tmp;    
  end

end