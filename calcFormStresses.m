%%%
%%% calcFormStresses.m
%%% 
%%% Computes isopycnal and topographic form stresses from ECCO output.
%%% 
function calcFormStresses (n)


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




  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CALCCULATION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%

  tstart = tic();

  disp(n);

  %%% Current date as a datenum (for v4r4 daily data)
  thedate = startdate + n - 1;
  fileidstr = ['_',datestr(thedate,'yyyy_mm_dd')];
  disp(fileidstr);

  %%% Load density and pressure
  DENS = read_nctiles([density_dir 'DENS' fileidstr],'DENS');  
  PHIHYDcR = read_nctiles_daily([myenv.nctilesdir 'PHIHYDcR' filesep 'PHIHYDcR' fileidstr '.nc'],'PHIHYDcR');

  %%% Convert to arrays
  DENS = convert2array(DENS);
  PHIHYDcR = convert2array(PHIHYDcR);
  XC = convert2array(mygrid.XC);
  YC = convert2array(mygrid.YC);
  RC = mygrid.RC;
  DRF = mygrid.DRF;
  hFacC = convert2array(mygrid.hFacC);
  hFacW = convert2array(mygrid.hFacW);
  hFacS = convert2array(mygrid.hFacS);
  AngleCS = convert2array(mygrid.AngleCS);
  DXC = convert2array(mygrid.DXC);
  DYC = convert2array(mygrid.DYC);

  %%% Modify matrices based on rotation of LLC grid. Note that this only
  %%% works for latitudes between 70S and 57N.
  AngleCS3D = repmat(AngleCS,[1 1 length(RC)]);
  hFacW(AngleCS3D<0.5) = hFacS(AngleCS3D<0.5);
  DXC(AngleCS<0.5) = DYC(AngleCS<0.5);
  hFacE = hFacW([2:end 1],:,:);

  %%% Compute TFS
  msk = ones(size(PHIHYDcR));
  msk(hFacC==0) = NaN;
  TFS_z = rhoConst*squeeze(nansum(PHIHYDcR.*msk.*(hFacE-hFacW),1));
  TFS = sum(TFS_z.*repmat(reshape(DRF,[1 length(DRF)]),[size(TFS_z,1) 1]),2);
  secLen = sum(DXC,1)'; %%% Zonal section length
  secLats = mean(YC,1);
  
  %%% Compute IFS
  IFS = zeros(size(YC,2),length(dens_bnds));
  for m = 1:length(dens_bnds)

    %%% Mask for pressure points, excluding regions where density is higher
    %%% than the current bin
    msk = ones(size(PHIHYDcR));
    msk(hFacC==0) = NaN;
    msk(DENS>dens_bnds(m)) = NaN;

    %%% Modified hFacs based on "topography" of this density surface
    mskW = hFacW;
    mskW(DENS>dens_bnds(m)) = 0;
    mskW(DENS([end 1:end-1],:,:)>dens_bnds(m)) = 0;
    mskE = mskW([2:end 1],:,:);

    %%% Integrate IFS zonally and vertically
    IFS_z = rhoConst*squeeze(nansum(PHIHYDcR.*msk.*(mskE-mskW),1));
    IFS(:,m) = sum(IFS_z.*repmat(reshape(DRF,[1 length(DRF)]),[size(IFS_z,1) 1]),2);

  end
  
  %%% Write to output file
  save([ifs_dir 'IFS' fileidstr],'IFS','TFS','secLen','dens_bnds','dens_levs','secLats');

  toc(tstart)
  
end