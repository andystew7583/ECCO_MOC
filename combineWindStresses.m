%%% 
%%% combineWindStresses.m
%%%
%%% Loads all wind stress data, rotates, and saves as a single .mat file.
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

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




%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

for n=1:Nt

  tstart = tic();

  disp(n);

  %%% Current date as a datenum (for v4r4 daily data)
  thedate = startdate + n - 1;
  fileidstr = ['_',datestr(thedate,'yyyy_mm_dd')];
  disp(fileidstr);

  %%% Load density and pressure
  oceTAUX = read_nctiles_daily([myenv.nctilesdir 'oceTAUX' filesep 'oceTAUX' fileidstr '.nc'],'oceTAUX');
  oceTAUY = read_nctiles_daily([myenv.nctilesdir 'oceTAUY' filesep 'oceTAUY' fileidstr '.nc'],'oceTAUY');

  %%% Convert to arrays
  oceTAUX = convert2array(oceTAUX);
  oceTAUY = convert2array(oceTAUY);
  XC = convert2array(mygrid.XC);
  YC = convert2array(mygrid.YC);
  AngleCS = convert2array(mygrid.AngleCS);
  DXC = convert2array(mygrid.DXC);
  DYC = convert2array(mygrid.DYC);

  %%% Modify matrices based on rotation of LLC grid. Note that this only
  %%% works for latitudes between 70S and 57N.
  DXC(AngleCS<0.5) = DYC(AngleCS<0.5);
  oceTAUX(AngleCS<0.5) = oceTAUY(AngleCS<0.5);
  
  %%% Initialize wind stress storage matrix
  if (n == 1)
    WS = zeros(size(YC,2),Nt);
  end

  %%% Compute TFS
  WS(:,n) = squeeze(sum(oceTAUX.*DXC,1));
  secLen = sum(DXC,1)'; %%% Zonal section length
  secLats = mean(YC,1);
  
  toc(tstart)

end

%%% Write to output file
save([products_dir 'WS.mat'],'WS','secLen','secLats');

