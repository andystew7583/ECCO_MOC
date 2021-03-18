%%%
%%% combineFormStresses.m
%%%
%%% Combines individual IFS .mat files into a single .mat file.
%%%

clear all;

%%% Load static parameters
isopDefinitions;

%%% Loop over snapshots
for n = 1:Nt
  
  disp(n);
  tstart = tic();

  %%% Current date as a datenum (for v4r4 daily data)
  thedate = startdate + n - 1;
  fileidstr = ['_',datestr(thedate,'yyyy_mm_dd')];
  disp(fileidstr);

  %%% Load form stresses from .mat file
  IFS_tmp = load([ifs_dir 'IFS',fileidstr,'.mat'],'IFS');
  IFS_tmp = IFS_tmp.IFS;
  TFS_tmp = load([ifs_dir 'IFS',fileidstr,'.mat'],'TFS');
  TFS_tmp = TFS_tmp.TFS;
  
  %%% Initialize storage matrix
  if (n == 1)
    load([ifs_dir 'IFS',fileidstr,'.mat']); %%% Load full .mat file to get coordinate vectors too
    IFS = zeros(size(IFS_tmp,1),size(IFS_tmp,2),Nt);
    TFS = zeros(length(TFS_tmp),Nt);
  end  
  
  %%% Add streamfunction snapthot to storage matrix
  IFS(:,:,n) = IFS_tmp;
  TFS(:,n) = TFS_tmp;
  
  toc(tstart)
  
end

%%% Write to .mat file
save([products_dir,'IFS','.mat'],'dens_levs','dens_bnds','secLats','secLen','IFS','TFS','-v7.3');