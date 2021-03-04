%%%
%%% combineStreamfunction.m
%%%
%%% Combines individual streamfunction .mat files into a single .mat file.
%%%

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

  %%% Load streamfunction from .mat file
  PSI_tmp = load([psi_dir 'PSI',psitype,fileidstr,'.mat'],'PSI');
  PSI_tmp = PSI_tmp.PSI;
  
  %%% Initialize storage matrix
  if (n == 1)
    load([psi_dir 'PSI',psitype,fileidstr,'.mat']); %%% Load full .mat file to get coordinate vectors too
    PSI = zeros(size(PSI_tmp,1),size(PSI_tmp,2),Nt);
  end  
  
  %%% Add streamfunction snapthot to storage matrix
  PSI(:,:,n) = PSI_tmp(:,:,n);
  
  toc(tstart)
  
end

%%% Write to .mat file
save([products_dir,'PSI',psitype,'.mat'],'dens_levs','dens_bnds','lat','PSI')