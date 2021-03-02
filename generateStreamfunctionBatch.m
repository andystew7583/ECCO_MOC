%%%
%%% generateIsopFluxBatch.m
%%%
%%% Generates a batch of cluster jobs to compute streamfunctions from isopycnal 
%%% fluxes, derived from ECCO output.
%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Start by clearing memory
clear all;

%%% Load global variables
isopDefinitions;

%%% Output directories
mkdir(psi_dir);

%%% To store script files
script_dir = 'StreamfunctionScripts';
mkdir(fullfile('.',script_dir));

%%% Open batch file
fid = fopen(fullfile('.',script_dir,'run_calcStreamfunction.sh'),'w');

%%% Loop through snapshots and compute fluxes in density classes. 
for n = 81:Nt
  
  %%% Create SGE submission file
  sge_fname = createSGEfile(fullfile('.',script_dir),2,'all.q','astewart@atmos.ucla.edu','calcStreamfunctionSnapshot',n);
  
  %%% Add execution line to batch file
  fprintf(fid,'qsub %s\n',sge_fname);
  
end

%%% Close batch file
fclose(fid);
