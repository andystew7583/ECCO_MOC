%%%
%%% generateIsopFluxBatch.m
%%%
%%% Generates a batch of cluster jobs to compute isopycnal fluxes from ECCO
%%% output.
%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Start by clearing memory
clear all;

%%% Load global variables
isopDefinitions;

%%% Output directories
mkdir(ufluxes_dir);
mkdir(vfluxes_dir);

%%% To store script files
script_dir = 'IsopFluxScripts';
mkdir(fullfile('.',script_dir));

%%% Open batch file
fid = fopen(fullfile('.',script_dir,'run_calcIsopFluxes.sh'),'w');

%%% Loop through snapshots and compute fluxes in density classes. 
for n = 1:Nt
  
  %%% Create SGE submission file
  sge_fname = createSGEfile('./IsopFluxScripts',2,'all.q','astewart@atmos.ucla.edu','calcIsopFluxSnapshot',n);
  
  %%% Add execution line to batch file
  fprintf(fid,'qsub %s\n',sge_fname);
  
end

%%% Close batch file
fclose(fid);