%%%
%%% generateIFSBatch.m
%%%
%%% Generates a batch of cluster jobs to compute form stresses from ECCO output.
%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Start by clearing memory
clear all;

%%% Load global variables
isopDefinitions;

%%% Output directories
mkdir(ifs_dir);

%%% To store script files
script_dir = 'IFSScripts';
mkdir(fullfile('.',script_dir));

%%% Open batch file
fid = fopen(fullfile('.',script_dir,'run_calcIFS.sh'),'w');

%%% Loop through snapshots and compute IFSs
for n = 1:Nt
  
  %%% Create SGE submission file
  sge_fname = createSGEfile(fullfile('.',script_dir),2,'all.q','astewart@atmos.ucla.edu','calcFormStresses',n);
  
  %%% Add execution line to batch file
  fprintf(fid,'qsub %s\n',sge_fname);
  
end

%%% Close batch file
fclose(fid);
