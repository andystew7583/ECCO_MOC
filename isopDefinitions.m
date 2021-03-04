%%% 
%%% isopDefinitions.m 
%%%
%%% Defines global variables required for all calculations.
%%%

%%% Directories
basedir = ['Version4' filesep 'Release4'];
% basedir = ['Version4' filesep 'Release3'];
% basedir = ['Version5' filesep 'Alpha'];
products_dir = [basedir filesep 'myproducts_daily' filesep]; %%% Directory in which to store output of calculations
% products_dir = [basedir filesep 'myproducts_monthly' filesep]; %%% Directory in which to store output of calculations
density_dir = [products_dir 'DENS' filesep]; %%% Name of directory in which to store computed density variable
ufluxes_dir = [products_dir 'UFLUX' filesep]; %%% Name of directory in which to store computed u-fluxes in density space
vfluxes_dir = [products_dir 'VFLUX' filesep]; %%% Name of directory in which to store computed u-fluxes in density space
psi_dir = [products_dir 'PSI' filesep]; %%% Name of directory in which to store computed streamfunction in density space
ifs_dir = [products_dir 'IFS' filesep]; %%% Name of directory in which to store computed form stresses in density space
ECCO_grid_dir = [basedir filesep 'nctiles_grid']; %%% Directory holding ECCO data
ECCO_data_dir = [basedir filesep 'nctiles_daily']; %%% Directory holding ECCO data
% ECCO_data_dir = [basedir filesep 'nctiles_monthly']; %%% Directory holding ECCO data

%%% Number of snapshots
%%% N.B. 312 snapshots are available in Version5, but 288 may be preferable
%%% for comparison with Version4
startdate = datenum('1992-01-01');
enddate = datenum('2016-12-31');
% Nt = 288;
Nt = enddate-startdate+1;

%%% Time in months
tt = 1:Nt;

%%% Reference pressure for potential density calculation
P_ref = 2000;

%%% Reference density
rhoConst = 1029;

%%% Gravitational acceleration
gravity = 9.81;

%%% Density levels to use for MOC calculation
% dens_levs = 1000 + [24:1:27 27.5:.5:34 34.25:.25:35 35.1:.1:36 36.05:0.05:36.9 36.92:0.02:37.2 37.3:.1:37.8 38.3 38.8];
dens_levs = 1000 + [24:.5:27 27.25:.25:34 34.125:.125:35 35.05:.05:36 36.025:0.025:36.9 36.91:0.01:37.2 37.25:.05:37.8 38.05:.25:38.8];

%%% Boundaries of integrals in density space (see layers_remap.m)
dens_bnds = [dens_levs(1)-(dens_levs(2)-dens_levs(1))/2 ...
             (dens_levs(1:end-1)+dens_levs(2:end))/2 ... 
             dens_levs(end)+(dens_levs(end)-dens_levs(end-1))/2];

%%% Number of density levels
Nd = length(dens_levs);

%%% Number of vertical resolution doublings to use in flux calculation
nDblRes = 3;

%%% Select streamfunction to calculate
%%% 'eul' -> Eulerian-mean
%%% 'bol' -> Eddy bolus
%%% 'tot' -> Total, i.e. sum of E-M and bolus
psitype = 'tot';