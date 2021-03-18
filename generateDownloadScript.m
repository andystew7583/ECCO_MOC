%%%
%%% generateDownloadScript.m
%%%
%%% Generates a bash script to download fields from the ECCOv4r4
%%% repository.
%%%

isopDefinitions;

%%% Options
varname = 'oceTAUY';
thedate = startdate;

%%% Open script file
fid = fopen(['Version4/Release4/nctiles_daily/',varname,'/download_',varname,'.sh'],'w');

%%% Loop through all days in specified range
while (thedate <= enddate)

  %%% Current year and day of the year
  yearnum = str2num(datestr(thedate,'yyyy'));
  yearday = thedate - datenum([num2str(yearnum),'-01-01'])+1;
  
  %%% Construct target URL
  url = ['https://data.nas.nasa.gov/ecco/download_data.php?file=/eccodata/llc_90/ECCOv4/Release4/nctiles_daily/'];
  url = [url,varname,'/'];
  url = [url,datestr(thedate,'yyyy'),'/'];
  url = [url,num2str(yearday,'%.3d'),'/'];
  url = [url,varname,'_',datestr(thedate,'yyyy_mm_dd'),'.nc'];
  
  %%% Path to downloaded file
  fname = './Version4/Release4/nctiles_daily/';
  fname = [fname,varname,'/'];
  fname = [fname,varname,'_',datestr(thedate,'yyyy_mm_dd'),'.nc'];
  
  wgetopts = '--content-disposition --trust-server-names --no-parent --user astewart7583 --password=eOxotM3khABKLd@DMB9';
 
  %%% Increment by one day
  thedate = thedate + 1;
  
  %%% Uncomment to check for files that failed to download and
  %%% re-download
  if (exist(fname))
    continue;
  end
  
  %%% Add to list of wget commands in the script
  fprintf(fid,['wget ',wgetopts,' ',url,'\n']);     
  
end

%%% Close script
fclose(fid);