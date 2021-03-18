%%%
%%% read_nctiles_daily.m
%%%
%%% Convenience function to read daily ECCOv4r4 files, which can't be read
%%% using standard read_nctiles function.
%%%
function field_tiles_out = read_nctiles_daily (file_name,field_name)

  gcmfaces_global;
  
  field_data = ncread(file_name,field_name);
  gs = size(field_data,1);

  field_tiles_out = repmat(0*mygrid.XC,[1 1 size(field_data,4)]);
  field_tiles_out.f1(:,1:gs,:) = squeeze(field_data(:,:,1,:));  
  field_tiles_out.f1(:,gs+(1:gs),:) = squeeze(field_data(:,:,2,:));    
  field_tiles_out.f1(:,2*gs+(1:gs),:) = squeeze(field_data(:,:,3,:));    
  field_tiles_out.f2(:,1:gs,:) = squeeze(field_data(:,:,4,:));  
  field_tiles_out.f2(:,gs+(1:gs),:) = squeeze(field_data(:,:,5,:));    
  field_tiles_out.f2(:,2*gs+(1:gs),:) = squeeze(field_data(:,:,6,:));    
  field_tiles_out.f3(:,1:90,:) = squeeze(field_data(:,:,7,:));    
  field_tiles_out.f4(1:gs,:,:) = squeeze(field_data(:,:,8,:));  
  field_tiles_out.f4(gs+(1:gs),:,:) = squeeze(field_data(:,:,9,:));    
  field_tiles_out.f4(2*gs+(1:gs),:,:) = squeeze(field_data(:,:,10,:));    
  field_tiles_out.f5(1:gs,:,:) = squeeze(field_data(:,:,11,:));  
  field_tiles_out.f5(gs+(1:gs),:,:) = squeeze(field_data(:,:,12,:));    
  field_tiles_out.f5(2*gs+(1:gs),:,:) = squeeze(field_data(:,:,13,:));  
  
end

