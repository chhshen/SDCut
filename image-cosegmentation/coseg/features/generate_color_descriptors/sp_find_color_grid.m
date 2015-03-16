function features=sp_find_color_grid(grid_x,grid_y,I, patchSize,resolution)

num_patches = numel(grid_x);

features.x=grid_x(:) + patchSize/2 - 0.5;
features.y=grid_y(:) + patchSize/2 - 0.5;

% features.x=grid_x(:) + 16/2 - 0.5;
% features.y=grid_y(:) + 16/2 - 0.5;


features.data=zeros(num_patches,3*resolution);

r = patchSize/2;

n_per_bin=(1./resolution);


for i=1:num_patches
    
    cx=features.x(i);
    cy=features.y(i);   

    
    % find window of pixels that contributes to this descriptor
    x_lo = grid_x(i);
    x_hi = grid_x(i) + patchSize - 1;
    y_lo = grid_y(i);
    y_hi = grid_y(i) + patchSize - 1;
    
    
    try
    patch=I(y_lo:y_hi,x_lo:x_hi,:);
    catch 
        keyboard
    end
    
    for channel=1:3
        
       tmp=patch(:,:,channel);
       tmp=floor(tmp/n_per_bin)+1;
       tmp=tmp-(tmp==resolution+1);
       try
       tmp=accumarray(tmp(:),1);
       catch
           keyboard
       end
       tmp=[tmp;zeros(resolution-size(tmp,1),1)];
       
       features.data(i,1+(channel-1)*resolution:channel*resolution)=tmp';
       
       
    end

end

