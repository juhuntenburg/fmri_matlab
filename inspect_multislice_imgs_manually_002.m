% function inspect_multislice_imgs_manually_002(data,slices,slice_tseries,rows,wind,cmap)
% the function inspect_imgs_manually lets you inspect imgs corresponding to 
% timepoints in a mean time series. You can run through
% timepoints / images sequentially (chronologically). You can 
% also select specific points to inspect. 

% It combines two previous functions "inspect_timeseries_datapoint" 
% (see img of selected point in the time series) and 
% "inspect_video" (run through chronologically over the timeseries). 

% To run through images sequentially, position the cursor in the 
% empty area to the right of the timeseries plot and either press it
% or press the right keyboard arrow

% To select a specific point from the timeseries, position the cursor
% on top of the data point and click.

% Inputs: 
% data:                 4d matrix with nx,ny,slice,frame
% slices:               which slices (dimension 3) you want to check 
%                       if empty ([]) plot all slices
% slice_tseries:        the slice used for calculating the mean 
%                       timeseries (mean over voxels). if [] use 
%                       avg over whole brain 
% wind:                 window, optional, if exists: xlim will be 
%                       [this_frame-wind this_frame+wind]. if empty, xlim
%                       will be [1 total_nr_of_frame]

% Outputs:
% a figure with the images on top, and the mean timeseries
% mean over all non-nan voxels) in the bottom. 

% Madalena Fonseca 
% 05-07-2017

% ----------------
% call create montage directly
% a = create_montage_directly(img,[],4,1)
% ----------------------
study_folder = '20180717_120447_CL_180718_visualstim_FFF_1_1';
scan_folder  = 11; %1st run
other_folder = 'Processed';
nrepetitions = 340;
for i = 1:nrepetitions
    V_im            =   spm_vol(['/Users/francisca/Documents/Data/' study_folder filesep num2str(scan_folder) filesep other_folder filesep 'Image_' num2str(scan_folder) '_' num2str(i,'%04d') '.nii']);
    data(:,:,:,i)   =   flipdim(imrotate(spm_read_vols(V_im),-90),2);
end
slices=[];
slice_tseries=2;
rows=3;
wind=[];
cmap='gray';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dims ax1 ax2 lim h x frame hh lim3  

% params 
%rows            =   5;
%slice_tseries   =   10;

% get data matrix dimensions
dims    =   size(data);

% optional input variables
if ~exist('slices','var') || isempty(slices);
    slices = 1:dims(3);
end

% prepare figure;
fig_pos =[91  115   1155  680];
figure('position',fig_pos,'color','w','name','video'); hold on;
hold on;
%ax1     =   subplot(3,1,1:2);hold on;
ax2     =   subplot(4,1,4);hold on;

% find out position
% temp = figure;
% pos     = get(subplot(4,1,1:3),'position');
% close(temp);
pos = [ 0.1300    0.3291    0.7750    0.5959];

% initialize frame at al and loop over frames
frame   =   1;
while ~isempty(frame)
    
    %  multislice colage 
    % axes(ax1);
    ax1  =  create_montage_directly_int(data,slices,rows,frame);
    %hh  =   imagesc(data(:,:,slice,frame));
    set(ax1,'position',pos);
    set(gca,'yticklabel','','xticklabel',''); box off; axis off; axis image;
    %set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);
    colormap(cmap);
    
    
    % timeseries
    if exist('slice_tseries','var') && ~isempty(slice_tseries);
        roi_temp         =  reshape(squeeze(data(:,:,slice_tseries,:)),[dims(1)*dims(2) dims(4)]);
        roi_mean         =  nanmean(roi_temp,1);
        this_title       =  ['Mean signal avged across voxels (slice ' num2str(slice_tseries) ')'];
    else
        clear roi_temp mean_over_voxels mean_over_slices roi_mean
        roi_temp         =  reshape(data,[dims(1)*dims(2) dims(3) dims(4)]);
        mean_over_voxels =  squeeze(nanmean(roi_temp,1));
        mean_over_slices =  nanmean(mean_over_voxels,1);
        roi_mean         =  mean_over_slices;
        this_title       =  ['Mean signal avged across voxels (whole brain)'];
    end
    
    plot(ax2,1:dims(4),roi_mean,'.-k','linewidth',1);
    lim     =   get(ax2,'ylim');
    h       =   plot(ax2,[frame frame],lim,'r');
    if exist('wind','var') && ~isempty(wind);
        xlim(ax2,[frame-wind frame+wind]);
    else
        xlim(ax2,[1 dims(4)]);
    end
    ylabel(ax2,'Signal (a.u.)');
    xlabel(ax2,'Volume');
    title(ax2,[this_title ' [frame: ' num2str(frame) ']']);
    % set(ax2,'FontName','Myriad Pro');
    
    [x,y]       =   ginput(1);
    x           =   round(x);
    lim3        =   get(gca,'xlim');
    if x>lim3; 
        frame   =   frame+1;
    else
        frame   =    x;
    end;
    
    delete(h);
end

% end


function ax = create_montage_directly_int(data,slice_range,rows,frame)
% create imags montages 
% Inputs:
% data: 4d matrix
% slice_range that you want to plot 
% how many rows in the final panel?

% eg.
% data          = img(x,y,nslices,frames)   % 4D matrix
% slice_range   = 11:20;                    % plot from slice 11 to 20
% rows          = 2;                        % make a 2 x 5 mosaic
% frame         = 10;                       % 10th images
% ax = create_montage_directly(data,slice_range,rows,frame)

%close all
%clear dir

dims        = size(data);

% optional input variables
if ~exist('slice_range','var') || isempty(slice_range);
    slice_range = 1:dims(3);
end
if ~exist('rows','var') || isempty(rows);
    rows = 2;  
end
if ~exist('frame','var') || isempty(frame);
    frame = 1;
    warning('Using 1st frame');
end

% calculate cols 
cols        = ceil(numel(slice_range)/rows);

% prepare
SuperImage = cell(rows, cols);

% initialize variables
r   = 1; 
c   = 1;
idx = 0;

% Keep looping until go over the number of mosaics 
while c<=cols && r<=rows 
    idx = idx+1;
    if ~(idx>numel(slice_range)) % otherwise leave empty
        SuperImage{r,c} = data(:,:,slice_range(idx),frame);
    else
        SuperImage{r,c} = NaN(dims(1),dims(2));
    end
    if c == cols
        c = 1;
        r = r+1;
    else
        c = c+1;
        r = r;
    end
    if idx>(cols*rows)
        c = cols+1; % stop loop
    end
end

% convert to mat
SuperImage2 = cell2mat(SuperImage);

% plot in montage format 
%figure, 
ax = axes;
h=imagesc(SuperImage2(:,:),[min(data(:)) max(data(:))].*0.9);  
axis image; axis off; axis fill; colormap(gray(256)) %%% without mask
axis off; set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);
%colormap(hot)
end
      
