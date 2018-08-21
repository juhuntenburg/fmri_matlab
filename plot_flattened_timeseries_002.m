% function [ax] = plot_flattened_timeseries_002(mat,mask,event,choose_colormap)
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% This function plots a 2D heatmap of every timeseries within an
% image, with frames in the xaxis, voxels in the y-axis and intensity 
% as color. 

% Input:
%   mat:        4-D, data matrix with dimensions: nx,ny,nslice,frames 
%   optional:   mask: brain mask with dimensions (nx,ny,nslices),
%               vector with event as zeros and ones colormap (eg. 'jet')

% Output: 
% Figure with voxels in xaxis, frames in yaxis and intensity as color
%   (trend and mean removal, voxel by voxel)

% Example with all optional inpusts:
%   plot_flattened_timeseries_002(struct.data,struct.mask,struct.events.Licks,'jet');

% Example with all inputs except brain mask: 
%   plot_flattened_timeseries_002(struct.data,[],struct.events.Licks,'jet');

% Based on Power 2016 NeuroImage
% MF, Mainen/Shemesh Lab, 2017
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
main_folder = '/home/shemesh/julia'
study_folder = '20180720_123744_JH_training_01_1_1';
scan_folder  = 22; %1st run
image_folder = [main_folder filesep study_folder filesep num2str(scan_folder) filesep 'Processed']
brain_mask = [main_folder filesep study_folder filesep num2str(scan_folder) filesep 'Masks' filesep 'brain_mask.mat']
nrepetitions = 340;
mask = load(brain_mask);
mask = mask.roi_to_save;
% event = [zeros(1,40) ones(1,20) zeros(1,40) ones(1,20) zeros(1,40) ones(1,20) zeros(1,40) ones(1,20) zeros(1,40) ones(1,20) zeros(1,40)];
choose_colormap = 'jet';
for i = 1:nrepetitions
    V_im            =   spm_vol([image_folder filesep 'Image_' num2str(i,'%04d') '.nii']);
    mat(:,:,:,i)   =   spm_read_vols(V_im);
end
% resp_rate = load('/Users/francisca/Documents/Data/FFF_SAM-Data/August_2017/170807r_rate.mat');
% r_rate = resp_rate.r_rate(:,1); %1st run
% temp = load('/Users/francisca/Documents/Data/FFF_SAM-Data/August_2017/170807tem.mat');
% tem = temp.tem(:,1); %1st run
% 
clear dims flat idx ll 

% dimensions
dims        =   size(mat);

% apply mask if available. 
%if exist('event','var') && ~isempty(mask)
mask(mask==0)   = NaN;
mat             = mat.*repmat(mask,[1 1 1 dims(4)]);
%end

% flatten matrix 
flat        =   reshape(mat,[dims(1)*dims(2)*dims(3) dims(4)]);

% get the non-nan voxels 
idx         =   find(~isnan(flat(:,dims(4)/2))); % non-nan voxels; arbitrary timepoint

% build slice label
ll          =   NaN(dims(1),dims(2),dims(3),1);
for i = 1:dims(3)
    ll(:,:,i)  =  repmat(i,dims(1),dims(2));
end
ll          =  reshape(ll,[dims(1)*dims(2)*dims(3) 1]);

% slice transition
idx_slice_transition  =   find(diff(ll(idx,1))>0);

% mean and trend regressors
r0      =   ones(1,dims(4));
r1      =   linspace(0,1,dims(4));

% set up a temporal mask ignoring the first x volumes
exclude_first               =    0;
exclude_last                =    0;
tmask                       =   r0;
tmask(1:exclude_first)      =   0;
tmask(end+1-exclude_last:end) =   0;

% take out mean and trend (as in Power 2016, NeuroImage)
clear flat1
flat1               =   flat;
% presumes vox x time input variable structure
r=[r0;r1];
b=r(:,~~tmask)'\flat(:,~~tmask)'; % right matrix / left matrix
pred=r'*b;
resid=flat-pred';
flat1 = [resid pred' b'];

% figure
fig_pos = [138    77   994   628]; 
figure('position',fig_pos,'color','w');
% 
% % % plot an event (paradigm, motion parameters, odor, licks, etc)
% if exist('event','var') && ~isempty(event)
%      subplot(7,1,1); hold on; 
%      plot([find(event>0);find(event>0)],repmat([0;1],[1,100]),'color','k','linewidth',2);
%      ylabel('Events');
%      set(gca,'xlim',[0 size(flat1,2)],'TickDir', 'Out', 'TickLength', [0.01 0.025]); box off; 
% end
% 
% subplot(7,1,2); hold on; 
% plot(r_rate,'color','b','linewidth',2);
% ylabel('Resp rate (Resp/m)');
% set(gca,'xlim',[0 size(flat1,2)],'TickDir', 'Out', 'TickLength', [0.01 0.025]); box off; 
% 
% subplot(7,1,3); hold on; 
% plot(tem,'color','b','linewidth',2);
% ylabel('T (\circC)');
% set(gca,'xlim',[0 size(flat1,2)],'TickDir', 'Out', 'TickLength', [0.01 0.025]); box off; 

% voxel x timeseries 
ax = subplot(7,1,4:7);
hold on;
imagesc(flat1(idx,:));
ylabel('Voxel','Fontsize',12); 
xlabel('Frame','Fontsize',12);
set(gca,'TickDir', 'Out', 'TickLength', [0.01 0.025]); 
box off; colorbar; %colorbar('position',[0.13*7 0.11 0.02 0.53])

for iT = 1:numel(idx_slice_transition);
    hold on;
    plot(1:dims(4),repmat(idx_slice_transition(iT),1,dims(4)),'color',[0.6 0.6 0.6]);
end
xlim([0 dims(4)]);ylim([0 size(ll(idx,1),1)]);

if exist('choose_colormap','var') 
    colormap(choose_colormap);
else
    colormap(jet);
end
caxis([0 100]);
% saveas(gcf,'global_view','fig');

% end


% function [resid pred b]=myregress2(r,tc,varargin)
% 
% % presumes vox x time input variable structure
% clear b pred resid
% 
% if isempty(varargin)
%     % use all timepoints
%     b=r'\tc';
%     pred=r'*b;
%     resid=tc-pred';
%     
% else
%     % use only specified timepoints
%     tmask=varargin{1,1};
%     b=r(:,tmask)'\tc(:,tmask)'; % right matrix / left matrix
%     pred=r'*b;
%     resid=tc-pred';
% end
% end