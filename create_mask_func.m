close all
clear all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

main_folder = '/home/shemesh/julia';
study_folder = '20180720_123744_JH_training_01_1_1';
scan_folder  = [22];
name_of_mask = 'brain_mask';
% name_of_mask = 'rbrain_mask';
% name_of_mask = 'brain_mask_normalized';
numROI       = [1,1,1,1,1,1,1,1]; % ROIs per slice, adapt for number of slices

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

%%% Load images %%%
for p = 1:length(scan_folder)
    % One Nifti with the same geometry (this is just to use the same orientation matrix in the output Nifti header)
    if strcmp(name_of_mask,'brain_mask')
        V_im(p) = spm_vol([main_folder filesep study_folder filesep num2str(scan_folder(p)) filesep 'Processed' filesep 'mean_image.nii']);
        image(:,:,:,p)   =   spm_read_vols(V_im(p));
    elseif strcmp(name_of_mask,'rbrain_mask')
        V_im(p) = spm_vol([main_folder filesep study_folder filesep num2str(scan_folder(p)) filesep 'Processed' filesep 'meanaiImage_' num2str(scan_folder(p)) '_0001.nii']);
        image(:,:,:,p)   =   spm_read_vols(V_im(p));
    elseif strcmp(name_of_mask,'brain_mask_normalized')
        V_im(p) = spm_vol([main_folder filesep study_folder filesep num2str(scan_folder(p)) filesep 'Processed' filesep 'wwmeanaiImage_' num2str(scan_folder(p)) '_0001.nii']);
        image(:,:,:,p)   =   flip(spm_read_vols(V_im(p)),2);
    end
end

%%% Create mask %%%
nx      =   V_im(1).dim(1);
ny      =   V_im(1).dim(2);
nsl     =   V_im(1).dim(3);
roi     =   zeros(nx,ny,nsl);

for i=1:nsl
    one_roi = zeros(nx,ny,numROI(i));
    for j=1:numROI(i)
        for k=1:length(scan_folder)
            figure(1), imagesc(flipdim(imrotate(image(:,:,i,k),-90),2)); colormap(gray); axis image; axis off;
            if ~exist('pos','var')       % Check if position was loaded
                headerMask = impoly;     % Draw ROI       h = impoly begins interactive placement of a polygon on the current axes. The function returns h, a handle to an impoly object.
            else
                headerMask = impoly(gca, pos); shg; % Use existing ROI
            end
            pause;
            % Create mask
            pos = headerMask.getPosition; % Saves points
            one_roi(:,:,j,k) = flipdim(imrotate(headerMask.createMask,90),1); % Creates the mask
        end
    end
    clear pos;
    for k=1:length(scan_folder)
        roi(:,:,i,k) = sum(one_roi(:,:,:,k),3);
    end
end

%%% Save mask %%%
for p = 1:length(scan_folder)
    output_folder   =   [main_folder filesep study_folder filesep num2str(scan_folder(p)) filesep 'Masks' filesep];
    if ~exist('filename','dir')  
    	mkdir(output_folder);        % mkdir folderName creates the folder output_folder.
    end
    if strcmp(name_of_mask,'brain_mask_normalized')
        roi_to_save = flip(roi(:,:,:,p),2);
    else
        roi_to_save = roi(:,:,:,p);
    end
    save([output_folder,name_of_mask,'.mat'],'roi_to_save');
    V_out           =   V_im(p);    
    V_out.fname     =   [output_folder,name_of_mask,'.nii'];
    spm_write_vol(V_out, roi_to_save);
end