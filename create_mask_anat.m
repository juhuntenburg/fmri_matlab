close all
clear all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

main_folder = '/Users/francisca/Documents/Data';
study_folder = '20180717_120447_CL_180718_visualstim_FFF_1_1';
scan_folder  = 6;
name_of_mask = 'brain_mask';
numROI       = ones(1,20); % ROIs per slice

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

%%% Load images %%%
Bruker2nifti([main_folder filesep study_folder filesep num2str(scan_folder) filesep 'pdata' filesep '1' filesep '2dseq']);
V_im            =   spm_vol([main_folder filesep study_folder filesep num2str(scan_folder) filesep 'pdata' filesep '1' filesep 'Image_0001.nii']); % One Nifti with the same geometry (this is just to use the same orientation matrix in the output Nifti header)  
image   =   spm_read_vols(V_im);

%%% Create mask %%%
nx      =   V_im.dim(1);
ny      =   V_im.dim(2);
nsl     =   V_im.dim(3);
roi     =   zeros(nx,ny,nsl);

for i=1:nsl
    one_roi = zeros(nx,ny,numROI(i));
    for j=1:numROI(i)
        figure(1), imagesc(flipdim(imrotate(image(:,:,i),-90),2)); colormap(gray); axis image; axis off; colorbar;
        one_roi(:,:,j)    =   flipdim(imrotate(roipoly, 90),1);  
    end
    roi(:,:,i) = sum(one_roi,3);
    figure(2), imagesc(flipdim(imrotate(roi(:,:,i),-90),2)); colormap(gray); axis image; axis off;
end

%%% Save mask %%%
output_folder   =   [main_folder filesep study_folder filesep num2str(scan_folder) filesep];
save([output_folder,name_of_mask,'.mat'],'roi');
V_out           =   V_im;    
V_out.fname     =   [output_folder,name_of_mask,'.nii'];
spm_write_vol(V_out, roi);