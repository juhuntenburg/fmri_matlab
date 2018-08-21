clear all
close all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

main_folder = '/Users/francisca/Documents/Data';
study_folder = '20180717_120447_CL_180718_visualstim_FFF_1_1';
scan_folder = 6;
nslices = 20;

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

% Read images
input_folder = [main_folder filesep study_folder filesep num2str(scan_folder) filesep 'pdata' filesep '1'];
V = spm_vol([input_folder filesep 'Image_0001.nii']);
vols = flipdim(imrotate(spm_read_vols(V),-90),2);

% Load mask
m = spm_vol([main_folder filesep study_folder filesep num2str(scan_folder) filesep 'brain_mask.nii']);
mask = logical(flipdim(imrotate(spm_read_vols(m),-90),2));

% Apply mask to realigned images
masked = vols;
masked(~mask) = nan;
for j = 1:nslices
    figure(1),imagesc(masked(:,:,j)); axis image; axis off; colormap(gray(256)); colorbar; shg; pause(1);
end

% Save masked image
V_out = V;
V_out.fname = fullfile(input_folder,'cImage_0001.nii');
spm_write_vol(V_out, flipdim(imrotate(squeeze(masked),90),1));