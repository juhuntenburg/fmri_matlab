clear all
close all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

main_folder = '/Users/francisca/Documents/Data';
study_folder = '20180717_120447_CL_180718_visualstim_FFF_1_1';
nrepetitions = 340;
nslices = 13;
scan_folder = [11];
search_name='aiImage';
view_mask_applied = 0;
apply_to_mean_image=1;
apply_to_all_functional_images=0;

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

%%% For the mean image %%%

if apply_to_mean_image
    for m=1:length(scan_folder)
        scan_f=scan_folder(m);
        % Read image
        input_folder = [main_folder filesep study_folder filesep num2str(scan_f) filesep 'Processed'];
        V = spm_vol([input_folder filesep 'mean' search_name '_' num2str(scan_f) '_0001.nii']);
        vols = flipdim(imrotate(spm_read_vols(V),-90),2);

        % Load mask
        m = spm_vol([main_folder filesep study_folder filesep num2str(scan_f) filesep 'Masks' filesep 'rbrain_mask.nii']);
        mask = logical(flipdim(imrotate(spm_read_vols(m),-90),2));

        % Apply mask to mean image
        masked = vols;
        masked(~mask) = nan;

        if view_mask_applied
            for j = 1:nslices
                figure(1),imagesc(masked(:,:,j)); axis image; axis off; colormap(gray(256)); colorbar; shg; pause(1);
            end
        end

        % Save masked image
        V_out = V;
        V_out.fname = fullfile(input_folder,['cmean' search_name '_' num2str(scan_f) '_0001.nii']);
        spm_write_vol(V_out, flipdim(imrotate(squeeze(masked),90),1));
    end
end

%%% For the single volumes %%%

if apply_to_all_functional_images
    for m=1:length(scan_folder)
        scan_f=scan_folder(m);
        %%% Read realigned images %%%
        input_folder = [main_folder filesep study_folder filesep num2str(scan_f) filesep 'Processed'];
        for i = 1:nrepetitions
            file = fullfile(input_folder,['r' search_name '_' num2str(scan_f) '_' num2str(i,'%04d') '.nii']);
            V = spm_vol(file);
            vols(:,:,:,i) = flipdim(imrotate(spm_read_vols(V),-90),2);
        end

        %%% Load mask %%%
        m = spm_vol([main_folder filesep study_folder filesep num2str(scan_f) filesep 'Masks' filesep 'rbrain_mask.nii']);
        mask = logical(flipdim(imrotate(spm_read_vols(m),-90),2));
        mask = repmat(mask,1,1,1,nrepetitions);

        %%% Apply mask to realigned images %%%
        masked = vols;
        masked(~mask) = nan;
        
        if view_mask_applied
            for j = 1:nslices
                figure(1),imagesc(vols(:,:,j,1)); axis image; axis off; colormap(gray(256)); colorbar; shg; pause(1);
                figure(1),imagesc(masked(:,:,j,1)); axis image; axis off; colormap(gray(256)); colorbar; shg; pause(1);
            end
        end

        %%% Save masked images %%%
        V_out = V;
        for i = 1:nrepetitions
            V_out.fname = fullfile(input_folder,['cr' search_name '_' num2str(scan_f) '_' num2str(i,'%04d') '.nii']);
            spm_write_vol(V_out, flipdim(imrotate(squeeze(masked(:,:,:,i)),90),1));
        end
    end
end