clear all
close all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

study_folder = {'/Users/francisca/Documents/Data/20180717_120447_CL_180718_visualstim_FFF_1_1'};
scan_folder = [11];
masking_folder = [11];
name_of_mask = 'brain_mask';
% name_of_mask = 'rbrain_mask';
no_folders = [1];
before_correction = 0;
during_correction = 0;
after_correction = 1;
dbstop at 198 if during_correction==1
dbstop at 306 if during_correction==1
view_mask_applied = 0;
use_image_detrending = 0;
change_std_threshold = 0;
apply_std_threshold_for_correction = 0;
nrepetitions = 340;
nslices = 13;
event = [zeros(1,40) ones(1,20) zeros(1,40) ones(1,20) zeros(1,40) ones(1,20) zeros(1,40) ones(1,20) zeros(1,40) ones(1,20) zeros(1,40)];
onsets = 41:60:281;
dur = repmat(20,[1,5]);
shift = 20;

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

k=1;
for f = 1:length(study_folder)
    for s = 1:no_folders(f)
        %%% Images and mask %%%
        input_folder = [study_folder{f} filesep num2str(scan_folder(k)) filesep 'Processed'];
        mask_folder = [study_folder{f} filesep num2str(masking_folder(k)) filesep 'Masks' filesep name_of_mask '.nii'];

        %%% Read images %%%
        for i = 1:nrepetitions
            if strcmp(name_of_mask,'brain_mask')
                if after_correction
                    file = [input_folder  filesep 'iImage_' num2str(scan_folder(k)) '_' num2str(i,'%04d') '.nii'];
                else
                    file = [input_folder  filesep 'Image_' num2str(scan_folder(k)) '_' num2str(i,'%04d') '.nii'];
                end
            elseif strcmp(name_of_mask,'rbrain_mask')
                file = [input_folder  filesep 'raiImage_' num2str(scan_folder(k)) '_' num2str(i,'%04d') '.nii'];
            end
            V = spm_vol(file);
            vols(:,:,:,i) = flipdim(imrotate(spm_read_vols(V),-90),2);
        end

        %%% Load mask %%%
        m = spm_vol(mask_folder);
        mask = logical(flipdim(imrotate(spm_read_vols(m),-90),2));
        mask = repmat(mask,1,1,1,nrepetitions);

        %%% Apply mask to images %%%
        masked = vols;
        masked(~mask) = nan;
        if view_mask_applied
            for j = 1:nslices
               figure(1),imagesc(masked(:,:,j,60)); axis image; axis off; colormap(gray(256)); colorbar; shg; pause(1);
            end
        end

        permuteOrder = [3 4 1 2];
        nx = size(vols,1);
        ny = size(vols,2);
        vols_new=vols;

        %%% Show voxels timecourse (in chosen ROI) %%%
        help_matrix = permute(masked,permuteOrder);
        final_matrix = help_matrix(:,:,:);
        which = find(sum(~isnan(squeeze(final_matrix(:,1,:))),2)>0); % slices of roi

        vols_matrix = permute(vols,permuteOrder);  % not masked
        final_vols = vols_matrix(:,:,:);

        for sl = which'
            %%% Show voxels average temporal evolution (in chosen ROI) %%%
            continuing=1;
            to_correct=[];
            while continuing
                % Apply mask to new images
                masked = vols_new;
                masked(~mask) = nan;

                % Show voxels timecourse (in chosen ROI)
                help_matrix = permute(masked,permuteOrder);
                final_matrix = help_matrix(:,:,:);

                close all
                t=figure(sl);
                set(t,'position',[1 56 1280 649]);

                % Brain with ROI contour
                brain_with_roi = subplot(5,8,[1 2 9 10]); 
                imagesc(vols(:,:,sl,100)); colormap(gray(256));set(gca,'yticklabel','','xticklabel',''); axis image; axis off; hold on;
                zed_mask = squeeze(mask(:,:,:,1));
                [~,h1] = contour(zed_mask(:,:,sl),[1 1],'r');set(h1,'LineWidth',2); hold off;

                % Calculate avg over voxels inside roi
                roi_temp                =  squeeze(final_matrix(sl,:,:));
                roi_mean(sl,:)          =  nanmean(roi_temp,2);

                % Separate cycles
                trange   =   (onsets(1):nrepetitions)-shift;
                subplot(5,8,[4 5 12 13]);
                cycle    =   reshape(roi_mean(sl,trange),[numel(trange)/numel(onsets),numel(onsets)]);
                plot((0:size(cycle,1)-1)-shift,cycle,'.-','linewidth',1);
                lim = get(gca,'ylim');
                patch([0;dur(1);dur(1);0],[lim(1);lim(1);lim(2);lim(2)],[0.3 0.5 0.9], 'EdgeColor','none','FaceAlpha',0.2);
                title('Cycle timecourse of voxels inside ROI');
                ylabel('Signal (a.u.)'); xlabel('Volume');

                % Average cycle
                subplot(5,8,[7 8 15 16]);
                avg_cycles = nanmean(cycle,2); % average cycle
                sem_cycles = nanstd(cycle,[],2)./sqrt(sum(~isnan(cycle),2)); % standard error of the mean
                xxx  = [(0:size(cycle,1)-1)'-shift; flipud((0:size(cycle,1)-1)'-shift)];
                yyy  = [(avg_cycles+sem_cycles); flipud(avg_cycles-sem_cycles)];
                plot((0:size(cycle,1)-1)-shift,avg_cycles,'color','k','linewidth',1.2);
                patch(xxx,yyy,[0.4 0.1 0.1], 'EdgeColor','none','FaceAlpha',0.2);
                patch([0;dur(1);dur(1);0],[lim(1);lim(1);lim(2);lim(2)],[0.3 0.5 0.9], 'EdgeColor','none','FaceAlpha',0.2);
                title('Average cycle timecourse of voxels inside ROI');
                ylabel('Signal (a.u.)'); xlabel('Volume');

                % Plot avg over voxels inside roi
                tplot = subplot(5,8,[25:32 33:40]);
                plot(roi_mean(sl,:),'color',[0.6 0.6 0.6],'linewidth',0.5); xlim([0 nrepetitions]); hold on;
                title('Timecourse of voxels inside ROI');
                ylabel('Signal (a.u.)'); xlabel('Volume');
                lim =   [min(roi_mean(sl,:)) max(roi_mean(sl,:))];              % patch
                x   =   [onsets-1;onsets-1+dur]; xx =     [x;flipud(x)];        % patch
                yy  =   repmat([lim(1);lim(1);lim(2);lim(2)],1,numel(onsets));  % patch
                if size(xx)~=size(yy)
                    yy  =   repmat([lim(1);lim(1);lim(2);lim(2)],numel(onsets),1);% patch
                end
                patch(xx,yy,[0.3 0.5 0.9], 'EdgeColor','none','FaceAlpha',0.2);   % patch
                set(gca,'xlim',[0 nrepetitions],'ylim',lim);                         % correct axis
                if before_correction
                    saveas(gcf,[study_folder{f} filesep num2str(scan_folder(k)) filesep 'Processed' filesep 'slice_avg_' num2str(sl) '_no_correction'],'tif')
                elseif after_correction
                    saveas(gcf,[study_folder{f} filesep num2str(scan_folder(k)) filesep 'Processed' filesep 'slice_avg_' num2str(sl) '_correction'],'tif')
                elseif strcmp(name_of_mask,'rbrain_mask')
                    saveas(gcf,[study_folder{f} filesep num2str(scan_folder(k)) filesep 'Processed' filesep 'slice_avg_' num2str(sl) '_realigned'],'tif')
                end

                if isempty(to_correct)
                    % Adjust curve and standard deviation
                    ROI_mean_smooth = smooth(roi_mean(sl,:),5)';
                    x=1:nrepetitions;
                    p = polyfit(x,ROI_mean_smooth,2);
                    plot(x,polyval(p,x));
                    thres   =   1;
                    a       =   NaN(2,1);
                    a(1)    =   plot(polyval(p,x)+(thres*std(roi_mean(sl,:))));
                    a(2)    =   plot(polyval(p,x)-(thres*std(roi_mean(sl,:))));

                    if change_std_threshold
                        go      =   input('Happy with standard deviation threshold = 1? Press enter. Otherwise, enter new threshold. \n');
                        while ~isempty(go)
                            thres   =  go;
                            delete(a);
                            a       =   NaN(2,1);
                            a(1)    =   plot(polyval(p,x)+(thres*std(roi_mean(sl,:))));
                            a(2)    =   plot(polyval(p,x)-(thres*std(roi_mean(sl,:))));
                            go      =   input('Happy? Press enter. Otherwise, enter new threshold. \n');
                        end
                    end

                    if apply_std_threshold_for_correction
                        points = find(roi_mean(sl,:)<polyval(p,x)-thres*std(roi_mean(sl,:)) | roi_mean(sl,:)>polyval(p,x)+thres*std(roi_mean(sl,:)));
                    else
                        points = [];
                    end
                else
                    points = [];
                end

                clear h_indexedit
                b='h_indexedit';

                axes_lim = get(gca, 'YLim');
                axes_height = [axes_lim(1) axes_lim(2)];
                h_indextext = uicontrol(gcf, 'Units', 'characters', 'Position', [15 22 12 1],...
                'String', 'Frames to correct: ', 'Style', 'text', ...
                'HorizontalAlignment', 'left');
                pos = [30 22 100 1];
                h_selectframe = uicontrol(gcf, 'Units', 'characters', 'Position', [135 22 12 1],...
                    'String', 'Select frame', 'Style', 'pushbutton', ...
                    'Callback', 'if exist(b) points=str2num(h_indexedit.String); end; [points,h_indexedit]=updatee(points,axes_height,pos,to_correct)');
                h_finish = uicontrol(gcf, 'Units', 'characters', 'Position', [150 22 12 1],...
                    'String', 'Finish', 'Style', 'pushbutton', ...
                    'Callback', 'if exist(b) points=str2num(h_indexedit.String); end');

                if length(points)<=length(to_correct)
                    continuing=0;
                else
                    to_correct=points;
                    good_volumes = setdiff(1:nrepetitions,to_correct);
                    vq = interp1(good_volumes,squeeze(final_vols(sl,good_volumes,:)),1:nrepetitions,'spline');
                    fl = reshape(vq(to_correct,:),[size(to_correct,2),nx,ny]);
                    fl = permute(fl,[2,3,1]);
                    vols_new(:,:,sl,to_correct) = fl;
                    save([study_folder{f} filesep num2str(scan_folder(k)) filesep 'Processed' filesep 'to_correct' num2str(sl)],'to_correct');
                end
            end

            %%% Show voxels temporal evolution (in chosen ROI) %%%
            continuing=1;
            while continuing
                % Apply mask to new images
                masked = vols_new;
                masked(~mask) = nan;

                % Show voxels timecourse (in chosen ROI)
                help_matrix = permute(masked,permuteOrder);
                final_matrix = help_matrix(:,:,:);

                % Choose voxels that belong to the mask
                temp_matrix = squeeze(final_matrix(sl,:,:))';
                f_m = temp_matrix(find(sum(isnan(temp_matrix),2)==0),:);
                % f_m = temp_matrix(~isnan(temp_matrix(:,1)),:);   % Other way to create f_m

                % Sort voxels according to their intensity in a specific TR (in this case, the 1st)
                [~,I] = sort(f_m(:,1));
                f_m_order = f_m(I,:);

                close all
                t=figure(sl);
                set(t,'position',[1 56 1280 649]);

                subplot(6,1,2:6),
                imagesc(f_m_order); colormap(jet); colorbar; shg % ROI brain
                xlabel('Volume'); ylabel('Voxels');
                title('Voxels timecourse (raw)');
                if before_correction
                    saveas(gcf,[study_folder{f} filesep num2str(scan_folder(k)) filesep 'Processed' filesep 'slice_' num2str(sl) '_no_correction'],'tif')
                elseif after_correction
                    saveas(gcf,[study_folder{f} filesep num2str(scan_folder(k)) filesep 'Processed' filesep 'slice_' num2str(sl) '_correction'],'tif')
                elseif strcmp(name_of_mask,'rbrain_mask')
                    saveas(gcf,[study_folder{f} filesep num2str(scan_folder(k)) filesep 'Processed' filesep 'slice_' num2str(sl) '_realigned'],'tif')
                end

                clear h_indexedit
                b='h_indexedit';
                points = [];
                axes_lim = get(gca, 'YLim');
                axes_height = [axes_lim(1) axes_lim(2)];
                h_indextext = uicontrol(gcf, 'Units', 'characters', 'Position', [15 37 12 1],...
                'String', 'Frames to correct: ', 'Style', 'text', ...
                'HorizontalAlignment', 'left');
                pos = [30 37 80 1];
                h_selectframe = uicontrol(gcf, 'Units', 'characters', 'Position', [135 37 12 1],...
                    'String', 'Select frame', 'Style', 'pushbutton', ...
                    'Callback', 'if exist(b) points=str2num(h_indexedit.String); end; [points,h_indexedit]=updatee(points,axes_height,pos,to_correct)');
                h_finish = uicontrol(gcf, 'Units', 'characters', 'Position', [150 37 12 1],...
                    'String', 'Finish', 'Style', 'pushbutton', ...
                    'Callback', 'if exist(b) points=str2num(h_indexedit.String); end');

                if use_image_detrending
                    % 1: just removing the drift and mean by a least-squares fit of a straight line to the data and the subtraction of the resulting function from the data
                    r0          =   ones(1,nrepetitions);
                    r1          =   linspace(0,1,nrepetitions);
                    b=[r0;r1]'\f_m_order';
                    pred=[r0;r1]'*b;
                    resid1=f_m_order-pred';

                    %     for i = 1:size(f_m_order,1)            % same
                    %         resid1(i,:)= detrend(f_m_order(i,:));
                    %     end

                    % 2: removing the drift and mean using a moving average filter with
                    % a span of nrepetitions, and transforming values into a standard normal distribution
                    for i = 1:size(f_m_order,1)
                        resid2(i,:)= zscore(f_m_order(i,:)-smooth(f_m_order(i,:),nrepetitions)');
                    end

                    % 3: just removing the drift and mean using a moving average filter with
                    % a span of nrepetitions
                    for i = 1:size(f_m_order,1)
                        resid3(i,:)= f_m_order(i,:)-smooth(f_m_order(i,:),nrepetitions)';
                    end

                    figure(sl+1),
                    imagesc(resid1); colormap(jet); colorbar; shg % ROI brain
                    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                    xlabel('Volume'); ylabel('Voxels');
                    title('Voxels timecourse after mean and trend removal (using least-squares linear fit)');

                    figure(sl+2);
                    imagesc(resid2); colormap(jet); colorbar; shg % ROI brain
                    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                    xlabel('Volume'); ylabel('Voxels');
                    title('Voxels timecourse after mean and trend removal (using moving average filter) and z-score calculation');

                    figure(sl+3);
                    imagesc(resid3); colormap(jet); colorbar; shg % ROI brain
                    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                    xlabel('Volume'); ylabel('Voxels');
                    title('Voxels timecourse after mean and trend removal (using moving average filter)');
                end

                if length(points)<=length(to_correct)
                    continuing=0;
                else
                    to_correct=points;
                    good_volumes = setdiff(1:nrepetitions,to_correct);
                    vq = interp1(good_volumes,squeeze(final_vols(sl,good_volumes,:)),1:nrepetitions,'spline');
                    fl = reshape(vq(to_correct,:),[size(to_correct,2),nx,ny]);
                    fl = permute(fl,[2,3,1]);
                    vols_new(:,:,sl,to_correct) = fl;
                    save([study_folder{f} filesep num2str(scan_folder(k)) filesep 'Processed' filesep 'to_correct' num2str(sl)],'to_correct');
                end
            end
        end

        %%% Save corrected images
        if during_correction
            V_out           =   V;
            for i = 1:nrepetitions
                V_out.fname = fullfile(input_folder,['iImage_' num2str(scan_folder(k)) '_' num2str(i,'%04d') '.nii']);
                spm_write_vol(V_out, flipdim(imrotate(squeeze(vols_new(:,:,:,i)),90),1));
            end
        end
        k=k+1;
    end
end