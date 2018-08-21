close all
clear all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

main_folder = '/home/shemesh/julia';
study_folder = '20180720_123744_JH_training_01_1_1';
scan_folder = [22];
move_to_processed = 1;

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

for i=1:length(scan_folder)
    Bruker2nifti([main_folder filesep study_folder filesep num2str(scan_folder(i)) filesep 'pdata' filesep '1' filesep '2dseq']);
end

if move_to_processed == 1
    origin = [main_folder filesep study_folder filesep num2str(scan_folder(i)) filesep 'pdata' filesep '1' filesep '*.nii'];
    destin = [main_folder filesep study_folder filesep num2str(scan_folder(i)) filesep 'Processed'];
    movefile(origin,destin);
end