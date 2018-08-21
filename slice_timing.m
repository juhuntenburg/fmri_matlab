clear all
close all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

files   =   { ...
     '/Users/francisca/Documents/Data/20180717_120447_CL_180718_visualstim_FFF_1_1/11/Processed'
};
search_name = 'iImage';
TR      =   1; % (seconds)
sliceorder  =   [0 2 4 6 8 10 12 1 3 5 7 9 11]+1;       
refslice    =   1;
backup = 1;

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

%%% Read first dataset %%%
V       =   spm_vol(spm_select('FPList',files{1},['^' search_name '.*\d\d\d\d.nii']));
nx      =   V(1,1).dim(1);
ny      =   V(1,1).dim(2);
nsl     =   V(1,1).dim(3);
NR      =   size(V,1);

%%% Slice timing %%%
 
for i=1:size(files,2)
    % with SPM
    V           =   spm_vol(spm_select('FPList',files{i},['^' search_name '.*\d\d\d\d.nii']));
    vols_pre    =   spm_read_vols(V);
    spm_slice_timing(V, sliceorder, refslice, [TR/nsl TR/nsl], 'a');
end

%%% Backup original files %%%
if backup
    for i=1:size(files,2)
        V           =   spm_vol(spm_select('FPList',files{i},['^a' search_name '.*\d\d\d\d.nii']));
        bkpdir      =   fullfile(fileparts(files{i}),'Processed',['a' search_name]);
        mkdir(bkpdir);
        vols        =   zeros([V(1,1).dim,size(V,1)]);
        for f=1:size(V,1)
           vols(:,:,:,f)    =   spm_read_vols(V(f));
           [path,file,ext]  =   fileparts(V(f).fname);
           V(f).fname       =   fullfile(bkpdir,[file ext]);
           spm_write_vol(V(f),vols(:,:,:,f));
        end
    end
end