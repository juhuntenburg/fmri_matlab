clear all
close all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

files   =   { ...
     '/Users/francisca/Documents/Data/20180717_120447_CL_180718_visualstim_FFF_1_1/11/Processed'
};
search_name = 'aiImage';

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

%%% Realign %%%
for i = 1:size(files,2)
     V                  =   spm_vol(spm_select('FPList',files{i},['^' search_name '.*\d\d\d\d.nii']));
     [vox fov]          =   vox_calc(V(1));
     list               =   '';
     for v=1:size(V)
        list    =   [list;V(v).fname]; 
     end
     defaults.realign.estimate.sep         =   (vox(1)*2);
     defaults.realign.estimate.quality     =   1;
     defaults.realign.estimate.fwhm        =   vox(1);
     defaults.realign.estimate.rtm         =   1;
     % Realign
     spm_realign(list,defaults.realign.estimate);
     defaults.realign.write.which          =   2;
     % Reslice
     spm_reslice(list,defaults.realign.write);
end