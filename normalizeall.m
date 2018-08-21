clear all
close all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

spmmouse('load',['/Users/francisca/Documents/MATLAB/code' filesep 'modified' filesep 'mouse-C57.mat']);
global defaults;
file_fun = '/Users/francisca/Documents/Data/20180717_120447_CL_180718_visualstim_FFF_1_1/11/Processed/wraiImage_11_0001.nii'; % One functional image that aligns all the others!
file_static  = '/Users/francisca/Documents/Data/Atlas_ROI/2012_version/Atlas_axial.nii';  % RARE/Atlas that was selected to align all animals to
nrepetitions = 340;

proc = {
    '/Users/francisca/Documents/Data/20180717_120447_CL_180718_visualstim_FFF_1_1/11/Processed'
    };

file_mov = {
    '/Users/francisca/Documents/Data/20180717_120447_CL_180718_visualstim_FFF_1_1/6/pdata/1/cImage_0001.nii'
    };

file_anat = {
    '/Users/francisca/Documents/Data/20180717_120447_CL_180718_visualstim_FFF_1_1/6/pdata/1/Image_0001.nii'
    };
    
no_folders = [1];
search_name='aiImage';

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

image_static = spm_vol(file_static);
image_fun = spm_vol(file_fun);

[bbstatic,~] = bbvox_from_V(image_static);
[bbfun,vxfun] = bbvox_from_V(image_fun);
rxfun = vxfun(1);
ryfun = vxfun(2);
rzfun = vxfun(3);

m = 1;

for i = 1:length(file_mov)
    image_mov = spm_vol(file_mov{i});
    image_anat = spm_vol(file_anat{i});
    
    [bbmov,vxmov] = bbvox_from_V(image_mov);
    rx = vxmov(1); 
    ry = vxmov(2);
    rz = vxmov(3);

    dnrm = defaults.normalise;
    dnrm.write.vox = [rx,ry,rz];
    % [bb2,~] = spm_get_bbox(image_mov,'fv',eye(4));  
    % center = (bb2(1,3)+bb2(2,3))/2;
    % center2 = (bbmov(1,3)+bbmov(2,3))/2;
    dnrm.write.bb = bbmov;
    % dnrm.write.bb(:,3) = bbmov(:,3)+repmat(center-center2,2,1);
    dnrm.write.bb(1,2) = -bbmov(2,2);
    dnrm.write.bb(2,2) = -bbmov(1,2);
%     dnrm.write.bb(1,3) = -bbmov(2,3);      % Also needed if using 2012 brain atlas
%     dnrm.write.bb(2,3) = -bbmov(1,3);      % Also needed if using 2012 brain atlas

    
    %%%%%%%%%% Rigid Body Transformation %%%%%%%%%%
    flags = defaults.coreg.estimate;
    flags.sep = [abs(rx)*2 abs(rx)*1];
    flags.fwhm = [abs(rx)*4 abs(rx)*4];
    x = spm_coreg(image_static,image_mov,flags);
    % x = zeros(1,6); % If I don't want to apply any transformation
    Affine = image_mov.mat\spm_matrix(x(:)')*image_static.mat;
    Tr = [];
    VG = image_static; 
    VF = image_mov;
    cd(fileparts(file_mov{i}));
    save([fileparts(file_mov{i}) filesep 'coreg1norm.mat'],'Affine','VG', 'VF', 'Tr');
    % spm_write_sn(image_mov.fname, [fileparts(file_mov{i}) filesep 'coreg1norm.mat'], dnrm.write);
    % spm_write_sn(image_anat.fname, [fileparts(file_mov{i}) filesep 'coreg1norm.mat'], dnrm.write); 
    
    %%%%%%%%%% Affine Transformation %%%%%%%%%%
    flags.WF = [];
    flags.WG = [];
    flags.sep = flags.sep(1);
    flags.regtype = 'subj';
    flags.debug = 0;
    [xaff1,~] = spm_affreg(image_static,image_mov,flags,inv(spm_matrix(x(:)')),1); 
    % xaff1 = eye(4); % If I don't want to apply any transformation
    Affine = image_mov.mat\inv(xaff1)*image_static.mat;
    Tr = [];
    VG = image_static; 
    VF = image_mov;
    cd(fileparts(file_mov{i}));
    save([fileparts(file_mov{i}) filesep 'affreg1norm.mat'],'Affine','VG', 'VF', 'Tr');
    spm_write_sn(image_mov.fname, [fileparts(file_mov{i}) filesep 'affreg1norm.mat'], dnrm.write);
    spm_write_sn(image_anat.fname, [fileparts(file_mov{i}) filesep 'affreg1norm.mat'], dnrm.write);

    %%%%%%%%%% Apply transformation to functional images %%%%%%%%%%
    dnrm = defaults.normalise;
    dnrm.write.vox = [rxfun,ryfun,rzfun];
    % [bb2,~] = spm_get_bbox(image_fun,'fv',eye(4));  
    % center = (bb2(1,3)+bb2(2,3))/2;
    % center2 = (bbfun(1,3)+bbfun(2,3))/2;
    dnrm.write.bb = bbfun;
    % dnrm.write.bb(:,3) = bbfun(:,3)+repmat(center-center2,2,1);
    dnrm.write.bb(1,2) = -bbfun(2,2);
    dnrm.write.bb(2,2) = -bbfun(1,2);
%     dnrm.write.bb(1,3) = -bbfun(2,3);      % Also needed if using 2012 brain atlas
%     dnrm.write.bb(2,3) = -bbfun(1,3);      % Also needed if using 2012 brain atlas

    for k = 1:no_folders(i) 
        out_folder = proc{m};
        V_im = spm_vol(spm_select('FPList',out_folder,['^wr' search_name '.*\d\d\d\d.nii']));
        V_mean = spm_vol(spm_select('FPList',out_folder,['^wmean' search_name '.*\d\d\d\d.nii']));
        save([out_folder filesep 'affreg1norm.mat'],'Affine','VG', 'VF', 'Tr');
        for j=1:nrepetitions
            spm_write_sn(V_im(j),[proc{m} filesep 'affreg1norm.mat'], dnrm.write);
        end
        spm_write_sn(V_mean,[proc{m} filesep 'affreg1norm.mat'], dnrm.write);
        m=m+1;
    end
end