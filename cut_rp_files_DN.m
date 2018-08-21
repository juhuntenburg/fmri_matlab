
folder      =   '/Users/DanielNunes/Documents/MATLAB/A_DATA_Daniel/dfMRI/2ndRound_CryoCoil/Anisotropy_Average_all/dfMRI_261016';
exps        =   [13:16];
vols_skip   =   20;   



for i=1:length(exps)
    t_dir           =   fullfile(folder,num2str(exps(i)),'Processed');
    file            =   spm_select('FPList',t_dir,'^rp.*\d\d\d\d.txt');  
    [rr name ext]   =   fileparts(file);   
    fid             =   fopen(file,'r');
    dat             =   importdata(file);    
    out_dat         =   dat((vols_skip+1):end,:);
    f_out           =   fopen(fullfile(t_dir,[name '_mod.txt']),'w');
    s               =   size(out_dat,1);
    i               =   1;
    while (i<=s)
        fprintf(f_out,'%s\r\n',num2str(out_dat(i,:)));
        i           =   i+1;
    end
    fclose(f_out);
end