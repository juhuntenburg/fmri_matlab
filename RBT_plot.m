clear all
close all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

files   =   { ...
     '/Users/francisca/Documents/Data/20180717_120447_CL_180718_visualstim_FFF_1_1/11/Processed'
};
scan_folder=[11];
search_name='rp_aiImage_';

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

for i=1:length(files)
    cd(files{i});
    A = fscanf(fopen([search_name num2str(scan_folder(i)) '_0001.txt'],'r'),'%f',[6 Inf])';
    h=figure(1);plot(A)
    set(h,'position',[178   192   895   513]);
    title('Rigid-body transformation parameters');
    xlabel('Volume');
    legend('Tx (mm)','Ty (mm)','Tz (mm)','Rx (rad)','Ry (rad)','Rz (rad)');
    set(legend,'position',[0.0173    0.7851    0.0821    0.1413]);
    saveas(gcf,'rbt_plot','tif');
end