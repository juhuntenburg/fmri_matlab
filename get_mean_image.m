clear all
close all

%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%

name_of_image = 'mean_image';
% name_of_image = 'mean_image_all';

%%
%%%%%%%%%%%%%%%%%% Computations %%%%%%%%%%%%%%%%%%

F = spm_select(Inf,'image','select images');
nFsize = size(F,1);
for iscan = 1:nFsize
    f =deblank(F(iscan,:));
    V = spm_vol(f);
    Y = spm_read_vols(V);
    if exist('Ymean') ~= 1  Ymean = zeros(V.dim); end
    Ymean = Ymean+Y;
end
Ymean = Ymean/nFsize;   % divide YMean to get the mean
[dirname, xname, xext] = fileparts(V.fname);
fref = fullfile(dirname,[name_of_image xext]);
V.fname = fref;
spm_write_vol(V,Ymean);