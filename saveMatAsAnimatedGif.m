%% SAVEMATASANIMATEDGIF 3D matrix as animated gif, third dimension is time
%
% USAGE
%   [true]=saveMatAsAnimatedGif(data,fileName,cmap,clim,resizeFactor,delayTime)
%
%   default values are
%   cmap = gray(256);
%   clim = [min(data(:)), max(data(:))]
%   resizeFactor = 1;
%   delayTime = 0.1;
%

% 2014-11-13, create by Till Huelnhagen

function [true]=saveMatAsAnimatedGif(data,fileName,cmap,clim,resizeFactor,delayTime)

    if ndims(data) ~= 3
        fprintf('\nERROR in saving animated gif: Data must have three dimensions\n\n');
        return
    end

    if nargin < 3
        cmap = gray(256);
        clim = [min(data(:)), max(data(:))];
        resizeFactor = 1;
        delayTime = 0.1;
    elseif nargin < 4
        clim = [min(data(:)), max(data(:))]
        resizeFactor = 1;
        delayTime = 0.1;
    elseif nargin < 5
        resizeFactor = 1;
        delayTime = 0.1;
    elseif nargin < 6
        delayTime = 0.1;
    end

    for i=1:size(data,3)
        % T2* map from registered filtered images
        t=data(:,:,i);
        %t(t<clim(1))=clim(1); t(t>clim(2))=clim(2); t=(t-clim(1))/clim(2)*size(cmap,1)/(clim(2)-clim(1));
        t(t<clim(1))=clim(1); t(t>clim(2))=clim(2); t=(t-clim(1))/(clim(2)-clim(1))*size(cmap,1);
        if resizeFactor ~= 1
            t=imresize(t,resizeFactor,'nearest');
        end
        t=uint8(t);
        if i == 1;
            imwrite(t,cmap,fileName,'gif','LoopCount',Inf,'DelayTime',delayTime);
        else
            imwrite(t,cmap,fileName,'gif','WriteMode','append','DelayTime',delayTime);
        end
    end

end