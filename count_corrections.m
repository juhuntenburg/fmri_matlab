clear all

nslices=13;
study=[13:15,18:23];
folder = '/Users/francisca/Documents/Data/20171010_094318_FFF_med_3_loop_visual_freq_1_1/';

for i=1:length(study)
    k=1;
    for j=1:nslices
        clear g;
        try
            g=load([folder num2str(study(i)) '/Processed/to_correct' num2str(j)]);
        end
        if exist('g')
            m(i,k)=length(g.to_correct);
            k=k+1;
        end
    end
end
points=sum(m,2)/(340*nslices);