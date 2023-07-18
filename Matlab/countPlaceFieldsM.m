% MM=1のときpeakRate/5,MM=2のときpeakRate*3/5かつline29→line30~33, line36 off

function [out,B,peak_Point,b]=countPlaceFieldsM(rate_map,binside,ThF,MM)
b=[];peak_Point=[];B=[];o=1;

if nargin==2
    ThF=50;    
end
out=[];
xy=[];
verbose=0;

rm = rate_map;
if size(rate_map,1)==1
    binside2=binside;
else
    binside2=binside^2;
end

peakRate=max(rm(:));
if peakRate>0.7%mieno

    if MM==1
        sp_thr = peakRate/5;
    elseif MM==2
        sp_thr = peakRate*3/5;
    end
    field = (rm>=sp_thr);

    a = regionprops(field,{'Area' 'PixelList'});
    

    b = [];

    c=1;d=1;
    for i = 1:length(a)
        if ((a(i).Area * binside2) > ThF)
            b{c} = a(i).PixelList;
            [M,I]=max(rate_map(b{c}(:,1)));
            if MM==1
            elseif MM==2
                M=[];
                for j=1:size(b{c},1)
                M(j)=rate_map(b{c}(j,2),b{c}(j,1));
                end
            end
            if any(M>0.7)%mieno
                B{d}=b{c};
                peak_Point(d)=I+b{c}(1,1)-1;
                d=d+1;
                out=[out a(i).Area * binside2];
            end
            c=c+1;
        end
    end

    if verbose
        imagesc(rate_map);colormap jet;
        for k=1:(c-1)
            xy=b{k};
            hold on;
            plot(xy(:,1),xy(:,2),'.');
        end
    end
end

% for i=1:size(b,2)
%     [~,I]=max(rate_map(b{i}(:,1)));
%     peak_Point(i)=I+b{i}(1,1)-1;
% end

% if nargout==3
%     B=[];
%     BB=[];
%     for i=1:length(b)
%         BB=median(b{i});
%         B=[B BB(1)];
%     end
% end

return;