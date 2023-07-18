% cellの順でpeak index(NI)とfield width(Fwidth)
% countPlaceFieldMのFwidthのサイズ確認しとく
function [N_seq,Norder,NI,Fwidth]=Nseq(seq,cell)
    if nargin==2
        seq=seq(cell,:);
    else    
    end
    [~,NI]=max(seq,[],2);
    [~,Norder]=sort(NI);
    x=1;
    for j=1:size(seq,1)
%         if max(seq(j,:))>1
        N_seq(x,:)=seq(j,:)/max(seq(j,:));
        MM=2;
        [out,xy]=countPlaceFieldsM(seq(j,:),1,5,MM);
        if isempty(out)
            Fwidth(j)=0;
        else
            for k=1:size(xy,2)
                if intersect(xy{k}(:,1),NI(j))
                    Fwidth(j)=out(k);
                else
                end
            end
        end
        x=x+1;
%         end
    end

return;