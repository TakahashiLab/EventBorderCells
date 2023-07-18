%event_num1(複数可)~event_num2のspk出したい event_Num1<event_Num2
% event_num2までの距離でsequence化
%ensemble2{OK_ens(1),3}使う
%TMtimes:トレッドミル回転中のevent_NumをPosT換算したもの
%
function [raster,HistR,DIFF]=plotRasterMMD(spks,event_Num1,event_Num2,Traj,PosT,TMtimes)
DIFF=[];DIFF2=[];
% BinWidthCm=2.5;
spks=double(spks);
kHz=25;
binWidth=20;
velocity=52.83;%cm/sec
LbinWidth=0.25;%400cm/max(LTraj)
distancePerbin=(velocity/(1000/binWidth))/LbinWidth;%TM中に20msあたりLTraj何bin分進んでいるか

[~,I1]=min(abs(event_Num1-PosT));
% [~,I2]=min(abs(event_Num2-PosT));

SpkInIndex=[];
l=size(event_Num1,2);
    duration=max(event_Num2-event_Num1);
    eventLength=max(floor(duration/(kHz*binWidth))+2);

        raster=zeros(l,eventLength);
%         figure;
        for i=1:l
            SpkInTrial=spks-(event_Num1(i));
            SpkInIndex=[SpkInIndex find(SpkInTrial>=0 & SpkInTrial<=duration)];
            SpkInTrial=SpkInTrial(SpkInTrial>=0 & SpkInTrial<=duration);
            SpkSeq=floor(SpkInTrial/(kHz*binWidth));
            SpkSeq(SpkSeq==0)=[];
            SpkSeq(SpkSeq==eventLength+1)=[];
            raster(i,:)=full(sparse(1,SpkSeq,1,1,eventLength));
            for j=I1(i)+2:I1(i)+eventLength
                if j<size(PosT,1)
                    if ~isempty(intersect(j,TMtimes))
                        DIFF(i,j-I1(i))=distancePerbin;
                        
                    else
                        DIFF(i,j-I1(i))=sqrt(diff(Traj(j-1:j,3))^2+diff(Traj(j-1:j,4))^2);
                    end
                end
            end
            DIFF2(i,:)=floor(cumsum(DIFF(i,:)));
%             DIFF2(i,:)=floor(cumsum(DIFF(i,:))/BinWidthCm);%結局glmするときに圧縮するのでbin幅は無視
%             hold on;plot(DIFF2(i,:))
           
        end
         for k=0:max(DIFF2(:))
                if find(DIFF2==k)
                    m=length(find(DIFF2==k));
                    histR(k+1)=sum(raster(find(DIFF2==k)))/m;
                end
         end
        for m=26:length(histR)-25
            HistR(1,m-25)=sum(histR(1,m-25:m+25));
        end

return;