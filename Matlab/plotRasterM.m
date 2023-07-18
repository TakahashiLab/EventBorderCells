
%verbose figureを出したくない時は0
%load event.mat
%plotRaster(kkOut{4,3},event,[1 2],'binWidth',100);
%%% Meanがevent=[1 2]で作った時と各event別で作った時で数字が合わない
%Nsta,Nen:解析したい周回数のstart,end
%SpkInIndex:TM中の活動

%event_num前後のspk出したい [264328 828239 3061386 7149731]みたいに入れる
%ensemble2{OK_ens{1},3}使う
function [histR,raster,HistR]=plotRasterM(spks,event_num,jitter)


spks=double(spks);
kHz=25;
binWidth=20;%msec


SpkInIndex=[];
l=size(event_num,2);
%     jitter=20;
    duration=jitter*1000*kHz*2;
    eventLength=floor(duration/(kHz*binWidth));

        raster=zeros(l,eventLength);
        
        for i=1:l
            if event_num(i)-jitter*1000*kHz>0
                SpkInTrial=spks-(event_num(i)-ceil(jitter*1000*kHz));
                SpkInIndex=[SpkInIndex find(SpkInTrial>=0 & SpkInTrial<=duration)];
                SpkInTrial=SpkInTrial(SpkInTrial>=0 & SpkInTrial<=duration);
                SpkSeq=floor(SpkInTrial/(kHz*binWidth));
                SpkSeq(SpkSeq==0)=[];
                SpkSeq(SpkSeq==eventLength+1)=[];
                raster(i,:)=full(sparse(1,SpkSeq,1,1,eventLength));

            end
        end
    histR=sum(raster,1)/l;
    for m=26:length(histR)-25
        HistR(1,m-25)=sum(histR(1,m-25:m+25));
    end

   

return;
%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Up,Do,l]=splitUpDo(Up,Do,splitNum,ThS,kHz)
% duration=(Do-Up)/(kHz*1000);
% splitP=find(abs(diff(duration))>ThS )+1;
% splitP=[1 splitP];
% %duration
% if length(splitP)==splitNum
%   r=splitP(splitNum):length(Do);
% else
%   r=splitP(splitNum):(splitP(splitNum+1)-1);
% end
% 
% Up=Up(r);
% Do=Do(r);
% l=length(Up);
% 
% return;