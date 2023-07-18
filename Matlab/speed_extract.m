% event_num1(複数可)~event_num2のsp出したい event_Num1<event_Num2
% event_Num1=(DSpoint{i,j}-(size(DSzone4{i,j},2)/50/2)*25000);
% event_Num2=(DSpoint{i,j}+(size(DSzone4{i,j},2)/50/2)*25000);
% 似たようなspeed point,全然似てないspeed pointも出す...現在不使用
% all_sp:TM中の速度も加味したsession全体のspeed
%TMtimes:トレッドミル回転中のevent_NumをPosT換算したもの Traj Num
%DzoneTimes:TM回転時間+DS/DEzone Traj Num
function [SP,sp,otherSP,noSP,otherSPpoint,noSPpoint]=speed_extract(event_Num1,event_Num2,all_sp,Traj,PosT,TMtimes,DzoneTimes)
kHz=25;
binWidth=20;
velocity=52.83;%cm/sec
LbinWidth=0.25;%400cm/max(LTraj)
distancePerbin=(velocity/(1000/binWidth))/LbinWidth;%TM中に20msあたりLTraj何bin分進んでいるか

    [~,I1]=min(abs(event_Num1-PosT));%I1:Traj Num
    l=size(event_Num1,2);
    duration=max(event_Num2-event_Num1);%duration:event Num
    eventLength=max(floor(duration/(kHz*binWidth))+2);%eventLength:Traj Num
    for i=1:l
    for j=I1(i)+2:I1(i)+eventLength
        if j<size(PosT,1)
            if ~isempty(intersect(j,TMtimes))
                sp(i,j-I1(i))=distancePerbin;
            else
                sp(i,j-I1(i))=sqrt(diff(Traj(j-1:j,3))^2+diff(Traj(j-1:j,4))^2);
            end
        end
    end
    end
    for m=26:length(sp)-25
        SP(1,m-25)=sum(sp(1,m-25:m+25));
    end
    sp=sp(:,26:end-25);
    SP=SP/50;
    xx=1;index=[];otherSP=[];yy=1;noSP=[];otherSPpoint=[];noSPpoint=[];
    for k=(size(SP,2)/2+1):(length(all_sp)-size(SP,2)/2)
        Ind=(k-size(SP,2)/2+1):k+size(SP,2)/2;
        if ceil(Ind(1))==Ind(1) || floor(Ind(1))==Ind(1)
        else
            mieno
        end
        if isempty(intersect(DzoneTimes,Ind))
            if isempty(intersect(index,Ind))
                exSP=all_sp(Ind);
                if corr(exSP',SP')>0.5
                    otherSP(xx,:)=exSP;
                    xx=xx+1;
                    index=[index Ind];
                    otherSPpoint=[otherSPpoint ceil(k)];
    %                 figure;plot(SP)
    %                 hold on;plot(exSP)
                elseif corr(exSP',SP')<0.01
                    noSP(yy,:)=exSP;
                    yy=yy+1;
                    index=[index Ind];
                    noSPpoint=[noSPpoint ceil(k)];
                end
            end
        end
    end
return;
%先にやっとく
% dly1=[];dly2=[];DzoneTimes=[];  
%     [~,I3]=min(abs(event_Num{3}-PosT));
%     [~,I5]=min(abs(event_Num{5}-PosT));
%     [~,I4]=min(abs(event_Num{4}-PosT));
%     [~,I6]=min(abs(event_Num{6}-PosT));
% for i=1:size(event_Num{3},2)
%     if intersect(11:20,i)
%         h=3;
%         pre1=(250-min(DSzone4{3,h}));
%         pre2=(250-min(DSzone4{2,h}));
%         post1=(max(DEzone4{3,h})-1250);
%         post2=(max(DEzone4{2,h})-1250);
%     elseif intersect(21:30,i)
%         h=4;
%         pre1=(250-min(DSzone4{3,h}));
%         pre2=(250-min(DSzone4{2,h}));
%         post1=(max(DEzone4{3,h})-750);
%         post2=(max(DEzone4{2,h})-750);
%     elseif intersect(1:10,i)
%         h=2;
%         pre1=(250-min(DSzone4{3,h}));
%         pre2=(250-min(DSzone4{2,h}));
%         post1=(max(DEzone4{3,h})-750);
%         post2=(max(DEzone4{2,h})-750);
%     end
%     if isempty(I3:I5) || isempty(I4:I6)
%     mieno
%     else
%     dly1=[dly1 floor(I3(i)-pre1):floor(I5(i)+post1)];
%     dly2=[dly2 floor(I4(i)-pre2):floor(I6(i)+post2)];
%     end
% end
% DzoneTimes=[];
% DzoneTimes=[dly1 dly2];
% for i=length(DzoneTimes):-1:1
% if DzoneTimes(i)>length(Traj)
% DzoneTimes(:,i)=[];
% end
% end
% 
% [all_sp]=speed_extract(1,length(event),speed,Traj,PosT,TMtimes,DzoneTimes);
        

