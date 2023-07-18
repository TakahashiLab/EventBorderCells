% ALL_timeseqC2のDS/DE zone内のpeak firing timingを中心にしてspをカラーマッピングする
% sorted by sp peak timing でsp sequence作る
% cell_Num:DS/DE cellのNum, seq:ALL_DS/DEseqC2, sp:ALL_sp 20ms/1bin 150cells×sp(30s/40s=1499/1999bin)
% zone:DS/DE zone
% corr_Flag=1のとき　corr histogram出す　speed vs seq(seqを1000回shffledしたもの)
% corr_Flag=0のとき　seq_oriに[]いれとく
% speedのimagescの上にtime fieldを半透明で乗せる

% 20230712 index＝DS/E pointを細胞数分入れて、 for l=cell_Numの中でindex(l)で回せるようにする
function [SP,N_order,Field_list,TESTS,AA]=sp_seq(cell_Num,sp,seq,seq_ori,point,corr_Flag,Title1,Title2,Title3)
if nargin<7
    corr_Flag=0;
    Title1=1;
    Title2=1;
    Title3=1;
end
% index=floor(median(zone));20230712

        Field_list=zeros(size(cell_Num,2),500);
        x=1;SP=[];
%         FH1=figure;
        xx=1;TESTS=[];AA=[];
        for l=cell_Num
        index=floor(median(point(xx)));
        xx=xx+1;
        [~,I]=max(seq(l,:));
        SPI=500+index+I-250;
        SP(x,:)=sp(l,SPI-249:SPI+250);
%         SP(x,:)=sp{i,j}(:,SPI-149:SPI+150);
        
        test=NaN(1,500);
        if I-249<0
            test((251-I):500)=seq(l,1:(250+I));
        elseif I+250>500
            test(1:(750-I))=seq(l,(I-249):500);
        else
            test=seq(l,:);
        end
        
        Test=test;
        Test(isnan(test))=0;
        ThF=50;MM=1;
%         save 20230711.mat -append Test
        [~,xy]=countPlaceFieldsM(Test,1,ThF,MM);
        if isempty(xy)
            TT=seq(l,:);
%             save toriaezu.mat -append Test TT
        else
        for i=1:size(xy,2)
            if intersect(xy{i}(:,1),225:275)
                Field_list(x,xy{i}(:,1))=1;
%                 Field_list(x,xy{i}(:,1))=test(1,xy{i}(:,1));
                break
            else

            end
        end
        end

        test=test/max(test);
        

        seq1=seq(l,:);
        seq_ori1=seq_ori(l,:);
        sp1=SP(x,:);
        shN=1000;
        [testS,A]=corr_histo(seq1,seq_ori1,sp1,shN,l,Title1,Title2,Title3,corr_Flag);
        TESTS=[TESTS testS];
        AA=[AA A];
        x=x+1;
        end
        
        

        [~,N_order]=Nseq(SP);
        SP=SP*12.255;%bin/ms→cm/sec
        FH2=figure;imagesc(SP(N_order,:));
        colormap('pink');
        cc=colorbar;
        cc.Label.String = 'running speed (cm/s)';
        c=[0 0.4470 0.7410];
        y=1;
        for j=N_order'
            A=find(Field_list(j,:)==1);
            if isempty(A)
                TT=seq(l,:);
            save toriaezu.mat -append Test TT
            else
            hold on;patch([A(1) A(end) A(end) A(1)],[y-.5 y-.5 y+.5 y+.5],c,'FaceAlpha',.5,'EdgeColor','none')% ほんとは発火率に応じて濃淡を変えたい
            y=y+1;
            end
        end
        
return;


%%%%%%%%%%%%%%%%%%%%%%%%%
function [testS,A]=corr_histo(seq1,seq_ori1,sp1,shN,cell_Num,Title1,Title2,Title3,corr_Flag)
    A=corr(seq1',sp1');
    testS=[];
    for i=1:shN
        seqS=seq_ori1(:,randperm(size(seq_ori1,2)));
        seqSS=[];
        for j=26:length(seqS)-25
            seqSS(1,j-25)=sum(seqS(1,j-25:j+25));
        end 
        B=corr(seqSS',sp1');
        testS=[testS B];
    end
    if corr_Flag==1
        C=sort(testS);
        B=sort([testS A]);
        p=(find(B==A)/shN)*100;
        Pval=1000*.025;
        figure;
        M=histcounts(testS,size(seq1,2)/20);
        c=[0.9290 0.6940 0.1250];
        XX=[C(Pval) C(shN-Pval) C(shN-Pval) C(Pval)];
        YY=[0 0 max(M*1.2) max(M*1.2)];
        hold on;patch(XX,YY,c,'FaceAlpha',.5,'EdgeColor','none')
        hold on;xline(C(Pval),'--','2.5%','LabelOrientation','horizontal','LabelHorizontalAlignment','center')
        hold on;xline(C(shN-Pval),'--','97.5%','LabelOrientation','horizontal','LabelHorizontalAlignment','center')
        hold on;histogram(testS,size(seq1,2)/20)
        hold on;xline(A,'r','LineWidth',2)
    %     hold on;plot(A,linspace(realmin,max(M*1.2),500),'r','Marker','|','LineWidth',2)
        xlabel('correlation between speed and activity')
        ylabel('fraction of shuffles')
        title2=sprintf('p=%d',p);
        title(title2)
        ylim([0 max(M*1.2)])
        axis square
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        title1=[];
        if C(Pval)>A
            title1=sprintf('corr(sp vs FR) cell%d under5per',cell_Num);
        elseif  C(shN-Pval)<A
            title1=sprintf('corr(sp vs FR) cell%d over95per',cell_Num);
        else
            title1=sprintf('corr(sp vs FR) cell%d',cell_Num); 
        end
        Title=sprintf('%s %s %s %s',Title1,Title2,Title3,title1);
        printFig(gcf,Title);
    end
return;


% 先にやっとく
% 要修正
% 各試行ごとにspeed sequenceとその上にtime fieldをのっけた図を出す
% 細胞活動が加速行動に関係しているのかイベントボーダーに関係しているのか
% 各試行ごとに速度変化と細胞活動のピアソン相関と細胞活動をシャッフリングした時の相関を計算する
% 観察データがシャッフリングデータ内でどの順位にいるか＝パーセンタイル相関
% シャッフリングデータの上位もしくは下位５％以下に入っていれば異なる分布とみなしてcell Num控えとく
% シャッフリングよりも観察データが上回っているか下回っているかのグループ分けはしとく
% それらの細胞がDS/DE zoneに多いかどうか
% 各個体毎にduration±(jitter+5s=10s)のsp出す 20ms/1bin 6個体×sp(20s/30s)が2TM×3stage
% [all_sp]=speed_extract(1,length(event),speed,Traj,PosT,TMtimes,DzoneTimes);
% event_Num1=[];event_Num2=[];
% for i=1:3
%     for j=1:2
%         if j==1
%             event_Num1=event_Num{3}(((i-1)*10+1):(i*10))-25000*15;
%             event_Num2=event_Num{5}(((i-1)*10+1):(i*10))+25000*15;
%         elseif j==2
%             event_Num1=event_Num{4}(((i-1)*10+1):(i*10))-25000*15;
%             event_Num2=event_Num{6}(((i-1)*10+1):(i*10))+25000*15;
%         end
%         [All_sp{j,i}]=speed_extract(event_Num1,event_Num2,all_sp,Traj,PosT,TMtimes,DzoneTimes);
%         figure;plot(All_sp{j,i})
%     end
% end
% % 全部合体　ALL_sp 20ms/1bin 150cells×sp(30s/40s=1499/1999bin)
% % X:loadの時できたやつ
% ALL_sp=cell(2,3);
% for i=1:6
% addpath(PathName{i})
% load 20220318.mat All_sp
% for k=1:2
%     for l=1:3
%         for j=1:X{i}
%             ALL_sp{k,l}=vertcat(ALL_sp{k,l},All_sp{k,l});
%         end
%     end
% end
% clear All_sp
% end