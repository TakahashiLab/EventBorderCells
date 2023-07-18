%%% 221118 個体O抜きで解析
%%% 220511 series F-7 InfoSparce, countFieldにMM分岐追加
%%%各々 
%最初の1-3周の比較 time_seqF1:1-3周　time_seqF2:4-10周 time_seqF3:4-6周　time_seqF4:7-10周
%fig title : time_seqF1
time_seqF1=cell(3,4);%time_seq{2,:}のときevent2
time_seqF1_ori=cell(3,4);
HistR=cell(1,2);
for h=2:4
Nsta=10*(h-1)+1;
Nen=10*(h-1)+4;%
for i=1:size(OK_ens,2)
[~,histR]=plotRasterMM2(ensemble2{OK_ens(i),3},event,[1 2],PORK3,Nsta,Nen,'verbose',0,'jitter',5.5,'binwidth',20);
for j=1:2
for m=26:length(histR{j})-25
    HistR{j}(1,m-25)=sum(histR{j}(1,m-25:m+25));
end
end
time_seqF1{2,h}(i,:)=HistR{2};
time_seqF1{3,h}(i,:)=HistR{1};
time_seqF1_ori{2,h}(i,:)=histR{2};
time_seqF1_ori{3,h}(i,:)=histR{1};
HistR=[];HistR{2}=[];
end
end
%fig title : time_seqD2
time_seqF2=cell(3,4);%time_seq{2,:}のときevent2
time_seqF2_ori=cell(3,4);
HistR=[];HistR{2}=[];
for h=2:4
Nsta=10*(h-1)+4;
Nen=10*h;%
for i=1:size(OK_ens,2)
[~,histR]=plotRasterMM2(ensemble2{OK_ens(i),3},event,[1 2],PORK3,Nsta,Nen,'verbose',0,'jitter',5.5,'binwidth',20);
for j=1:2
for m=26:length(histR{j})-25
    HistR{j}(1,m-25)=sum(histR{j}(1,m-25:m+25));
end
end
time_seqF2{2,h}(i,:)=HistR{2};
time_seqF2{3,h}(i,:)=HistR{1};
time_seqF2_ori{2,h}(i,:)=histR{2};
time_seqF2_ori{3,h}(i,:)=histR{1};
HistR=[];HistR{2}=[];
end
end
%fig title : time_seqD3
time_seqF3=cell(3,4);%time_seq{2,:}のときevent2
time_seqF3_ori=cell(3,4);
HistR=[];HistR{2}=[];
for h=2:4
Nsta=10*(h-1)+4;
Nen=10*(h-1)+7;%
for i=1:size(OK_ens,2)
[~,histR]=plotRasterMM2(ensemble2{OK_ens(i),3},event,[1 2],PORK3,Nsta,Nen,'verbose',0,'jitter',5.5,'binwidth',20);
for j=1:2
for m=26:length(histR{j})-25
    HistR{j}(1,m-25)=sum(histR{j}(1,m-25:m+25));
end
end
time_seqF3{2,h}(i,:)=HistR{2};
time_seqF3{3,h}(i,:)=HistR{1};
time_seqF3_ori{2,h}(i,:)=histR{2};
time_seqF3_ori{3,h}(i,:)=histR{1};
HistR=[];HistR{2}=[];
end
end
%fig title : time_seqD4
time_seqF4=cell(3,4);%time_seq{2,:}のときevent2
time_seqF4_ori=cell(3,4);
HistR=[];HistR{2}=[];
for h=2:4
Nsta=10*(h-1)+7;
Nen=10*h;%
for i=1:size(OK_ens,2)
[~,histR]=plotRasterMM2(ensemble2{OK_ens(i),3},event,[1 2],PORK3,Nsta,Nen,'verbose',0,'jitter',5.5,'binwidth',20);
for j=1:2
for m=26:length(histR{j})-25
    HistR{j}(1,m-25)=sum(histR{j}(1,m-25:m+25));
end
end
time_seqF4{2,h}(i,:)=HistR{2};
time_seqF4{3,h}(i,:)=HistR{1};
time_seqF4_ori{2,h}(i,:)=histR{2};
time_seqF4_ori{3,h}(i,:)=histR{1};
HistR=[];HistR{2}=[];
end
end

%%%各々
% DlyCell:duration±jitter5s内で活動(InfoSparse>shuffle95%,peakFR>1Hz,1s<Fsize<10sのfieldが一個以上ある)していた細胞
% countPlaceFieldsMのline 16 peakRate/5を確認　MM＝１
% TimeCell:DlyCellの条件にwithin corr>.3とGLMを加えたもの
DlyCellF=cell(3,4);
for h=1:3
for i=1:4
    if ~isempty(time_seqF2{h,i})
        sequence3=time_seqF2_ori{h,i};
        for j=1:size(OK_ens,2)
        I1=calcInfoM(time_seqF2{h,i}(j,:));
        testS1=[];
        for g=1:1000
            seqS2=sequence3(j,randperm(size(sequence3,2)));
            seqS3=[];
            for k=26:length(seqS2)-25
                seqS3(1,k-25)=sum(seqS2(1,k-25:k+25));
            end
            It=calcInfoM(seqS3);
            testS1=[testS1 It];
        end
        A1=sort(testS1);
        MM=1;ThF=50;
        [out]=countPlaceFieldsM(time_seqF2{h,i}(j,:),1,ThF,MM);
        T_index=[];
        if ~isempty(out) && I1<A1(50) && any(out<500)
            DlyCellF{h,i}=[DlyCellF{h,i} j];
        end
        end
    end
end
end
save 20220318.mat -append time_seqF1 time_seqF2 time_seqF3 time_seqF4 time_seqF1_ori time_seqF2_ori time_seqF3_ori time_seqF4_ori DlyCellF


%%%全部
PathName{1}=('G:\OneDrive - 同志社大学\20220307\2\2_201208');
PathName{2}=('G:\OneDrive - 同志社大学\20220307\5');
PathName{3}=('G:\OneDrive - 同志社大学\20220307\9\210105');
PathName{4}=('G:\OneDrive - 同志社大学\20220307\P');
PathName{5}=('G:\OneDrive - 同志社大学\20220307\S');
X=[];
for i=1:5
addpath(PathName{i})
load 20220318.mat OK_ens
X{i}=size(OK_ens,2);
end
x=0;
ALL_timeseqF1=cell(3,4);ALL_timeseqF2=cell(3,4);ALL_timeseqF3=cell(3,4);ALL_timeseqF4=cell(3,4);
ALL_timeseqF1_ori=cell(3,4);ALL_timeseqF2_ori=cell(3,4);ALL_timeseqF3_ori=cell(3,4);ALL_timeseqF4_ori=cell(3,4);
ALL_DlyCellF=cell(3,4);
for i=1:5
addpath(PathName{i})
load 20220318.mat time_seqF1 time_seqF2 time_seqF3 time_seqF4 time_seqF1_ori time_seqF2_ori time_seqF3_ori time_seqF4_ori DlyCellF
if i==1
x=0;
else
x=x+X{i-1};
end
for j=1:4
    for k=1:3
        if ~isempty(time_seqF1{k,j})
            ALL_timeseqF1{k,j}=vertcat(ALL_timeseqF1{k,j},time_seqF1{k,j});
            ALL_timeseqF1_ori{k,j}=vertcat(ALL_timeseqF1_ori{k,j},time_seqF1_ori{k,j});
        end
        if ~isempty(time_seqF2{k,j})
            ALL_timeseqF2{k,j}=vertcat(ALL_timeseqF2{k,j},time_seqF2{k,j});
            ALL_timeseqF2_ori{k,j}=vertcat(ALL_timeseqF2_ori{k,j},time_seqF2_ori{k,j});
        end
        if ~isempty(time_seqF3{k,j})
            ALL_timeseqF3{k,j}=vertcat(ALL_timeseqF3{k,j},time_seqF3{k,j});
            ALL_timeseqF3_ori{k,j}=vertcat(ALL_timeseqF3_ori{k,j},time_seqF3_ori{k,j});
        end
        if ~isempty(time_seqF4{k,j})
            ALL_timeseqF4{k,j}=vertcat(ALL_timeseqF4{k,j},time_seqF4{k,j});
            ALL_timeseqF4_ori{k,j}=vertcat(ALL_timeseqF4_ori{k,j},time_seqF4_ori{k,j});
        end
        if ~isempty(DlyCellF{k,j})
            ALL_DlyCellF{k,j}=[ALL_DlyCellF{k,j} DlyCellF{k,j}+x];
        end
    end
end
clear time_seqF1 time_seqF2 time_seqF3 time_seqF4 time_seqF1_ori time_seqF2_ori time_seqF3_ori time_seqF4_ori DlyCellF
end
% compress ver.
for g=2:3
for h=1:89
for i=1:499
ALL_timeseqF1{g+2,3}(h,i)=(ALL_timeseqF1{g,3}(h,250+2*i-1)+ALL_timeseqF1{g,3}(h,250+2*i))/2;
ALL_timeseqF1_ori{g+2,3}(h,i)=(ALL_timeseqF1_ori{g,3}(h,275+2*i-1)+ALL_timeseqF1_ori{g,3}(h,275+2*i))/2;
ALL_timeseqF2{g+2,3}(h,i)=(ALL_timeseqF2{g,3}(h,250+2*i-1)+ALL_timeseqF2{g,3}(h,250+2*i))/2;
ALL_timeseqF2_ori{g+2,3}(h,i)=(ALL_timeseqF2_ori{g,3}(h,275+2*i-1)+ALL_timeseqF2_ori{g,3}(h,275+2*i))/2;
ALL_timeseqF3{g+2,3}(h,i)=(ALL_timeseqF3{g,3}(h,250+2*i-1)+ALL_timeseqF3{g,3}(h,250+2*i))/2;
ALL_timeseqF3_ori{g+2,3}(h,i)=(ALL_timeseqF3_ori{g,3}(h,275+2*i-1)+ALL_timeseqF3_ori{g,3}(h,275+2*i))/2;
ALL_timeseqF4{g+2,3}(h,i)=(ALL_timeseqF4{g,3}(h,250+2*i-1)+ALL_timeseqF4{g,3}(h,250+2*i))/2;
ALL_timeseqF4_ori{g+2,3}(h,i)=(ALL_timeseqF4_ori{g,3}(h,275+2*i-1)+ALL_timeseqF4_ori{g,3}(h,275+2*i))/2;
end
end
ALL_timeseqF1{g+2,3}=[ALL_timeseqF1{g,3}(:,1:250) ALL_timeseqF1{g+2,3} ALL_timeseqF1{g,3}(:,end-249:end)];
ALL_timeseqF1_ori{g+2,3}=[ALL_timeseqF1_ori{g,3}(:,1:275) ALL_timeseqF1_ori{g+2,3} ALL_timeseqF1_ori{g,3}(:,end-274:end)];
ALL_timeseqF2{g+2,3}=[ALL_timeseqF2{g,3}(:,1:250) ALL_timeseqF2{g+2,3} ALL_timeseqF2{g,3}(:,end-249:end)];
ALL_timeseqF2_ori{g+2,3}=[ALL_timeseqF2_ori{g,3}(:,1:275) ALL_timeseqF2_ori{g+2,3} ALL_timeseqF2_ori{g,3}(:,end-274:end)];
ALL_timeseqF3{g+2,3}=[ALL_timeseqF3{g,3}(:,1:250) ALL_timeseqF3{g+2,3} ALL_timeseqF3{g,3}(:,end-249:end)];
ALL_timeseqF3_ori{g+2,3}=[ALL_timeseqF3_ori{g,3}(:,1:275) ALL_timeseqF3_ori{g+2,3} ALL_timeseqF3_ori{g,3}(:,end-274:end)];
ALL_timeseqF4{g+2,3}=[ALL_timeseqF4{g,3}(:,1:250) ALL_timeseqF4{g+2,3} ALL_timeseqF4{g,3}(:,end-249:end)];
ALL_timeseqF4_ori{g+2,3}=[ALL_timeseqF4_ori{g,3}(:,1:275) ALL_timeseqF4_ori{g+2,3} ALL_timeseqF4_ori{g,3}(:,end-274:end)];
end
%%%全部
% time seq
% vsStages  fig file title : Ntime_seq(sorted by stage 2 trial)
% DEzoneなしver.
for i=2:3 
    for j=2:4
        if j==3
            XX=1250;
        else
            XX=750;
        end
        [~,N_order]=Nseq(ALL_timeseqF2{i,2}(:,1:750),ALL_DlyCellF{i,2});
        [N_seq]=Nseq(ALL_timeseqF2{i,j}(:,1:XX),ALL_DlyCellF{i,2});
        figure;imagesc(N_seq(N_order,:))
        title3=sprintf('stage%d TM%d',j,i-1);
        title3_1=sprintf('Stage%d Treadmill%d',j,i-1);
        colorbar;
        title(title3_1)
        xlabel('Elapsed time (s)')
        ylabel({'Number of cells';'(sorted by stage 2 trial)'})
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        cc= [1.0000    0.8500    0.8500];
        hold on;xline(250,'Color',cc,'LineWidth',5)
        xticks([50:100:length(N_seq)])
        printFig(gcf,title3);
    end
end
% DEzoneありver.
for i=2:3
    [~,N_order]=Nseq(ALL_timeseqF2{i,2},ALL_DlyCellF{i,2});
    for j=2:4
    [N_seq]=Nseq(ALL_timeseqF2{i,j},ALL_DlyCellF{i,2});
    figure;imagesc(N_seq(N_order,:))
    title3=sprintf('stage%d TM%d',j,i-1);
    title3_1=sprintf('Stage%d Treadmill%d',j,i-1);
    colorbar;
    title(title3_1)
    xlabel('Elapsed time (s)')
    ylabel({'Number of cells';'(sorted by stage 2 trial)'})
    set(gca, 'FontName', 'Times New Roman','FontSize',16)
    cc= [1.0000    0.8500    0.8500];
    hold on;xline(250,'Color',cc,'LineWidth',5)
    cc= [0.8500 0.3250 0.0980];
    hold on;xline(length(N_seq)-250,'Color',cc,'LineWidth',5)
    xticks([50:100:length(N_seq)])
    printFig(gcf,title3);
    end
end
% vsTMs  fig file title : Ntime_seq(sorted by TM1 trial)
% DEzoneなしver.
for j=2:4
    if j==3
        XX=1250;
    else
        XX=750;
    end
    for i=2:3 
        [~,N_order]=Nseq(ALL_timeseqF2{2,j}(:,1:XX),ALL_DlyCellF{2,j});
        [N_seq]=Nseq(ALL_timeseqF2{i,j}(:,1:XX),ALL_DlyCellF{2,j});
        figure;imagesc(N_seq(N_order,:))
        title3=sprintf('stage%d TM%d',j,i-1);
        title3_1=sprintf('Stage%d Treadmill%d',j,i-1);
        colorbar;
        title(title3_1)
        xlabel('Elapsed time (s)')
        ylabel({'Number of cells';'(sorted by treadmill 1 trial)'})
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        cc= [1.0000    0.8500    0.8500];
        hold on;xline(250,'Color',cc,'LineWidth',5)
        xticks([50:100:length(N_seq)])
        printFig(gcf,title3);
    end
end
% DEzoneありver.
for j=2:4
    for i=2:3
    [~,N_order]=Nseq(ALL_timeseqF2{2,j},ALL_DlyCellF{2,j});
    [N_seq]=Nseq(ALL_timeseqF2{i,j},ALL_DlyCellF{2,j});
    figure;imagesc(N_seq(N_order,:))
    title3=sprintf('stage%d TM%d',j,i-1);
    title3_1=sprintf('Stage%d Treadmill%d',j,i-1);
    colorbar;
    title(title3_1)
    xlabel('Elapsed time (s)')
    ylabel({'Number of cells';'(sorted by treadmill 1 trial)'})
    set(gca, 'FontName', 'Times New Roman','FontSize',16)
    cc= [1.0000    0.8500    0.8500];
    hold on;xline(250,'Color',cc,'LineWidth',5)
    cc= [0.8500 0.3250 0.0980];
    hold on;xline(length(N_seq)-250,'Color',cc,'LineWidth',5)
    xticks([50:100:length(N_seq)])
    printFig(gcf,title3);
    end
end
%%%全部
% seq比較　dot plot ver.
% fig file title : dot plot(VS stages)
% countPlaceFieldsMのline 16 peakRate*1/5を確認
% vsStages
ALL_dot_VSstage2=cell(3,3);
for i=2:3
    for j=1:3 
        if j==1
            j1=2;j2=3;
            title3=sprintf('stage2vs3 TM%d',i-1);
            title3_1=sprintf('Stage2vs3 Treadmill %d',i-1);
        elseif j==2
            j1=4;j2=3;
            title3=sprintf('stage3vs4 TM%d',i-1);
            title3_1=sprintf('Stage3vs4 Treadmill %d',i-1);
        elseif j==3
            j1=2;j2=4;
            title3=sprintf('stage2vs4 TM%d',i-1);
            title3_1=sprintf('Stage2vs4 Treadmill %d',i-1);
        end
        X1=length(ALL_timeseqF2{i,j1});
        X2=length(ALL_timeseqF2{i,j2});
        jitter=250;
        XXX=seq_dot(ALL_timeseqF2{i,j1}(:,(251-jitter):X1-(250-jitter)),ALL_timeseqF2{i,j2}(:,(251-jitter):X2-(250-jitter)),intersect(ALL_DlyCellF{i,j1},ALL_DlyCellF{i,j2}),0,50,jitter);
        ALL_dot_VSstage2{i,j}=XXX;
        FH=gca;
        title(title3_1)
%         close
%         printFig(gcf,title3);
        % density color map ver.
        % 1sec/bin SD=0.8
        D=[];
        for k1=1:ceil((X1-(250-jitter)*2)/50)
            for k2=1:ceil((X2-(250-jitter)*2)/50)
                Y=[];
                Y=find((50*(k1-1)+1)<XXX(:,1) & XXX(:,1)<50*k1 & (50*(k2-1)+1)<XXX(:,2) & XXX(:,2)<50*k2);
                D(k1,k2)=size(Y,1)/size(XXX,1);
            end
        end
        D=D*100;
        SD=0.8;
        D=imgaussfilt(D,SD);
        FH=figure;imagesc(D');colorbar
        J=jitter/50;
        cc= [1.0000    0.8500    0.8500];
        hold on;xline(J,'Color',cc,'LineWidth',5)
        hold on;yline(J,'Color',cc,'LineWidth',5)
        cc= [0.8500 0.3250 0.0980];
        hold on;xline(size(D,1)-J,'Color',cc,'LineWidth',5)
        hold on;yline(size(D,2)-J,'Color',cc,'LineWidth',5)
        axis equal
        xlim([1 size(D,1)])
        ylim([1 size(D,2)])
        axis xy
        xticks([0:J:size(D,1)])
        yticks([0:J:size(D,2)])
%         xticks([0.5:2:size(D,1)])
%         yticks([0.5:2:size(D,2)])
        xlabel('Elapsed time (s)')
        ylabel('Elapsed time (s)')
        title(title3_1)
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        title4=sprintf('%s D',title3);
%         printFig(FH,title4);
    end
end
% fig file title : dot plot(VS TMs)
% vsTMs
ALL_dot_VSTM2=cell(1,4);
for j=2:4
    XX=length(ALL_timeseqF2{2,j});
    title3=sprintf('TM1vs2 stage%d',j);
    title3_1=sprintf('Treadmill 1vs2 stage%d',j);
    jitter=250;
    XXX=seq_dot(ALL_timeseqF2{2,j}(:,(251-jitter):XX-(250-jitter)),ALL_timeseqF2{3,j}(:,(251-jitter):XX-(250-jitter)),intersect(ALL_DlyCellF{2,j},ALL_DlyCellF{3,j}),0,50,jitter);
    ALL_dot_VSTM2{1,j}=XXX;
    FH=gca;
    title(title3_1)
%     close
    printFig(gcf,title3);
    % density color map ver.
    % 1sec/bin SD=0.8
    D=[];
    for k1=1:ceil((XX-(250-jitter)*2)/50)
        for k2=1:ceil((XX-(250-jitter)*2)/50)
            Y=[];
            Y=find((50*(k1-1)+1)<XXX(:,1) & XXX(:,1)<50*k1 & (50*(k2-1)+1)<XXX(:,2) & XXX(:,2)<50*k2);
            D(k1,k2)=size(Y,1)/size(XXX,1);
        end
    end
    D=D*100;
    SD=0.8;
    D=imgaussfilt(D,SD);
    FH=figure;imagesc(D');colorbar
    J=jitter/50;
    cc= [1.0000    0.8500    0.8500];
    hold on;xline(J,'Color',cc,'LineWidth',5)
    hold on;yline(J,'Color',cc,'LineWidth',5)
    cc= [0.8500 0.3250 0.0980];
    hold on;xline(size(D,1)-J,'Color',cc,'LineWidth',5)
    hold on;yline(size(D,2)-J,'Color',cc,'LineWidth',5)
    axis equal
    xlim([1 size(D,1)])
    ylim([1 size(D,2)])
    axis xy
    xticks([0:J:size(D,1)])
    yticks([0:J:size(D,2)])
    title4=sprintf('%s D',title3);
    xlabel('Elapsed time (s)')
    ylabel('Elapsed time (s)')
    title(title3_1)
    set(gca, 'FontName', 'Times New Roman','FontSize',16)
    printFig(FH,title4);
end

%%%全部
%fig title : duration Hist(allTMs&stages)
% accumulation of time fields
% time seqをnormalizedして
% ここのcountPlaceFieldsM変更有　line16 1/5→3/5,line29→line30~33, line36 off
DSzone7=cell(5,4);
DEzone7=cell(5,4);
DSpoint7=cell(5,4);
DEpoint7=cell(5,4);
noDzone7=cell(5,4);
AAA=[];BBB1=[];BBB2=[];
options = statset('MaxIter',1000);
z=4;M=[];
for g=2:3
for h=2:4
test=[];
AIC = zeros(1,z);
GMModels = cell(1,z);
S=250;
if h==3
    x=150;
    E=1350;
else
    x=100;
    E=850;
end
test=[];
if ~isempty(ALL_timeseqF2{g,h})
for i=ALL_DlyCellF{g,h}
out=[];xy=[];
MM=1;ThF=50;%20220608変更10→50
[out,xy,peak_Point]=countPlaceFieldsM(ALL_timeseqF2{g,h}(i,:),1,ThF,MM);
T=size(ALL_timeseqF2{g,h}(i,:),2);
if ~isempty(out)
for j=1:size(xy,2)
I=zeros(1,T);
I(1,xy{j}(:,1))=ALL_timeseqF2{g,h}(i,xy{j}(:,1))/max(ALL_timeseqF2{g,h}(i,xy{j}(:,1)));
test=[test xy{j}(:,1)'];
end
end
end
for k=1:z
GMModels{k} = fitgmdist(test',k,'Options',options,'CovarianceType','diagonal');
AIC(k)= GMModels{k}.AIC;
end
[minAIC,numComponents] = min(AIC);
y=numComponents;
GMModel = fitgmdist(test',y);
xgrid = linspace(min(test),max(test),x)';
[~,A]=min(abs(S-GMModel.mu));
[~,B]=min(abs(E-GMModel.mu));
DSpoint7{g,h}=GMModel.mu(A);
DEpoint7{g,h}=GMModel.mu(B);
DSzone7{g,h}=max(26,floor(GMModel.mu(A)-sqrt(GMModel.Sigma(A)*3))):min(length(ALL_timeseqF2{g,h})-26,ceil(GMModel.mu(A)+sqrt(GMModel.Sigma(A)*3)));
DEzone7{g,h}=max(26,floor(GMModel.mu(B)-sqrt(GMModel.Sigma(B)*3))):min(length(ALL_timeseqF2{g,h})-26,ceil(GMModel.mu(B)+sqrt(GMModel.Sigma(B)*3)));
figure
yyaxis left
hold on;histogram(test)
yyaxis right
hold on; plot(xgrid,pdf(GMModel,xgrid),'k');
M=[M max(pdf(GMModel,xgrid))];
hold on; plot(DSzone7{g,h}',pdf(GMModel,DSzone7{g,h}'),'r-','LineWidth',3);
hold on; plot(DEzone7{g,h}',pdf(GMModel,DEzone7{g,h}'),'b-','LineWidth',3);
AAA=[AAA DSzone7{g,h}];
if h==3
BBB2=[BBB2 DEzone7{g,h}];
else
BBB1=[BBB1 DEzone7{g,h}];
end
title3=sprintf('TM%d stage%d',g-1,h);
title3_1=sprintf('Treadmill %d Stage %d',g-1,h);
title(title3_1)
xlabel('Elapsed time (s)')
yyaxis left
ylabel('field Num')
xticks([0:250:size(ALL_timeseqF2{g,h},2)])
cc= [1.0000    0.8500    0.8500];
hold on;xline(250,'Color',cc,'LineWidth',5)
cc= [0.8500 0.3250 0.0980];
if h==3
hold on;xline(1250,'Color',cc,'LineWidth',5)
else
    hold on;xline(750,'Color',cc,'LineWidth',5)
end
set(gca, 'FontName', 'Times New Roman','FontSize',16)
FH=gca;
% printFig(gcf,title3);
end
end
end
for g=2:3
    for h=2:4
        if ~isempty(DSzone7{g,h})
            noDzone7{g,h}=max(DSzone7{g,h}):min(DEzone7{g,h});
        end
    end
    noDzone7{g+2,3}=floor(min(noDzone7{g,3})/2):floor(max(noDzone7{g,3})/2);
    DEzone7{g+2,3}=DEzone7{g,3}-500;
end
xticks(linspace(0,1500,7))
c=[0.4660 0.6740 0.1880];
M=max(M);
hold on;plot(mean(AAA),linspace(realmin,max(M*1.2),50),'MarkerEdgeColor',c,'Marker','.')
hold on;plot(mean(BBB1),linspace(realmin,max(M*1.2),50),'MarkerEdgeColor',c,'Marker','.')
hold on;plot(mean(BBB2),linspace(realmin,max(M*1.2),50),'MarkerEdgeColor',c,'Marker','.')
legend off
legend('TM1 stage2','TM1 stage3','TM1 stage4','TM2 stage2','TM2 stage3','TM2 stage4')
xlabel('elapsed time (20ms/bin)')
ylabel('percent of cells')
set(gca, 'FontName', 'Times New Roman','FontSize',16)
title3=sprintf('ALL_duration hist');
% plot ver.
figure;
A=[];B=[];x=1;C=[];
for i=2:3
    xx=1;
for j=2:4
B=size(DSzone7{i,j},2);
A=[A size(DSzone7{i,j},2)];
    if i==2
        c=[0.1+0.25*xx 0.1+0.25*xx 1];
    else
        c=[1 0.1+0.25*xx 0.1+0.25*xx];
    end
hold on;plot(B,1,'.','MarkerEdgeColor',c,'MarkerSize',50)
x=x+1;
end
end
C=A;A=[];B=[];x=1;
for i=2:3
    xx=1;
for j=2:4
B=size(DEzone7{i,j},2);
C=[C size(DEzone7{i,j},2)];
    if i==2
        c=[0.1+0.25*xx 0.1+0.25*xx 1];
    else
        c=[1 0.1+0.25*xx 0.1+0.25*xx];
    end
hold on;plot(B,2,'.','MarkerEdgeColor',c,'MarkerSize',50)
xx=xx+1;
end
end
A=[];B=[];x=1;
for i=2:3
    xx=1;
for j=2:4
    if j==3
        B=size(noDzone7{i+2,j},2);
    else
        B=size(noDzone7{i,j},2);
    end
    C=[C B];
    if i==2
        c=[0.1+0.25*xx 0.1+0.25*xx 1];
    else
        c=[1 0.1+0.25*xx 0.1+0.25*xx];
    end
hold on;plot(B,3,'.','MarkerEdgeColor',c,'MarkerSize',50)
xx=xx+1;
end
end
xlim([min(C)-50 max(C)+50])
ylim([0 4])
yticks(linspace(0,3,4))
legend('TM1 stage2','TM1 stage3','TM1 stage4','TM2 stage2','TM2 stage3','TM2 stage4')
xlabel('Field size (20ms/bin)')
set(gca, 'FontName', 'Times New Roman','FontSize',16)
title3=sprintf('ALL_duration plot');
printFig(gcf,title3);
save 20221206.mat -append DSzone7 DEzone7 DSpoint7 DEpoint7 noDzone7

% fig file title : zone significant(VS TMs)
% random data: XXXのdot数と同じ数のrandom dot
% shuffled data: random data * 100
% 各zone内のdot数/全dot数をttest2 (S/O/E*S/O/E)の９通り vs shuffled data
ShNum=10000;
% Ori_data=ALL_dot_VSstage2;
Ori_data=ALL_dot_VSTM2;
Sh_data=cell(size(Ori_data));
for i=1:size(Ori_data,1)
    for j=1:size(Ori_data,2)
        if ~isempty(Ori_data{i,j})
            if max(Ori_data{i,j}(:,1))<1000
                XX1=999;
            else
                XX1=1499;
            end
            if max(Ori_data{i,j}(:,2))<1000
                XX2=999;
            else
                XX2=1499;
            end
            for k=1:ShNum
                Sh_data{i,j,k}(:,1)=randi([1 XX1],length(Ori_data{i,j}),1);
                Sh_data{i,j,k}(:,2)=randi([1 XX2],length(Ori_data{i,j}),1);
            end
        end
    end
end
% 密度図に線引いて有意差が出たとこだけpatch
ccS= [1.0000    0.8500    0.8500];
ccE= [0.8500 0.3250 0.0980];
ccO=[1 1 1];
X=[];S=[];ALL_D=[];
ALL_SD=cell(9,1);
for k=1:size(Ori_data,1)
    for l=1:size(Ori_data,2)
        if ~isempty(Ori_data{k,l})
            X=[];
            figure;% SOE zoneを色分けしてpacth zone内に密度を記載　密度を残しとく
%             VS TMs
            title1=sprintf('Treadmill 1vs2 stage %d',l);
%             % VS stages
%             if l==1
%                 title1=sprintf('Stage2vs3 Treadmill %d',k-1);
%             elseif l==2
%                 title1=sprintf('Stage3vs4 Treadmill %d',k-1);
%             elseif l==3
%                 title1=sprintf('Stage2vs4 Treadmill %d',k-1);
%             end
            title2='Density (%)';
            title(sprintf('%s\n%s',title1,title2))
            axis equal
            if max(Ori_data{k,l}(:,1))>1000
                XX=1499;
            else
                XX=999;
            end
            if max(Ori_data{k,l}(:,2))>1000
                YY=1499;
            else
                YY=999;
            end
            xlim([0 XX])
            ylim([0 YY])
            hold on;plot(Ori_data{k,l}(:,1),Ori_data{k,l}(:,2),'.k')
            xticks([0:250:XX])
            yticks([0:250:YY])
            xxx=1;
            for i=1:3 % S/O/E zone
                if i==1
                    zone1=DSzone7;
                    cc1=ccS;
                elseif i==2
                    zone1=noDzone7;
                    cc1=ccO;
                elseif i==3
                    zone1=DEzone7;
                    cc1=ccE;
                else 
                end
                for j=1:3 % S/O/E zone
                    if j==1
                        zone2=DSzone7;
                        cc2=ccS;
                    elseif j==2
                        zone2=noDzone7;
                        cc2=ccO;
                    elseif j==3
                        zone2=DEzone7;
                        cc2=ccE;
                    else 
                    end
%                     % VS stages
%                     if l==1
%                         zone3=zone1{k,2};
%                         zone4=zone2{k,3};                      
%                     elseif l==2
%                         zone3=zone1{k,4};
%                         zone4=zone2{k,3};                       
%                     elseif l==3
%                         zone3=zone1{k,2};
%                         zone4=zone2{k,4};
%                     end                    
%                     VS TMs
                    zone3=zone1{2,l};
                    zone4=zone2{3,l};  

                    XX2 = [zone3(1) zone3(1) zone3(end) zone3(end)];
                    YY2 = [zone4(1) zone4(1) zone4(end) zone4(end)];
                    if max(Ori_data{k,l}(:,1))<1000
                        XX1 = [0 999 999 0];
                    else
                        XX1 = [0 1499 1499 0];
                    end
                    if max(Ori_data{k,l}(:,2))<1000
                        YY1 = [0 999 999 0];
                    else
                        YY1 = [0 1499 1499 0];
                    end
                        A=find(Ori_data{k,l}(:,1)>zone3(1)-1 & Ori_data{k,l}(:,1)<zone3(end)+1);
                        B=find(Ori_data{k,l}(:,2)>zone4(1)-1 & Ori_data{k,l}(:,2)<zone4(end)+1);
                        x=length(intersect(A,B))/length(Ori_data{k,l});
                        X=vertcat(X,x);
                        S=[];
                        for m=1:ShNum
                            if m==1
                                if max(Ori_data{k,l}(:,1))<1000 && max(Sh_data{k,l,m}(:,1))<1000

                                elseif max(Ori_data{k,l}(:,1))>1000 && max(Sh_data{k,l,m}(:,1))>1000
                                else
                                    matigai
                                end
                            end
                            SA=find(Sh_data{k,l,m}(:,1)>zone3(1)-1 & Sh_data{k,l,m}(:,1)<zone3(end)+1);
                            SB=find(Sh_data{k,l,m}(:,2)>zone4(1)-1 & Sh_data{k,l,m}(:,2)<zone4(end)+1);
                            Sx=length(intersect(SA,SB))/length(Sh_data{k,l,m});
                            S=vertcat(S,Sx);
                        end
                        SS=sort(S);
                        ALL_SD{xxx}=vertcat(ALL_SD{xxx},S);
                        xxx=xxx+1;
                        AA=knnsearch(SS,x);
                        % /9はbonferroni
                        if AA>ShNum*(1-0.01/9)
                            formatSpec = '%0.1d***';
%                             formatSpec = '%0.1d***\n(%0.1d %d)';
                        elseif AA>ShNum*(1-0.1/9)
                            formatSpec = '%0.1d**';
%                             formatSpec = '%0.1d**\n(%0.1d %d)';
                        elseif AA>ShNum*(1-0.5/9)
                            formatSpec = '%0.1d*';
%                             formatSpec = '%0.1d*\n(%0.1d %d)';
                        else
                            formatSpec = '%0.1d';
%                             formatSpec = '%0.1d\n(%0.1d %d)';
                        end

                    hold on;patch(XX2,YY1,cc1,'FaceAlpha',0.1,'LineStyle','none')
                    hold on;patch(XX1,YY2,cc2,'FaceAlpha',0.1,'LineStyle','none')
                    txt=sprintf(formatSpec,x*100);
%                     txt=sprintf(formatSpec,x*100,mean(SS)*100,knnsearch(SS,x));
                    hold on;text(mean(XX2),mean(YY2),txt,'FontSize',14,'HorizontalAlignment','center','FontName', 'Times New Roman','Color',[0 0.4470 0.7410],'FontWeight','bold');
                end
            end
%             sum(ALL_D{k,l})
        h = get(gca,'Children');
        hg = findobj(h,'type','text');
        for n=1:length(hg)
            ind = (h == hg(end));
            newh = [h(ind); h(~ind)];
            set(gca,'children',newh)
            h = get(gca,'Children');
            hg = findobj(h,'type','text');
        end

        xlabel('Elapsed time (s)')
        ylabel('Elapsed time (s)')

        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        title3=sprintf('Dot density (%s)',title1);
        printFig(gcf,title3);
        ALL_D(:,end+1)=X;
        end
    end
end
figure
XX=999;YY=999;
xticks([0:250:XX])
yticks([0:250:YY])
ccS= [1.0000    0.8500    0.8500];
ccE= [0.8500 0.3250 0.0980];
XX2 = [150 150 350 350];
YY2 = [750 750 950 950];
XX1 = [0 999 999 0];
YY1 = [0 999 999 0];
hold on;patch(XX2,YY1,ccS,'FaceAlpha',0.2,'LineStyle','none')
hold on;patch(YY1,XX2,ccS,'FaceAlpha',0.2,'LineStyle','none')
hold on;patch(XX1,YY2,ccE,'FaceAlpha',0.2,'LineStyle','none')
hold on;patch(YY2,XX1,ccE,'FaceAlpha',0.2,'LineStyle','none')
title('Mean Density (%)')
XY={[250,250];[250,500];[250,850];[500,250];[500,500];[500,850];[850,250];[850,500];[850,850]};
for i=1:9
    [h,p]=ttest2(ALL_D(i,:),ALL_SD{i});
    p=p*9;%Bonferroni
    if h==1% 有意差があれば*追加
        if p>.001
            formatSpec = '%0.1d***';
        elseif p>.01
            formatSpec = '%0.1d**';
        elseif p<.05
            formatSpec = '%0.1d*';
        else
            matigai
        end
    end
    txt=sprintf(formatSpec,mean(ALL_D(i,:))*100);
    hold on;text(XY{i}(1),XY{i}(2),txt,'FontSize',14,'HorizontalAlignment','center','FontName', 'Times New Roman','Color',[0 0.4470 0.7410],'FontWeight','bold')
end
xlabel('Elapsed time (s)')
ylabel('Elapsed time (s)')
xlim([0 XX])
ylim([0 YY])
axis equal
set(gca, 'FontName', 'Times New Roman','FontSize',16)
% % vs stages
% title3=sprintf('ALL_Dot density (vs stages)');
% vsTMs
title3=sprintf('ALL_Dot density (vs TMs)');
        printFig(gcf,title3);




%%%各々
% DE/DS zoneのglm
% DS/DE pointを中心にDS/DE zone sizeのT,P,Dseqを作ってglm
event_Num{1}=PORK3;
event_Num{2}=PORK4;
[event_Num{3},event_Num{5},DS1]=gettimesM(event(1,:),2.5);
[event_Num{4},event_Num{6},DS2]=gettimesM(event(2,:),2.5);
dly1=[];dly2=[];TMtimes=[];
for i=1:size(event_Num{3},2)
[~,I3]=min(abs(event_Num{3}(i)-PosT));
[~,I5]=min(abs(event_Num{5}(i)-PosT));
[~,I4]=min(abs(event_Num{4}(i)-PosT));
[~,I6]=min(abs(event_Num{6}(i)-PosT));
if isempty(I3:I5) || isempty(I4:I6)
mieno
else
dly1=[dly1 I3:I5];
dly2=[dly2 I4:I6];
end
end
TMtimes=[dly1 dly2];
TMzone_x=cell(1,2);TMzone_y=cell(1,2);
TMzone_x{1}=[ceil(min(Traj(dly1,3))-100):ceil(max(Traj(dly1,3))+100)];
TMzone_x{2}=[ceil(min(Traj(dly2,3))-100):ceil(max(Traj(dly2,3))+100)];
TMzone_y{1}=[ceil(min(Traj(dly1,4))-50):ceil(max(Traj(dly1,4))+50)];
TMzone_y{2}=[ceil(min(Traj(dly2,4))-50):ceil(max(Traj(dly2,4))+50)];
% TMzone=[min(LTraj(dly1))-100:max(LTraj(dly1))+100 min(LTraj(dly2))-100:max(LTraj(dly2))+100];
[~,I1]=min(abs(event_Num{1}-PosT));
[~,I2]=min(abs(event_Num{2}-PosT));
Fzone_x=cell(1,2);Fzone_y=cell(1,2);
Fzone_x{1}=[ceil(min(Traj(I1,3))-100):ceil(max(Traj(I1,3))+100)];
Fzone_x{2}=[ceil(min(Traj(I2,3))-100):ceil(max(Traj(I2,3))+100)];
Fzone_y{1}=[ceil(min(Traj(I1,4))-50):ceil(max(Traj(I1,4))+50)];
Fzone_y{2}=[ceil(min(Traj(I2,4))-50):ceil(max(Traj(I2,4))+50)];
% Fzone=[min(LTraj(I1))-100:max(LTraj(I1))+100 min(LTraj(I2))-100:max(LTraj(I2))+100];
addpath('C:\Users\dddd1002\OneDrive - 同志社大学\20220307\storage6_2')
load 20221206.mat DSzone7 DEzone7 noDzone7 DSpoint7 DEpoint7
DSseqF=cell(3,4);DEseqF=cell(3,4);
DSseqF_ori=cell(3,4);DEseqF_ori=cell(3,4);
DSpoint=cell(3,4);DEpoint=cell(3,4);
DSpoint2=cell(3,4);DEpoint2=cell(3,4);
kekkaS=zeros(size(OK_ens,2),3);kekkaE=zeros(size(OK_ens,2),3);
GLMS=cell(3,4);GLME=cell(3,4);
within_cellS=cell(3,4);within_cellE=cell(3,4);within_cellNoSE=cell(3,4);
TimeCellF=cell(3,4);
DSCellF=cell(3,4);
DECellF=cell(3,4);
noDCellF=cell(3,4);
for i=1:3
    if i==2
        TT=3;
        EE=4;
        DD=6;
    elseif i==3
        TT=2;
        EE=3;
        DD=5;
    end
    for j=1:4
        if j==3
            MM=1250;
        elseif j==2 || j==4
            MM=750;
        end
        if ~isempty(DSpoint7{i,j})
            for k=1:10
                %DS/Epoint,DS/Epoint2:DS/DE pointのevent_Num ver.と Traj_Num ver.
                DSpoint{i,j}(k)=event_Num{EE}(10*(j-2)+k)+(DSpoint7{TT,j}-250)*25000/50;
                DEpoint{i,j}(k)=event_Num{DD}(10*(j-2)+k)+(DEpoint7{TT,j}-MM)*25000/50;
                DSpoint2{i,j}=[DSpoint2{i,j} find(PosT>(event_Num{EE}(10*(j-2)+k)+(min(DSzone7{TT,j})-250)*25000/50) & PosT<(event_Num{EE}(10*(j-2)+k)+(max(DSzone7{TT,j})-250)*25000/50))'];
                DEpoint2{i,j}=[DEpoint2{i,j} find(PosT>(event_Num{EE}(10*(j-2)+k)+(min(DEzone7{TT,j})-250)*25000/50) & PosT<(event_Num{EE}(10*(j-2)+k)+(max(DEzone7{TT,j})-250)*25000/50))'];
            end
            rS=DSpoint2{i,j};
            PosTm1S=PosT(rS);
            LTrajSm1S=LTrajS(rS,:);
            rE=DEpoint2{i,j};
            PosTm1E=PosT(rE);
            LTrajSm1E=LTrajS(rE,:);
            for k=1:size(OK_ens,2)
                [histR1,~,HistR1]=plotRasterM(ensemble2{OK_ens(k),3},DSpoint{i,j},5.5);
                [histR2,~,HistR2]=plotRasterM(ensemble2{OK_ens(k),3},DEpoint{i,j},5.5);
                    DSseqF{i,j}(k,:)=HistR1;
                    DEseqF{i,j}(k,:)=HistR2;
                    DSseqF_ori{i,j}(k,:)=histR1;
                    DEseqF_ori{i,j}(k,:)=histR2;
                    HistR1=[];HistR2=[];
                [~,rasterS,HistRTS]=plotRasterM(ensemble2{OK_ens(k),3},DSpoint{i,j},size(DSzone7{i,j},2)/50/2);
                [~,rasterE,HistRTE]=plotRasterM(ensemble2{OK_ens(k),3},DEpoint{i,j},size(DEzone7{i,j},2)/50/2);
                [~,rasterNoSE]=plotRasterM(ensemble2{OK_ens(k),3},mean([DSpoint{i,j};DEpoint{i,j}],1),ceil(((min(DEzone7{i,j})-max(DSzone7{i,j}))*20)/1000)/2);
                for l=1:3
                    if l==1
                        A=sum(rasterS(3:6,:));
                        B=sum(rasterS(7:10,:));
                    elseif l==2
                        A=sum(rasterE(3:6,:));
                        B=sum(rasterE(7:10,:));
                    elseif l==3
                        A=sum(rasterNoSE(3:6,:));
                        B=sum(rasterNoSE(7:10,:));
                    end
                    A1=[];B1=[];
                    for m=26:length(A)-25
                        A1(1,m-25)=sum(A(1,m-25:m+25));
                    end
                    for m=26:length(B)-25
                        B1(1,m-25)=sum(B(1,m-25:m+25));
                    end
                    if corr(A1',B1')>.3
                        if l==1
                            within_cellS{i,j}=[within_cellS{i,j} k];
                        elseif l==2
                            within_cellE{i,j}=[within_cellE{i,j} k];
                        elseif l==3
                            within_cellNoSE{i,j}=[within_cellNoSE{i,j} k];
                        end
                    end
                    A1=[];B1=[];
                end
                [~,HistRDS]=plotRasterMMD(ensemble2{OK_ens(k),3},(DSpoint{i,j}-(size(DSzone7{i,j},2)/50/2)*25000),(DSpoint{i,j}+(size(DSzone7{i,j},2)/50/2)*25000),Traj,PosT,TMtimes);
                [~,HistRDE]=plotRasterMMD(ensemble2{OK_ens(k),3},(DEpoint{i,j}-(size(DEzone7{i,j},2)/50/2)*25000),(DEpoint{i,j}+(size(DEzone7{i,j},2)/50/2)*25000),Traj,PosT,TMtimes);
%                 Len=max(LTraj);
                [HistRPS]=pmapK(ensemble2{OK_ens(k),3},LTrajSm1S,PosTm1S,0,'animal','rat');
                [HistRPE]=pmapK(ensemble2{OK_ens(k),3},LTrajSm1E,PosTm1E,0,'animal','rat');
                if size(HistRPS,2)==1
                    HistRPS=HistRPS';
                end
                if size(HistRPE,2)==1
                    HistRPE=HistRPE';
                end
            [kekka]=GLMm2(HistRTS,HistRDS,HistRPS);
            kekkaS(k,1:size(kekka,2))=kekka;
%             figure;
%             hold on;plot(HistRTS)
%             hold on;plot(HistRDS)
%             hold on;plot(HistRPS')           
            HistRTS=[];HistRDS=[];HistRPS=[];
            [kekka]=GLMm2(HistRTE,HistRDE,HistRPE);
            kekkaE(k,1:size(kekka,2))=kekka;
%             figure;
%             hold on;plot(HistRTE)
%             hold on;plot(HistRDE)
%             hold on;plot(HistRPE')
            HistRTE=[];HistRDE=[];HistRPE=[];
            MM=2;
            [out,xy,peak_Point]=countPlaceFieldsM(time_seqF2{i,j}(k,:),1,50,MM);
            T_index=[];
            if ~isempty(xy)
            for l=1:size(xy,2)
                if out(l)<500
                index=xy{l}(:,1)';
                T_index=[T_index peak_Point];
        %         T_index=[T_index index];
                end
            end
            if ~isempty(intersect(T_index,DSzone7{i,j})) && ~isempty(intersect(kekkaS(k,:),1)) && ~isempty(intersect(within_cellS{i,j},k)) && ~isempty(intersect(k,DlyCellF{i,j}))
                DSCellF{i,j}=[DSCellF{i,j} k];
            end
            if ~isempty(intersect(T_index,DEzone7{i,j})) && ~isempty(intersect(kekkaE(k,:),1)) && ~isempty(intersect(within_cellE{i,j},k)) && ~isempty(intersect(k,DlyCellF{i,j}))
                DECellF{i,j}=[DECellF{i,j} k];
            end
            if ~isempty(intersect(T_index,noDzone7{i,j})) && ~isempty(intersect(within_cellNoSE{i,j},k))  && ~isempty(intersect(k,DlyCellF{i,j}))
                noDCellF{i,j}=[noDCellF{i,j} k];
            end
            end
            end
            GLMS{i,j}=kekkaS;
            GLME{i,j}=kekkaE;
            kekkaS=zeros(size(OK_ens,2),3);kekkaE=zeros(size(OK_ens,2),3);
            TimeCellF{i,j}=union(DSCellF{i,j},union(DECellF{i,j},noDCellF{i,j}));
            if size(TimeCellF{i,j},1)>size(TimeCellF{i,j},2)
                TimeCellF{i,j}=TimeCellF{i,j}';
            end 
        end
    end
end
save 20220318.mat -append DSseqF DEseqF GLMS GLME within_cellS within_cellE within_cellNoSE DSseqF_ori DEseqF_ori TimeCellF DSCellF DECellF noDCellF
%%%各々 PlaceCellとnoDSEseq作る
% stage1でplace cell認定されたcell(maze, TM) within,info,peak,field
%fig title : pmap
PlaceCellF=cell(1,4);PlaceCellFonTM=cell(1,4);PlaceCellFonF=cell(1,4);
addpath('C:\Users\dddd1002\OneDrive - 同志社大学\20220307\storage6_2')
load 20221206.mat DSzone7 DEzone7 noDzone7
noDSEseqF=cell(1,4);PseqF=cell(1,4);
for h=1:4%stage
     r1=[];r2=[];r=[];
    Nsta1=[1 11 21 31];Nsta2=[6 16 26 36];Nsta=[1 11 21 31];
    Nen1=[6 16 26 36];Nen2=[11 21 31 size(PORK3,2)];Nen=[11 21 31 size(PORK3,2)];
    r1=find(PosT>PORK3(Nsta1(h)) & PosT<PORK3(Nen1(h)));
    r2=find(PosT>PORK3(Nsta2(h)) & PosT<PORK3(Nen2(h)));
    r=find(PosT>PORK3(Nsta(h)) & PosT<PORK3(Nen(h)));
    PosTm1=PosT(r1);PosTm2=PosT(r2);PosTm=PosT(r);
    LTrajSm1=LTrajS(r1,:);LTrajSm2=LTrajS(r2,:);LTrajSm=LTrajS(r,:);
    TrajSm=Traj(r,:);
    if h==1
    else
        % JPre,JPost DS/DEzoneのサイズ分だけのjitter
        if h==3
        JPost1=(max(DEzone7{3,h})-1250)*20/1000;%50fps eventとseqの順番は逆
        JPost2=(max(DEzone7{2,h})-1250)*20/1000;
        else
        JPost1=(max(DEzone7{3,h})-750)*20/1000;
        JPost2=(max(DEzone7{2,h})-750)*20/1000;
        end 
        JPre1=(250-min(DSzone7{3,h}))*20/1000;
        JPre2=(250-min(DSzone7{2,h}))*20/1000;
    end
    for i=1:size(OK_ens,2)
        HistR1=[];HistR2=[];HistR=[];oc_map=[];
        if h==1
        [HistR1]=pmapK(ensemble2{OK_ens(i),3},LTrajSm1,PosTm1,0,'animal','rat');
        [HistR2]=pmapK(ensemble2{OK_ens(i),3},LTrajSm2,PosTm2,0,'animal','rat');
        [LHistR]=pmapK(ensemble2{OK_ens(i),3},LTrajSm,PosTm,0,'animal','rat');
        %LHistRS　1000回シャッフルしてスムージングも終わったLHistR
        [LHistRS]=pmapK(ensemble2{OK_ens(i),3},LTrajSm,PosTm,0,'animal','rat','shuffle',1);
        [HistR, ~, ~, oc_map]=pmapK(ensemble2{OK_ens(i),3},TrajSm,PosTm,0,'animal','rat');
        noDSEseqF{1,h}(i,:)=LHistR;
        PseqF{1,h}(i,:)=LHistR;
        else
        [HistR1]=pmapK(extractDelayM(ensemble2{OK_ens(i),3},event,[1 2],'jitterpre1',JPre1,'jitterpre2',JPre2,'jitterpost1',JPost1,'jitterpost2',JPost2),LTrajSm1,PosTm1,0,'animal','rat');    
        [HistR2]=pmapK(extractDelayM(ensemble2{OK_ens(i),3},event,[1 2],'jitterpre1',JPre1,'jitterpre2',JPre2,'jitterpost1',JPost1,'jitterpost2',JPost2),LTrajSm2,PosTm2,0,'animal','rat');
        [LHistR]=pmapK(extractDelayM(ensemble2{OK_ens(i),3},event,[1 2],'jitterpre1',JPre1,'jitterpre2',JPre2,'jitterpost1',JPost1,'jitterpost2',JPost2),LTrajSm,PosTm,0,'animal','rat');
        [LHistR1]=pmapK(ensemble2{OK_ens(i),3},LTrajSm,PosTm,0,'animal','rat');
        [LHistRS]=pmapK(extractDelayM(ensemble2{OK_ens(i),3},event,[1 2],'jitterpre1',JPre1,'jitterpre2',JPre2,'jitterpost1',JPost1,'jitterpost2',JPost2),LTrajSm,PosTm,0,'animal','rat','shuffle',1);
        [HistR,~,~,oc_map]=pmapK(extractDelayM(ensemble2{OK_ens(i),3},event,[1 2],'jitterpre1',JPre1,'jitterpre2',JPre2,'jitterpost1',JPost1,'jitterpost2',JPost2),TrajSm,PosTm,0,'animal','rat');
        noDSEseqF{1,h}(i,:)=LHistR;
        PseqF{1,h}(i,:)=LHistR1;
        end
        if size(HistR1,1)==1
            HistR1=HistR1';
        else
        end
        if size(HistR2,1)==1
            HistR2=HistR2';
        else
        end
        R=[];P=[];I=[];xy=[];
        R=corr(HistR1,HistR2);
        P=max(LHistR);
        I=calcInfoM(LHistR);

        testS1=[];
        seqS3=LHistRS;
        for g=1:1000
            It=calcInfoM(seqS3(g,:));
            testS1=[testS1 It];
        end
        A=sort(testS1);
        %要確認　peakの3/5をfield(2.5cm/bin)　field size>50cm2であればPFとする
        %MM=1にして問題なければこれで
        MM=1;ThF=50;
        [~,xy]=countPlaceFieldsM(HistR,2.5,ThF,MM);
        x=[];y=[];
%         title0='';
        if I<A(50) && R>0.3 && P>1 && ~isempty(xy)
            PlaceCellF{h}=[PlaceCellF{h} i];
%             title0='PlaceCell';
            x=0;y=0;
            for j=1:size(xy,2)
                if intersect(xy{j}(:,1),[20:40])
                    PlaceCellFonTM{h}=[PlaceCellFonTM{h} i];
                    x=1;
%                     title0='PlaceCellonTM';
                end
                if  ~isempty(intersect(xy{j}(:,1),[1:10 size(HistR,2)-9:size(HistR,2)])) && ~isempty(intersect(xy{j}(:,2),ceil(size(HistR,1)/2)-4:ceil(size(HistR,1)/2)+5))
                    PlaceCellFonF{h}=[PlaceCellFonF{h} i];
                    y=1;
%                     title0='PlaceCellonF';
                end
            end
            PlaceCellFonTM{h}=union(PlaceCellFonTM{h},PlaceCellFonTM{h});
            PlaceCellFonF{h}=union(PlaceCellFonF{h},PlaceCellFonF{h});
        end
        if isempty(x) && isempty(y)
             title0='';
        else
            if x==1 && y==1
                title0='PlaceCellonTMandF';
            elseif x==1
                title0='PlaceCellonTM';
            elseif y==1
                title0='PlaceCellonF';
            elseif x==0 && y==0
                title0='PlaceCell';
            else
%                 mieno
            end
        end
%         figure;imagePmap(HistR,oc_map);colorbar
%         set(gca, 'FontName', 'Times New Roman','FontSize',16)
%         title3=sprintf('stage%d-cell%d %s',h,i,title0);
%         printFig(gcf,title3);
    end
end
save 20220318.mat -append PlaceCellF PlaceCellFonTM noDSEseqF PlaceCellFonF PseqF
% DS/DE cellsのPseqがstage1とstage2で変化するのか
cells=union(DSCellF{2,2},DSCellF{3,2});
if size(cells,1)<size(cells,2)
else
    cells=cells';
end
[~,~,TMarea]=fitLinearM(Traj);
TMarea=ceil(TMarea/10);
% cells=union(DECellF{2,2},DECellF{3,2});
for i=cells
    figure;
    hold on;plot(PseqF{1,1}(i,:))
    MM=2;
    ThF=30;
    [~,xy]=countPlaceFieldsM(PseqF{1,1}(i,:),2.5,ThF,MM);
    for j=1:size(xy,2)
        hold on;plot(xy{j}(:,1),PseqF{1,1}(i,xy{j}(:,1)),'LineWidth',2)
    end
    hold on;plot(PseqF{1,2}(i,:))
    [~,xy]=countPlaceFieldsM(PseqF{1,2}(i,:),2.5,ThF,MM);
    for j=1:size(xy,2)
        hold on;plot(xy{j}(:,1),PseqF{1,2}(i,xy{j}(:,1)),'LineWidth',2)
    end
    M=max([PseqF{1,1}(i,:) PseqF{1,2}(i,:)]);
    %TMがあるところを塗りつぶす
%     hold on;plot(TMarea(1),linspace(realmin,max(M*1.2),50),'r','Marker','|')
%     hold on;plot(TMarea(2),linspace(realmin,max(M*1.2),50),'r','Marker','|')
    XX=[TMarea(1) TMarea(2) TMarea(2) TMarea(1)];
    YY=[0 0 M*1.2 M*1.2];
    c=[0.4660 0.6740 0.1880];
    hold on;patch(XX,YY,c,'FaceAlpha',.5,'EdgeColor','none')
    XX=[TMarea(3) TMarea(4) TMarea(4) TMarea(3)];
    hold on;patch(XX,YY,c,'FaceAlpha',.5,'EdgeColor','none')
end

%各々
% DS/DE zoneの走行速度
event_Num{1}=PORK3;
event_Num{2}=PORK4;
[event_Num{3},event_Num{5},DS1]=gettimesM(event(1,:),2.5);
[event_Num{4},event_Num{6},DS2]=gettimesM(event(2,:),2.5);
dly1=[];dly2=[];TMtimes=[];
for i=1:size(event_Num{3},2)
[~,I3]=min(abs(event_Num{3}(i)-PosT));
[~,I5]=min(abs(event_Num{5}(i)-PosT));
[~,I4]=min(abs(event_Num{4}(i)-PosT));
[~,I6]=min(abs(event_Num{6}(i)-PosT));
if isempty(I3:I5) || isempty(I4:I6)
mieno
else
dly1=[dly1 I3:I5];
dly2=[dly2 I4:I6];
end
end
TMtimes=[dly1 dly2];
dly1=[];dly2=[];DzoneTimes=[];
[~,I3]=min(abs(event_Num{3}-PosT));
[~,I5]=min(abs(event_Num{5}-PosT));
[~,I4]=min(abs(event_Num{4}-PosT));
[~,I6]=min(abs(event_Num{6}-PosT));
for i=1:size(event_Num{3},2)
if intersect(11:20,i)
h=3;
pre1=(250-min(DSzone7{3,h}));
pre2=(250-min(DSzone7{2,h}));
post1=(max(DEzone7{3,h})-1250);
post2=(max(DEzone7{2,h})-1250);
elseif intersect(21:30,i)
h=4;
pre1=(250-min(DSzone7{3,h}));
pre2=(250-min(DSzone7{2,h}));
post1=(max(DEzone7{3,h})-750);
post2=(max(DEzone7{2,h})-750);
elseif intersect(1:10,i)
h=2;
pre1=(250-min(DSzone7{3,h}));
pre2=(250-min(DSzone7{2,h}));
post1=(max(DEzone7{3,h})-750);
post2=(max(DEzone7{2,h})-750);
end
if isempty(I3:I5) || isempty(I4:I6)
mieno
else
dly1=[dly1 floor(I3(i)-pre1):floor(I5(i)+post1)];
dly2=[dly2 floor(I4(i)-pre2):floor(I6(i)+post2)];
end
end
DzoneTimes=[dly1 dly2];
for i=length(DzoneTimes):-1:1
if DzoneTimes(i)>length(Traj)
DzoneTimes(:,i)=[];
end
end
All_sp=cell(2,3);
[all_sp]=speed_extract(1,length(event),speed,Traj,PosT,TMtimes,DzoneTimes);
event_Num1=[];event_Num2=[];
for i=1:3
    for j=1:2
        if j==1
            event_Num1=event_Num{3}(((i-1)*10+1):(i*10))-25000*15;
            event_Num2=event_Num{5}(((i-1)*10+1):(i*10))+25000*15;
        elseif j==2
            event_Num1=event_Num{4}(((i-1)*10+1):(i*10))-25000*15;
            event_Num2=event_Num{6}(((i-1)*10+1):(i*10))+25000*15;
        end
        [All_sp{j,i}]=speed_extract(event_Num1,event_Num2,all_sp,Traj,PosT,TMtimes,DzoneTimes);
        figure;plot(All_sp{j,i})
        c=[0.4660 0.6740 0.1880];
        M=max(All_sp{j,i});
        c=[0.4660 0.6740 0.1880];
        XX=[DSzone7{j+1,i+1}(1) DSzone7{j+1,i+1}(end) DSzone7{j+1,i+1}(end) DSzone7{j+1,i+1}(1)];
        YY=[0 0 M*1.2 M*1.2];
        hold on;patch(XX+500,YY,c,'FaceAlpha',.5,'EdgeColor','none')
        XX=[DEzone7{j+1,i+1}(1) DEzone7{j+1,i+1}(end) DEzone7{j+1,i+1}(end) DEzone7{j+1,i+1}(1)];
        YY=[0 0 M*1.2 M*1.2];
        hold on;patch(XX+500,YY,c,'FaceAlpha',.5,'EdgeColor','none')
        xlabel('elapsed time (20ms/bin)')
        ylabel('running speed(20ms/bin)')
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        title3=sprintf('stage%d TM%d',i+1,j);
        printFig(gcf,title3);
    end
end
save 20220318.mat -append All_sp

%%%全部
x=0;
ALL_DSseqF=cell(3,4);ALL_DEseqF=cell(3,4);
ALL_DSseqF_ori=cell(3,4);ALL_DEseqF_ori=cell(3,4);
ALL_GLMS=cell(3,4);ALL_GLME=cell(3,4);
ALL_within_cellS=cell(3,4);ALL_within_cellE=cell(3,4);ALL_within_cellNoSE=cell(3,4);
ALL_TimeCellF=cell(3,4);ALL_DSCellF=cell(3,4);ALL_DECellF=cell(3,4);ALL_noDCellF=cell(3,4);
for i=1:5
addpath(PathName{i})
load 20220318.mat DSseqF DEseqF GLMS GLME within_cellS within_cellE within_cellNoSE DSseqF_ori DEseqF_ori TimeCellF DSCellF DECellF noDCellF
if i==1
x=0;
else
x=x+X{i-1};
end
for j=1:4
    for k=1:3
        if ~isempty(DSseqF{k,j})
            ALL_DSseqF{k,j}=vertcat(ALL_DSseqF{k,j},DSseqF{k,j});
            ALL_DSseqF_ori{k,j}=vertcat(ALL_DSseqF_ori{k,j},DSseqF_ori{k,j});
        end
        if ~isempty(DEseqF{k,j})
            ALL_DEseqF{k,j}=vertcat(ALL_DEseqF{k,j},DEseqF{k,j});
            ALL_DEseqF_ori{k,j}=vertcat(ALL_DEseqF_ori{k,j},DEseqF_ori{k,j});
        end
        if ~isempty(GLMS{k,j})
            ALL_GLMS{k,j}=vertcat(ALL_GLMS{k,j},GLMS{k,j});
        end
        if ~isempty(GLME{k,j})
            ALL_GLME{k,j}=vertcat(ALL_GLME{k,j},GLME{k,j});
        end
        if ~isempty(within_cellS{k,j})
            ALL_within_cellS{k,j}=[ALL_within_cellS{k,j} within_cellS{k,j}+x];
        end
        if ~isempty(within_cellE{k,j})
            ALL_within_cellE{k,j}=[ALL_within_cellE{k,j} within_cellE{k,j}+x];
        end
        if ~isempty(within_cellNoSE{k,j})
            ALL_within_cellNoSE{k,j}=[ALL_within_cellNoSE{k,j} within_cellNoSE{k,j}+x];
        end
        if ~isempty(TimeCellF{k,j})
            ALL_TimeCellF{k,j}=[ALL_TimeCellF{k,j} TimeCellF{k,j}+x];
        end
        if ~isempty(DSCellF{k,j})
            ALL_DSCellF{k,j}=[ALL_DSCellF{k,j} DSCellF{k,j}+x];
        end
        if ~isempty(DECellF{k,j})
            ALL_DECellF{k,j}=[ALL_DECellF{k,j} DECellF{k,j}+x];
        end
        if ~isempty(noDCellF{k,j})
            ALL_noDCellF{k,j}=[ALL_noDCellF{k,j} noDCellF{k,j}+x];
        end
    end
end
end
save 20221206.mat -append ALL_DSseqF ALL_DSseqF_ori ALL_DEseqF ALL_DEseqF_ori ALL_GLMS ALL_GLME ALL_within_cellS ALL_within_cellE ALL_within_cellNoSE
save 20221206.mat -append ALL_TimeCellF ALL_DSCellF ALL_DECellF ALL_noDCellF
%noDseq:timeseqからDSpointとDEpointの平均値から±jitter 
ALL_noDseqF=cell(3,4);
for i=2:3
    for j=2:4
        noDpoint7=floor(mean(noDzone7{i,j}));
%         noDzone=(noDpoint7-249):(noDpoint7+250);
        ALL_noDseqF{i,j}=ALL_timeseqF2{i,j}(:,noDzone);
    end
end


% zone毎にwithin corrの比較
for i=2:3
    for j=2:4
        C=[];
        for l=1:3 %1:DSzone, 2:noDzone, 3:DEzone
            if l==1
                zone=DSzone7{i,j};
            elseif l==2
                zone=noDzone7{i,j};
            elseif l==3
                zone=DEzone7{i,j};
            end
            x=1;
            % ALL_TimeCellF{i,j}にすると条件によって結果がブレる
            for k=ALL_DlyCellF{i,j}
                c=corr(ALL_timeseqF3{i,j}(k,zone)',ALL_timeseqF4{i,j}(k,zone)');
                C(x,l)=c;
                x=x+1;
            end
        end
        [p,~,stats]=anova1(C);
        xticklabels({'ESzone','OEzone','EEzone'})
        ylabel('correlation')
        axis square
        title2=sprintf('p = %d',p);
        title(title2)
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        title3=sprintf('corr VSzone(TM%d stage%d)',i-1,j);
%         printFig(gcf,title3);
        figure;c=multcompare(stats);
        title3_1=sprintf('1vs2 %d \n 1vs3 %d \n 2vs3 %d',c(1,6),c(2,6),c(3,6));
        title(title3_1)
        yticklabels({'ESzone','OEzone','EEzone'})
        axis square
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        title4=sprintf('corr VSzone(TM%d stage%d) mult',i-1,j);
%         printFig(gcf,title4);
        % [A,B]=FRstability(ALL_timeseqF2{2,j},ALL_timeseqF2{3,j},ALL_timeseqF3{2,j},ALL_timeseqF4{2,j},ALL_TimeCellF{2,j},2);
    end
end

% zone毎にbetween corrの比較
% TM1 vs TM2 →有意差なし　DlycellをTimecellに変えても結果は同じ
for j=2:4
    C=[];
    for l=1:3 %1:DSzone, 2:noDzone, 3:DEzone
        if l==1
            zone=intersect(DSzone7{2,j},DSzone7{3,j});
        elseif l==2
            zone=intersect(noDzone7{2,j},noDzone7{3,j});
        elseif l==3
            zone=intersect(DEzone7{2,j},DEzone7{3,j});
        end
        x=1;
        for k=intersect(ALL_DlyCellF{2,j},ALL_DlyCellF{3,j})
            c=corr(ALL_timeseqF2{2,j}(k,zone)',ALL_timeseqF2{3,j}(k,zone)');
            C(x,l)=c;
            x=x+1;
        end
    end
    [p,~,stats]=anova1(C);
    xticklabels({'ESzone','OEzone','EEzone'})
    ylabel('correlation')
    axis square
    title2=sprintf('p = %d',p);
    title(title2)
    set(gca, 'FontName', 'Times New Roman','FontSize',16)
    title3=sprintf('corr VSzone(stage%d)',j);
        printFig(gcf,title3);
    figure;c=multcompare(stats);
    title3_1=sprintf('1vs2 %d \n 1vs3 %d \n 2vs3 %d',c(1,6),c(2,6),c(3,6));
    title(title3_1)
    yticklabels({'ESzone','OEzone','EEzone'})
    axis square
    set(gca, 'FontName', 'Times New Roman','FontSize',16)
    title4=sprintf('corr VSzone(stage%d) mult',j);
        printFig(gcf,title4);
end
% stage 2 vs stages →　有意差あり
for i=2:3
    for j=3:4
        C=[];
        if j==3
            I=i+2;
        else 
            I=i;
        end
        for l=1:3 %1:DSzone, 2:noDzone, 3:DEzone
            if l==1
                zone=intersect(DSzone7{i,2},DSzone7{i,j});
            elseif l==2
                zone=intersect(noDzone7{i,2},noDzone7{I,j});
            elseif l==3
                zone=intersect(DEzone7{i,2},DEzone7{I,j});
            end
            x=1;
            for k=intersect(ALL_DlyCellF{i,2},ALL_DlyCellF{i,j})
                c=corr(ALL_timeseqF2{i,2}(k,zone)',ALL_timeseqF2{I,j}(k,zone)');
                C(x,l)=c;
                x=x+1;
            end
        end
        [p,~,stats]=anova1(C);
        xticklabels({'ESzone','OEzone','EEzone'})
        ylabel('correlation')
        axis square
        title2=sprintf('p = %d',p);
        title(title2)
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        title3=sprintf('corr VSzone(TM%d stage%d)',i,j);
            printFig(gcf,title3);
        figure;c=multcompare(stats);
        title3_1=sprintf('1vs2 %d \n 1vs3 %d \n 2vs3 %d',c(1,6),c(2,6),c(3,6));
        title(title3_1)
        yticklabels({'ESzone','OEzone','EEzone'})
        axis square
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        title4=sprintf('corr VSzone(TM%d stage%d) mult',i,j);
            printFig(gcf,title4);
    end
end



%%%全部
X=[];
for i=1:5
addpath(PathName{i})
load 20220318.mat OK_ens
X{i}=size(OK_ens,2);
end
x=0;
ALL_PlaceCellF=cell(1,4);ALL_PlaceCellFonTM=cell(1,4);
ALL_noDSEseqF=cell(1,4);ALL_PlaceCellFonF=cell(1,4);ALL_PseqF=cell(1,4);
for i=1:5
    addpath(PathName{i})
    load 20220318.mat PlaceCellF PlaceCellFonTM noDSEseqF PlaceCellFonF PseqF
    if i==1
    x=0;
    else
    x=x+X{i-1};
    end
    for j=1:4
        if ~isempty(PlaceCellF{1,j})
            ALL_PlaceCellF{1,j}=[ALL_PlaceCellF{1,j} PlaceCellF{1,j}+x];
        end
        if ~isempty(PlaceCellFonTM{1,j})
            ALL_PlaceCellFonTM{1,j}=[ALL_PlaceCellFonTM{1,j} PlaceCellFonTM{1,j}+x];
        end
        if ~isempty(PlaceCellFonF{1,j})
            ALL_PlaceCellFonF{1,j}=[ALL_PlaceCellFonF{1,j} PlaceCellFonF{1,j}+x];
        end
        if ~isempty(noDSEseqF{1,j})
            ALL_noDSEseqF{1,j}=vertcat(ALL_noDSEseqF{1,j},noDSEseqF{1,j});
        end
        if ~isempty(PseqF{1,j})
            ALL_PseqF{1,j}=vertcat(ALL_PseqF{1,j},PseqF{1,j});
        end
    end
end
save 20221206.mat -append ALL_PlaceCellF ALL_PlaceCellFonTM ALL_noDSEseqF ALL_PlaceCellFonF ALL_PseqF
options = statset('MaxIter',1000);
z=10;
T=159;
ADC=[];ATC=[];
for i=2:3
    for j=2:4
        ADC=[ADC ALL_DlyCellF{i,j}];
        ATC=[ATC ALL_TimeCellF{i,j}];
    end
end
ADC=union(ADC,ADC);
ATC=union(ATC,ATC);
Cells=setdiff(ALL_PlaceCellF{1,1},ATC);
Cells=intersect(ATC,ALL_PlaceCellF{1,1});
Cells=ALL_PlaceCellF{1,1};
extractD=0;
for h=1:4
test=[];
AIC = zeros(1,z);
GMModels = cell(1,z);
test=[];test1=[];
    if ~isempty(ALL_noDSEseqF{1,h})
        
        for i=Cells
%         for i=ALL_PlaceCellF{1,h}
        out=[];xy=[];
        MM=1;%
        ThF=30;
            if extractD
            [out,xy,peak_Point]=countPlaceFieldsM(ALL_noDSEseqF{1,h}(i,:),2.5,ThF,MM);
                if ~isempty(out)
                    for j=1:size(xy,2)
                        I=zeros(1,T);
                        I(1,xy{j}(:,1))=ALL_noDSEseqF{1,h}(i,xy{j}(:,1))/max(ALL_noDSEseqF{1,h}(i,xy{j}(:,1)));
                        test=[test xy{j}(:,1)' (xy{j}(:,1)+T)'];%2周分
                    end
                end
            else
            %遅延ありver.
            [out,xy,peak_Point]=countPlaceFieldsM(ALL_PseqF{1,h}(i,:),2.5,ThF,MM);
                if ~isempty(out)
                    for j=1:size(xy,2)
                        I=zeros(1,T);
                        I(1,xy{j}(:,1))=ALL_PseqF{1,h}(i,xy{j}(:,1))/max(ALL_PseqF{1,h}(i,xy{j}(:,1)));
                        test=[test xy{j}(:,1)' (xy{j}(:,1)+T)'];%2周分
                    end
                end
            end
        end
        for k=1:z
            GMModels{k} = fitgmdist(test',k,'Options',options,'CovarianceType','diagonal');
            AIC(k)= GMModels{k}.AIC;
        end
    [minAIC,numComponents] = min(AIC);
    y=numComponents;
    GMModel = fitgmdist(test',y);
    xgrid = linspace(min(test),max(test),T*2)';%2周分にしてseqの両端にfieldがたまっててもわかるようにする
    figure
    yyaxis left
    hold on;histogram(test,T)
    yyaxis right
    hold on; plot(xgrid,pdf(GMModel,xgrid),'k');
    M=[M max(pdf(GMModel,xgrid))];
    title(GMModel.mu)
    % title3=sprintf('TM%d stage%d',g-1,h);
    % xlabel('elapsed time (20ms/bin)')
    % yyaxis left
    % ylabel('field Num')
    % set(gca, 'FontName', 'Times New Roman','FontSize',16)
    % FH=gca;
    % printFig(gcf,title3);
    end
end

test=[];
for i=ALL_PlaceCellF{1,1}
out=[];xy=[];
MM=1;ThF=50;%20220608変更10→50
[out,xy,peak_Point]=countPlaceFieldsM(ALL_PseqF{1,1}(i,:),1,ThF,MM);
T=size(ALL_PseqF{1,1}(i,:),2);
if ~isempty(out)
for j=1:size(xy,2)
I=zeros(1,T);
I(1,xy{j}(:,1))=ALL_PseqF{1,1}(i,xy{j}(:,1))/max(ALL_PseqF{1,1}(i,xy{j}(:,1)));
test=[test xy{j}(:,1)'];
end
end
end
figure;histogram(test,'BinWidth',10)
ylabel('field Num')
xticks([0:10:size(ALL_PseqF{1,1},2)])

h=1;
[N_seq,Norder,NI]=Nseq(ALL_PseqF{1,h},ALL_PlaceCellF{1,1});
for h=1:4
    [N_seq]=Nseq(ALL_PseqF{1,h},ALL_PlaceCellF{1,1});
    figure;imagesc(N_seq(Norder,:))
    title3=sprintf('stage%d',h);
    colorbar;
    xlabel('Linearized location(2.5cm/bin)')
    ylabel({'number of cells';'(sorted by stage 1 trial)'})
    set(gca, 'FontName', 'Times New Roman','FontSize',16)
    printFig(gcf,title3);
end 

% place cell dot plot stage 1 vs stages
ALL_dot_VSstage1=cell(1,4);
for j=2:4
    XX=length(ALL_PseqF{1,j});
    title3=sprintf('map dot plot (stage 1 vs %d)',j);
    title3_1=sprintf('stage 1 vs stage %d',j);
    jitter=0;
    XXX=seq_dot(ALL_PseqF{1,1},ALL_PseqF{1,j},ALL_PlaceCellF{1,1},2,4,jitter);
    ALL_dot_VSstage1{1,j}=XXX;
    figure;plot(XXX(:,1),XXX(:,2),'k.')
    FH=gca;
    title(title3_1)
%     close
    printFig(gcf,title3);
    % density color map ver.
    % 15cm/bin SD=0.8
    D=[];BTh=6;
    for k1=1:ceil(XX/BTh)
        for k2=1:ceil(XX/BTh)
            Y=[];
            Y=find((BTh*(k1-1)+.1)<XXX(:,1) & XXX(:,1)<BTh*k1+.1 & (BTh*(k2-1)+.1)<XXX(:,2) & XXX(:,2)<BTh*k2+.1);
            D(k1,k2)=size(Y,1)/size(XXX,1);
        end
    end
    D=D*100;
    SD=0.8;
    D=imgaussfilt(D,SD);
    FH=figure;imagesc(D');colorbar
    axis equal
    xlim([1 size(D,1)])
    ylim([1 size(D,2)])
    axis xy
%     xticks([0:J:size(D,1)])
%     yticks([0:J:size(D,2)])
    title4=sprintf('%s D',title3);
    xlabel('Linearized location(15cm/bin)')
    ylabel('Linearized location(15cm/bin)')
    title(title3_1)
    set(gca, 'FontName', 'Times New Roman','FontSize',16)
    printFig(FH,title4);
end
%place cell histogram & gaussfit
options = statset('MaxIter',1000);
z=4;
for g=1:4
test=[];
AIC = zeros(1,z);
GMModels = cell(1,z);
test=[];
for i=ALL_PlaceCellF{1,1}
out=[];xy=[];
MM=1;ThF=4;
[out,xy,peak_Point]=countPlaceFieldsM(ALL_PseqF{1,g}(i,:),1,ThF,MM);
T=size(ALL_PseqF{1,g}(i,:),2);
if ~isempty(out)
for j=1:size(xy,2)
I=zeros(1,T);
I(1,xy{j}(:,1))=ALL_PseqF{1,g}(i,xy{j}(:,1))/max(ALL_PseqF{1,g}(i,xy{j}(:,1)));
test=[test xy{j}(:,1)'];
end
end
end
for k=1:z
GMModels{k} = fitgmdist(test',k,'Options',options,'CovarianceType','diagonal');
AIC(k)= GMModels{k}.AIC;
end
[minAIC,numComponents] = min(AIC);
y=numComponents;
GMModel = fitgmdist(test',y);
xgrid = linspace(min(test),max(test),x)';
figure
yyaxis left
hold on;histogram(test,'BinWidth',5)
yyaxis right
hold on; plot(xgrid,pdf(GMModel,xgrid),'k');
title3=sprintf('stage%d',g);
title(title3)
xlabel('Linearized location(2.5cm/bin)')
yyaxis left
ylabel('field Num')
xticks([0:20:size(ALL_PseqF{1,g},2)])
set(gca, 'FontName', 'Times New Roman','FontSize',16)
FH=gca;
printFig(gcf,title3);
end



% pmap field sequence
for h=1:4
    cells=ALL_PlaceCellF{1,h};
    x=1;Seq=[];
    for i=cells
        sequence=ALL_PseqF{1,h}(i,:);
        Seq(x,:)=zeros(size(sequence));
        MM=2;ThF=30;
        [out,xy,peak_Point]=countPlaceFieldsM(sequence,2.5,ThF,MM);
        for j=1:size(xy,2)
            Seq(x,xy{j}(:,1))=sequence(:,xy{j}(:,1));
            x=x+1;
        end
    end
    [N_seq,Norder]=Nseq(Seq);
    figure;imagesc(N_seq(Norder,:))
end
%time field sequence
for i=2:3
    for j=2:4
        if j==3
            XX=1250;
        else
            XX=750;
        end
        cells=ALL_DlyCellF{i,2};
        x=1;Seq=[];
        for k=cells
            sequence=ALL_timeseqF2{i,j}(k,1:XX);
            Seq(x,:)=zeros(size(sequence));
            MM=1;ThF=50;
            [out,xy,peak_Point]=countPlaceFieldsM(sequence,1,ThF,MM);
            for l=1:size(xy,2)
                Seq(x,xy{l}(:,1))=sequence(:,xy{l}(:,1));
                x=x+1;
            end
        end
        [N_seq,Norder]=Nseq(Seq);
        figure;imagesc(N_seq(Norder,:))
    end
end

%%%全部　ALL_sp 20ms/1bin 150cells×sp(30s/40s=1499/1999bin)
% TM前後のrunning speed
% X:loadの時できたやつ
ALL_sp=cell(2,3);
for i=1:5
addpath(PathName{i})
load 20220318.mat All_sp
for k=1:2
    for l=1:3
        for j=1:X{i}
            ALL_sp{k,l}=vertcat(ALL_sp{k,l},All_sp{k,l});
        end
    end
end
clear All_sp
end
% TM速度が52.83cm/sec
for i=1:2%TM
    for j=1:3%stage
        figure
        M=[];
        for k=1:5
            if k==1
            x=0;
            else
            x=x+X{k-1};
            end
        SP=ALL_sp{i,j}(x+1,:)*(52.38/4.3109);
        hold on;plot(SP,'LineWidth',2)
        M=[M max(SP)];
        end
        M=max(M);
        c=[0.4660 0.6740 0.1880];
        XX=[DSzone7{i+1,j+1}(1) DSzone7{i+1,j+1}(end) DSzone7{i+1,j+1}(end) DSzone7{i+1,j+1}(1)];
        YY=[0 0 M*1.1 M*1.1];
        hold on;patch(XX+475,YY,c,'FaceAlpha',.5,'EdgeColor','none')
        c=[0.3660 0.4740 0.2880];
        XX=[DEzone7{i+1,j+1}(1) DEzone7{i+1,j+1}(end) DEzone7{i+1,j+1}(end) DEzone7{i+1,j+1}(1)];
        YY=[0 0 M*1.1 M*1.1];
        hold on;patch(XX+475,YY,c,'FaceAlpha',.5,'EdgeColor','none')
        xlabel('elapsed time (20ms)')
        ylabel('running speed (cm/sec)')
        xticks([0:250:length(ALL_sp{i,j})])
        legend('rat 1','rat 2','rat 3','rat 4','rat 5','ES zone','EE zone')
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        title3=sprintf('TM%d stage%d',i+1,j+1);
        printFig(gcf,title3);
    end
end
% 加速度ver.
% figure出す必要なし
% 5rats*6ocnsition分のmax(diff(SP))とそれ以外のzoneのmax(diff(SP))をsigned runk test
% ES/EE zoneの2個分　boxplot？のfigure出す
AS=[];AE=[];otherAS=[];otherAE=[];
for i=1:2%TM
    for j=1:3%stage
%         M=[];
        for k=1:5
            if k==1
            x=0;
            else
            x=x+X{k-1};
            end
        SP=ALL_sp{i,j}(x+1,:)*(52.38/4.3109);
%         hold on;plot(diff(SP),'LineWidth',2)
%         M=[M max(diff(SP))];
        zone=DSzone7{i+1,j+1}+475;%spのsizeと合わせる
        AS=[AS max(diff(SP(:,zone)))];
        L=length(zone);
        otherAS=[otherAS max(diff(SP(:,(zone-L):zone-1)))];
        zone=DEzone7{i+1,j+1}+475;
        AE=[AE max(diff(SP(:,zone)))];
        L=length(zone);
        otherAE=[otherAE max(diff(SP(:,(zone(end)+1):(zone(end)+L))))];
        end
    end
end
% 一応減速度ver.も有意差があった
% DS=[];DE=[];otherDS=[];otherDE=[];
% for i=1:2%TM
%     for j=1:3%stage
% %         M=[];
%         for k=1:5
%             if k==1
%             x=0;
%             else
%             x=x+X{k-1};
%             end
%         SP=ALL_sp{i,j}(x+1,:)*(52.38/4.3109);
% %         hold on;plot(diff(SP),'LineWidth',2)
% %         M=[M max(diff(SP))];
%         zone=DSzone7{i+1,j+1}+475;%spのsizeと合わせる
%         DS=[DS min(diff(SP(:,zone)))];
%         L=length(zone);
%         otherDS=[otherDS min(diff(SP(:,(zone-L):zone-1)))];
%         zone=DEzone7{i+1,j+1}+475;
%         DE=[DE min(diff(SP(:,zone)))];
%         L=length(zone);
%         otherDE=[otherDE min(diff(SP(:,(zone(end)+1):(zone(end)+L))))];
%         end
%     end
% end
% figure;boxplot([AS;otherAS]')
% p = signrank(AS,otherAS);
% xticklabels({'ESzone','other zone'})
% title1=sprintf('ES zone \n max rate of acceleration (cm/sec)');
% ylabel(title1)
% title2=sprintf('p = %d',p);
% title(title2)
% axis square
% set(gca, 'FontName', 'Times New Roman','FontSize',16)
% title3=sprintf('rate of acceleration (ES zone)');
% printFig(gcf,title3);
figure;boxplot([AE;otherAE]')
p = signrank(AE,otherAE);
xticklabels({'EEzone','other zone'})
title1=sprintf('EE zone \n max rate of acceleration (cm/sec)');
ylabel(title1)
title2=sprintf('p = %d',p);
title(title2)
axis square
set(gca, 'FontName', 'Times New Roman','FontSize',16)
title3=sprintf('rate of acceleration (EE zone)');
% printFig(gcf,title3);
% 相関関係とshuffledの比較出す
% 20230712 zoneでなくpointに変更　詳細はsp_seq　line9参照
for i=1:2%TM
    for j=1:3%stage
        for k=1:2
            if k==1
                zone=DSzone7{i+1,j+1};
                cells=ALL_DSCellF{i+1,j+1};
                seq=ALL_DSseqF{i+1,j+1};
                seq_ori=ALL_DSseqF_ori{i+1,j+1};
                title3='ES zone';
            elseif k==2
                zone=DEzone7{i+1,j+1};
                cells=ALL_DECellF{i+1,j+1};
                seq=ALL_DEseqF{i+1,j+1};
                seq_ori=ALL_DEseqF_ori{i+1,j+1};
                title3='EE zone';
            end
        title1=sprintf('FR(TM%d stage%d)',i,j+1);
        title2=sprintf('speed(TM%d stage%d)',i,j+1);
        sp_seq(cells,ALL_sp{i,j},seq,seq_ori,zone,1,title1,title2,title3);
        Title=sprintf('%s VS %s (%s)',title1,title2,title3);
        FH=gcf;
        xlabel('elapsed time (20ms/bin)')
        ylabel(sprintf('%s \n number of cells',title3))
        set(gca, 'FontName', 'Times New Roman','FontSize',16)
        printFig(FH,Title);
        end
    end
end
%全条件のまとめて出す→speed6個重ねたものにcellを合体したもの 合計89*6になれば正解
% cell_Num→足したものを後ろに追加していく　sp→ひたすら合体　seq,seq_ori→ひたすら合体
% point→各条件のcell_Num分だけ同じものを追加していく
% shuffleの分布と実測値の分布が同じかどうかの検定になるので5%の線はひかなくてもいいのか
SSS=cell(1,2);AAA=cell(1,2);
for k=1:2%DS/E point
    ALL_point7=[];ALL_SPEED=zeros(89*6,2451);ALL_CELLS=[];ALL_seq=cell(1,1);ALL_seq_ori=cell(1,1);x=0;
    for i=1:2%TM
        for j=1:3%stage
            if k==1
                point=DSpoint7{i+1,j+1};
                cells=ALL_DSCellF{i+1,j+1};
                seq=ALL_DSseqF{i+1,j+1};
                seq_ori=ALL_DSseqF_ori{i+1,j+1};
                title3='ES zone';
            elseif k==2
                point=DEpoint7{i+1,j+1};
                cells=ALL_DECellF{i+1,j+1};
                seq=ALL_DEseqF{i+1,j+1};
                seq_ori=ALL_DEseqF_ori{i+1,j+1};
                title3='EE zone';
            end
            if i==1 && j==1
                ALL_seq=seq;
                ALL_seq_ori=seq_ori;
            else
               ALL_seq=vertcat(ALL_seq,seq);
               ALL_seq_ori=vertcat(ALL_seq_ori,seq_ori);
            end
            ALL_SPEED((1+x*89):(x*89+89),1:length(ALL_sp{i,j}))=ALL_sp{i,j};
            ALL_CELLS=[ALL_CELLS cells+89*x];
           for l=cells
            ALL_point7=[ALL_point7 point];
           end
           x=x+1;
        end
    end
    [~,~,~,TESTS,AA]=sp_seq(ALL_CELLS,ALL_SPEED,ALL_seq,ALL_seq_ori,ALL_point7,0);
    Title=sprintf('ALL_FR vs ALL_SP (%s)',title3);
    FH=gcf;
    xlabel('elapsed time (20ms/bin)')
    ylabel(sprintf('%s \n number of cells',title3))
    set(gca, 'FontName', 'Times New Roman','FontSize',16)
    printFig(FH,Title);
    
    SSS{k}=TESTS;AAA{k}=AA;
    M1=histcounts(SSS{k},'BinWidth',.1);
    figure;histogram(SSS{k},'BinWidth',.1)
    yyaxis left
    ylabel('fraction of shuffles')
    ylim([0 max(M1*1.2)])
    C=sort(SSS{k});
    shN=length(SSS{k});
    Pval=shN*.025;
    c=[0.9290 0.6940 0.1250];
    XX=[C(Pval) C(shN-Pval) C(shN-Pval) C(Pval)];
    YY=[0 0 max(M1*1.2) max(M1*1.2)];
    hold on;patch(XX,YY,c,'FaceAlpha',.25,'EdgeColor','none')
    hold on;xline(C(Pval),'--','2.5%','LabelOrientation','horizontal','LabelHorizontalAlignment','center')
    hold on;xline(C(shN-Pval),'--','97.5%','LabelOrientation','horizontal','LabelHorizontalAlignment','center')
    yyaxis right
    hold on;histogram(AA,'BinWidth',.1)
    M2=histcounts(AA,'BinWidth',.1);
    ylabel('observed')
    ylim([0 max(M2*1.2)])
    xlabel('correlation between speed and activity')
    s=length(union(find(AA<C(Pval)),find(AA>C(shN-Pval))));
    p=.05;
    n=length(AA);
    pout=myBinomTest(s,n,p,'two');
    title(sprintf('p=%d',pout))
    axis square
    set(gca, 'FontName', 'Times New Roman','FontSize',16)
    title1=sprintf('corr(sp vs FR) (%s)',title3);
    printFig(gcf,title1); 
end