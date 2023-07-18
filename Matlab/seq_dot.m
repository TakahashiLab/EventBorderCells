% seq比較をField dot plot(%発火率上位３位までを抽出)で クラスター解析もする
% cell:intersect(ALL_TimeCellC1{i,j1},ALL_TimeCellC1{i,j2})
% flagに1を入れるとき：クラスターの中心だけ出したいとき 2を入れるとき：seqの開始線がいらない時
% title3:figureのtitle title3=sprintf('stage2vs3 TM%d',i-1);
% 出力 X:dotの座標
function [X,C]=seq_dot(seq1,seq2,cell_num,flag,ThF,jitter)
figFlag=0;
if figFlag
    if flag==1
        c={'rx','gx','bx','kx','yx','mx','cx'};
    else
        c={'r.','g.','b.','k.','y.','m.','c.'};
    end
end
if nargin<5
    ThF=50;
else
end

X=[];xx=1;
for k=cell_num
    MM=2;%place cellの時にはMM=1
    [out1,~,peak_Point1]=countPlaceFieldsM(seq1(k,:),1,ThF,MM);%line16 peakRate*3/5
    [out2,~,peak_Point2]=countPlaceFieldsM(seq2(k,:),1,ThF,MM);
    if ~isempty(out1) && ~isempty(out2)%fieldがある
        if length(out1)>2%fieldが２つ以上あるとき
        [~,I]=sort(seq1(k,peak_Point1));
        peak_Point1=[peak_Point1(find(I==1)) peak_Point1(find(I==2)) peak_Point1(find(I==3))];
        end
        if length(out2)>2
        [~,I]=sort(seq2(k,peak_Point2));
        peak_Point2=[peak_Point2(find(I==1)) peak_Point2(find(I==2)) peak_Point2(find(I==3))];
        end
        Out1=[];Out2=[];
        for l=1:length(peak_Point1)
        if out1(l)<500%field sizeが10s以内
        Out1=[Out1 peak_Point1(l)];
        end
        end
        for l=1:length(peak_Point2)
        if out2(l)<500
        Out2=[Out2 peak_Point2(l)];
        end
        end
        if ~isempty(Out1) && ~isempty(Out2)
        X=vertcat(X,combvec(Out1,Out2)');
        xx=xx+size(combvec(Out1,Out2),2);%これはいるのか？
        end
    end
end

xy=[];
for y=1:size(X,1)
[~,~,sumd] = kmeans(X,y);
% hold on;plot(y,sum(sumd),'b.')
xy=[xy sum(sumd)];
end
A=[];B=[];C=[];

% J=250-jitter
for l=1:size(X,1)-4%多分エルボー法になってるはず
A=xy(l)-xy(l+1);B=xy(l+2)-xy(l+3);
if A>B*2
% hold on;plot(l,xy(l),'c*')
else
break
end
end
y=l;
% [idx,C,sumd] = kmeans(X,y);mieno
if figFlag
figure%
if flag==2
else
cc= [1.0000    0.8500    0.8500];
hold on;xline(jitter,'Color',cc,'LineWidth',5)
hold on;yline(jitter,'Color',cc,'LineWidth',5)
cc= [0.8500 0.3250 0.0980];
hold on;xline(size(seq1,2)-jitter,'Color',cc,'LineWidth',5)
hold on;yline(size(seq2,2)-jitter,'Color',cc,'LineWidth',5)

%%% DS vs DE ver.
% cc= [1.0000    0.8500    0.8500];
% hold on;plot(1:size(seq1,2),size(seq1,2)/2,'MarkerEdgeColor',cc,'Marker','.')
% cc=[0.9500 0.4250 0.2980];
% hold on;plot(size(seq2,2)/2,1:size(seq2,2),'MarkerEdgeColor',cc,'Marker','.')
end
if flag==1
    for m=1:y
    hold on;plot(C(m,1),C(m,2),c{m},...
        'MarkerSize',15,'LineWidth',3)
    end
else
    for m=1:y
        hold on;plot(X(idx==m,1),X(idx==m,2),'k.','MarkerSize',12)%
%     hold on;plot(X(idx==m,1),X(idx==m,2),c{m},'MarkerSize',12)%
    % hold on;plot(C(m,1),C(m,2),'kx',...%
    %     'MarkerSize',15,'LineWidth',3)
    % hold on;plot(C(m,1),C(m,2),c{xx},...
    end
end
axis equal
xlim([0 size(seq1,2)])
ylim([0 size(seq2,2)])
% xticks([0:20:size(seq1,2)])
% yticks([0:20:size(seq2,2)])
xticks([0:jitter:size(seq1,2)])%time seqの場合
yticks([0:jitter:size(seq2,2)])
xlabel('Elapsed time (s)')
ylabel('Elapsed time (s)')
set(gca, 'FontName', 'Times New Roman','FontSize',16)
% axis square
end

return;

% % 使用例
% for i=2:3
%     X=[];
%     for j=1:3 
%         I=i;X1=750;X2=750;
%         if j==1
%             j1=2;j2=3;
%             title3=sprintf('stage2vs3 TM%d',i-1);
%             
%         elseif j==2
%             j1=4;j2=3;
%             title3=sprintf('stage3vs4 TM%d',i-1);
%             
%         elseif j==3
%             j1=2;j2=4;
%             title3=sprintf('stage2vs4 TM%d',i-1);
%             
%         end
%         if j2==3
%             I=I+2;X2=1250;
%         end
%         x=1;
%         seq_dot(ALL_timeseqD2{i,j1}(:,1:X1),ALL_timeseqD2{I,j2}(:,1:X2),intersect(ALL_TimeCellD1{i,j1},ALL_TimeCellD1{i,j2}),0);
%         FH=gca;
%         printFig(gcf,title3);
%     end
% end