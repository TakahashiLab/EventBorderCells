% 入力　seqT,seqD,seqP
% 出力　T=1,D=2,P=3
% glmは長いseqになるほど残差も大きくなるので一番小さいseq sizeに合わせて他を圧縮して比較

function [kekka]=GLMm2(seqT,seqD,seqP)
    kekka=[];
    warning ('off')

    m=min([size(seqT,2) size(seqD,2) size(seqP,2)]);
    if m==size(seqT,2)
        seqD=compress_Seq(seqT,seqD);
        seqP=compress_Seq(seqT,seqP);
    elseif m==size(seqD,2)
        seqT=compress_Seq(seqD,seqT);
        seqP=compress_Seq(seqD,seqP);
    elseif m==size(seqP,2)
        seqT=compress_Seq(seqP,seqT);
        seqD=compress_Seq(seqP,seqD);
    end
    seqT=seqT/max(seqT);
    seqD=seqD/max(seqD);
    seqP=seqP/max(seqP);
    
%     figure
%     hold on;plot(seqT)
%     hold on;plot(seqD)
%     hold on;plot(seqP)
%     legend('seqT','seqD','seqP')
    
    DT=1:size(seqT,2);
    mdlT=fitglm(DT',seqT','poly5','Distribution','poisson');
    Dev(1)=mdlT.Deviance;
    
    DD=1:size(seqD,2);
    mdlD=fitglm(DD',seqD','poly5','Distribution','poisson');
    Dev(2)=mdlD.Deviance;
    
    DP=1:size(seqP,2);
    mdlP=fitglm(DP',seqP','poly5','Distribution','poisson');
    Dev(3)=mdlP.Deviance;
    
    
    [M]=min(Dev,[],2);
    if intersect(1,find(Dev==M))
        kekka=[kekka 1];
    end
    if intersect(2,find(Dev==M))
        kekka=[kekka 2];
    end
    if intersect(3,find(Dev==M))
        kekka=[kekka 3];
    end
%     legend('seqT','seqD',(sprintf('%d',kekka)))
return;

%%%%%%%%%%%%%%%%%%%%%%%%
% seq1に合わせてseq2を圧縮　seq2_1として出力
function [seq2_1]=compress_Seq(seq1,seq2)
seq2_1=[];
b=size(seq2,2);c=size(seq1,2);
a=b/c;
for i=1:c
if a*i<b
seq2_1(1,i)=(sum(seq2(1,a*(i-1)+1:a*i)))/a;
end
end
