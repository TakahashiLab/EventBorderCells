%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%get times
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [riseT,fallT,DS]=gettimesM(x,Th)

P=find(x>Th);

%if event width==1
P(x(P+1)<Th)=[];

P=[0 P length(x)];


tmp=find((diff(P)~=1));
tmp=unique([tmp tmp+1]);
tmp(1)=[];

PosT=P(tmp);

%remove oversampled data(interval <5msec)
ro=find(diff(PosT)<125);%125:5msec 25000:1sec
ro=unique([ro ro+1]);

 
if size(ro,2)>1
  fprintf('removing oversampled data:%d...\n',size(ro,2));
end

if x(1) > Th
    riseT=PosT(2:2:end-1);
    fallT=PosT(3:2:end);
else
    riseT=PosT(1:2:end-1);
    fallT=PosT(2:2:end);
end

% dly時間が一秒以上変化した最初のriseTを検索
DS=[];
f=fallT(1)-riseT(1);
for i=2:size(riseT,2)
g=fallT(i)-riseT(i);
if abs(f-g)>2500
DS=[DS riseT(i)];
end
f=fallT(i)-riseT(i);
end

return;

