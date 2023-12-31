function info=calcInfoM(seq,oc_map)

if nargin==1
    [x,y]=size(seq);
    if x > y
        seq=seq';
    end
    
  PlaceMap=seq(1,:);
  oc_map=ones(size(PlaceMap));
end


for i=1:size(seq,1)
  PlaceMap=seq(i,:);
  [InfoPerSec,InfoPerSpk,InfoSparse]=calcInfoNew(PlaceMap,oc_map);
  info(i)=InfoSparse;% 変更　mieno
end
return;
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [InfoPerSec,InfoPerSpk,InfoSparse]=calcInfoNew(PlaceMap,OccupancyMap) 

tc=PlaceMap(:);
occ=OccupancyMap(:);
ix=tc~=0;
tc=tc(ix);
occ=occ(ix);

occ=occ/sum(occ);


f=sum(occ.*tc);
tc1=tc;
tc=tc/f;

if 0
SB=(occ(ix).*tc(ix)).*log2(tc(ix));
SB=sum(SB);

SS=(occ(ix).*tc1(ix)).*log2(tc(ix));
SS=sum(SS);

tau2=sum(tc1(ix).*occ(ix))^2;
SP=tau2/sum(occ(ix).*(tc1(ix).^2));
SP=sum(SP);
else
SB=(occ.*tc).*log2(tc);
SB=sum(SB);

SS=(occ.*tc1).*log2(tc);
SS=sum(SS);

tau2=sum(tc1.*occ)^2;
SP=tau2/sum(occ.*(tc1.^2));
SP=sum(SP);    
end


InfoPerSec=SS;
InfoPerSpk=SB;
InfoSparse=SP;



return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [InfoPerSec,InfoPerSpk,InfoSparse]=calcInfoC(PlaceMap,OccupancyMap) 

tc=PlaceMap(:);
occ=OccupancyMap(:);
occ=occ/sum(occ);
f=sum(occ.*tc);
tc1=tc;
tc=tc/f;
ix=tc~=0;
SB=(occ(ix).*tc(ix)).*log2(tc(ix));
SB=sum(SB);

tc1=tc;
SS=(occ(ix).*tc1(ix)).*log2(tc(ix));
SS=sum(SS);

tau2=sum(tc1(ix).*occ(ix))^2;
SP=tau2/sum(occ(ix).*(tc1(ix).^2));
SP=sum(SP);

%%%%%%%%%%%%%%%%
taux=PlaceMap;
id=find(taux~=0 & ~isnan(taux));
taux=taux(id);
OccupancyMap=OccupancyMap(id);

px=OccupancyMap/sum(OccupancyMap);
%tau=sum(taux.*px);
tau=mean(taux);
%px

%meaning taui=PlaceMap;

if tau==0
  InfoPerSec=NaN;
  InfoPerSpk=NaN;
  InfoSparse=NaN;
  return;
end

%Info=(PlaceMap/tau);
tau2=sum(taux.*px)^2;

%tmp=taux.*log2(taux/tau).*px;
%surf(reshape(tmp,64,64)');

%InfoPerSec=sum(taux.*log2(taux/tau).*px);
%InfoPerSpk=sum((taux/tau).*log2(taux/tau).*px);
%InfoSparse=tau2/sum(px.*(taux.^2));
InfoPerSec=SS;
InfoPerSpk=SB;
InfoSparse=SP;

%taux
%(taux/tau)
%log2(taux/tau)

return;



