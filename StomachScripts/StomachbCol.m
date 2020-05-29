function Stats=StomachbCol(czifile,s,tGSII)

% this function is used to quantify the amount of Trop2 staining in gastric
% tissue sections labelled with DAPI and collagenVI, as well as unknown marker (GSII?).
% CollagenVI signal is found in interstitial tissue but excluded from gland
% tissues and it is therefore used after signal complementation to segment
% glands.

% INPUT:
% czifile  - name of the .czi file to analyze
% s        - factor used for background correction of the  collagen signal (range [0.5 2])
% tGSII    - threshold value ([0 1]) derived from Otsu used to segment GSII
% labelling

% OUTPUT:
% Stats    - structure array containg the quantified metrics
% Example:
% czifile='1569_DAPI_CollagenVI_Unknown_Trop2.czi';
% Stats=ValeriebCol(czifile,0.8,1)


fn=czifile(1:end-4);
expr='_';
ss=split(fn,expr);
im=ss{end};


%% import data
IS=bfopen(czifile);

% change the indices for channel correspondence
% e.g. if DAPI is in channel 2,
% C1=IS{1,1}{2,1};

C1=IS{1,1}{1,1}; %DAPI
C2=IS{1,1}{2,1}; %collagen
C3=IS{1,1}{3,1}; %GSII
C4=IS{1,1}{4,1}; %Trop2


%% DAPI segmentation

h1=fspecial('average',4);
T=graythresh(C1);
C1W=imbinarize(C1,0.8*T);
C1Wf=imfilter(C1W,h1);

%% Collagen-based gland segmentation
C2g=imgaussfilt(C2,2);
C2m=medfilt2(C2g);
C2mt=imtophat(C2m,strel('disk',25));


q=median(C2mt(:));
C2mtb=C2mt-s*q;


% find best threshold value for collagen by scanning threshold values leading to highest
% numbers of objects
t=1:10:max(C2mtb(:));
nob=zeros(numel(t),1);
for i=1:numel(t)
FWW=C2mtb>=t(i);
C=bwconncomp(FWW);
nob(i)=C.NumObjects;
end

NOM=movmean(nob,100);
[hp, pp]=findpeaks(NOM);
[~, ip]=max(hp);
TH=pp(ip);

if isempty(TH)
TH=1;
end

FWW=C2mtb>=TH;
FWWa=bwareaopen(FWW,500);
h2=fspecial('average',[5 2]);
FWWaf=imfilter(FWWa,h2);

%closing the gaps of glands opened towards the lumen

FI=imfill(FWWaf,'holes');
hh=fspecial('average',10);
FII=imfilter(FI,hh);
FIIsk=bwskel(FII,'MinBranchLength',50);
Se=regionprops(bwmorph(FIIsk,'endpoints'),'Centroid');
Cs=vertcat(Se.Centroid);
p=pdist(Cs);
Z=squareform(p);
Z(logical(tril(ones(size(Z)))))=Inf;
[m, I]=min(Z,[],2);
np=1:numel(Se);
dval=[m I np'];
tf=m>200;
dval(tf,:)=[];
if ~isempty(dval)
dvals=sortrows(dval,1);
dvals(end,:)=[];
[~, Ia, ~]=unique(dvals(:,2));
dvalsu=dvals(Ia,:);
endpairs=dvalsu(:,2:3);
bw=false(size(FWWaf));
bw2=false(size(FWWaf));
figure, imshow(bw);
for i=1:size(endpairs,1)
    hline=drawline(gca,'Position',[Cs(endpairs(i,1),:); Cs(endpairs(i,2),:)]);
bw=createMask(hline); bw2=bw2 | bw;
end

bw2=bwmorph(bw2,'diag');
FWWaf=FWWaf | bw2;
close;
end
%

rFWWaf=imclearborder(~FWWaf);
rFo=imopen(rFWWaf,strel('disk',5));
rFoa=bwareaopen(rFo,1000);
rFoa=bwmorph(rFoa,'thicken',5);




%% GSII segmentation
%C3t=imtophat(C3,strel('disk',25));
%FRtm=medfilt2(FRt);
C3tm=medfilt2(C3);
%TF=graythresh(FRtm);

C3W=imbinarize(C3tm,0.5*tGSII);
C3Wa=bwareaopen(C3W,10);
C3Waf=imfill(C3Wa,'holes');

%% Trop2 segmentation
h4=fspecial('average',4);
C4t=imtophat(C4,strel('disk',10));
mC4t=mean(C4t(:));
stdC4t=std(double(C4t(:)));
T=mC4t+3*stdC4t;
C4W=C4t>T;
C4Wf=imfilter(C4W,h4);

%% Stats
Stats=regionprops(rFoa,'Area','Centroid','Solidity','Eccentricity','Circularity','MajorAxisLength','MinorAxisLength','PixelIdxList');
%SG=regionprops(rFoa,G,'MeanIntensity');
STrop2=regionprops(rFoa,C4,'MeanIntensity');
%SGW=regionprops(rFoa,GWf,'PixelValues');
STrop2W=regionprops(rFoa,C4Wf,'PixelValues');
SGSIIW=regionprops(rFoa,C3Waf,'PixelValues');

for i=1:numel(Stats)
    %Stats(i).PercG=sum(SGW(i).PixelValues)/Stats(i).Area;
    %Stats(i).SumG=sum(SGW(i).PixelValues);
    Stats(i).PercTrop2=sum(STrop2W(i).PixelValues)/Stats(i).Area; %fraction of gland area that is Trop2 positive
    %Stats(i).SumTrop2=sum(STrop2W(i).PixelValues);
    Stats(i).PercGSII=sum(SGSIIW(i).PixelValues)/Stats(i).Area; %fraction of the gland area that is GSII positive
end

%[Stats(1:numel(Stats)).MIG488]=deal(SG.MeanIntensity);
[Stats(1:numel(Stats)).MITrop2]=deal(STrop2.MeanIntensity); %mean intensity of the trop2 signal

% filter out small and weirdly shaped glands
A=[Stats.Area]';
S=[Stats.Solidity]';
tf=A<2000 | S<0.65;
Stats(tf)=[];


MASK=false(size(C1));
Pix=vertcat(Stats.PixelIdxList);
MASK(Pix)=true;


Stats=struct2table(Stats);

if ~isempty(Stats)
IMG=[repmat(cellstr(im),size(Stats,1),1) cellstr(int2str((1:size(Stats,1))'))];
IMG=table(IMG,'VariableNames',{'SampleID'});
Stats=[IMG Stats];
Stats.PixelIdxList=[];
else
    Stats=[];
end

%C1W=imbinarize(C1);
C1C3=C1Wf & C3Waf;
nD=numel(find(C1Wf));

% m=mean(C3(C1Wf));
% s=std(double(C3(C1Wf)));
% T=m+3*s;
% h=fspecial('average',4);
% C3W=C3>T;
% C3Wf=imfilter(C3W,h);
% rC3I=mean(C3(C3Wf))/nD;
% rnC3=numel(find(C1C3))/nD;
% C3S=[rC3I rnC3];
% C3S=rnC3;

%% Figures
%figure, imshow(imadjust(FRmtb));hold on
figure, imshow(C4Wf);hold on
if ~isempty(Stats)
B=bwboundaries(MASK);

for i=1:length(B)
    bb=B{i};
    plot(gca,bb(:,2),bb(:,1),'-r','Linewidth',1);
end

CC=[Stats.Centroid];

for k=1:size(Stats,1)
    text(CC(k,1),CC(k,2)-15,['#' int2str(k)],'color','w','FontWeight','Bold');
    text(CC(k,1),CC(k,2)+15,num2str(Stats.PercTrop2(k),'%0.3f'),'color','y');
end
end

% 
lrgb=label2rgb(bwlabel(rFoa),@jet,'k','shuffle');
hi=imshow(lrgb);
set(hi,'AlphaData',0.3);
drawnow;

% 
% figure, plot(NOM);hold on;plot(gca,nob);
%  figure, imshow(C4,[]);
%  figure, imshow(C3,[]);
%   figure, imshowpair(C1Wf,C1C3);
%   text(100,100,num2str(C3S,'%0.3f'),'color','y','FontWeight','Bold','FontSize',16);
%   figure, imshow(C1Wf);
% figure, imshow(rFoa);
