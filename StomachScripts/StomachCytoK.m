function Stats=StomachCytoK(pd,k)

%This function is used to to quantify the number of Ki67 nuclei and the amount of GSII labelling in gland
%tissues from gastric IF sections labelled with DAPI, GSII, Cytokeratin and a 4th marker

%INPUT:
% pd  -a structure array returned by the "dir" function that contains the
% list of all the .lsm files to analyze
% k   -the index of the file/image to analyze

%OUTPUT:
% Stats  -a structure array containing the quantified metrics in each
% field.

%% read the .lsm file - requires bioformats readers
IS=bfopen(pd(k).name);

%Switch the indices depending on the marker's channel, i.e. if DAPI is in
%channel 2, use D=IS{1,1}{2,1} instead.

D=IS{1,1}{1,1}; %DAPI

G=IS{1,1}{4,1}; %CytoK

R=IS{1,1}{3,1}; %Ki67

FR=IS{1,1}{2,1}; %GSII

Flag=1;
h=fspecial('average',4);

%% DAPI Segmentation
T=graythresh(D);
DW=imbinarize(D,0.8*T);
DWf=imfilter(DW,h);
if numel(find(DWf))<numel(D)/4
    Flag=0;
end


%% GSII segmentation
FRt=imtophat(FR,strel('disk',25));
FRtm=medfilt2(FRt);
TF=graythresh(FRtm);

FRW=imbinarize(FRtm,0.75*TF);
FRWa=bwareaopen(FRW,10);
FRWaf=imfill(FRWa,'holes');

%% CytoK Segmentation
GW=imbinarize(G);

GWc=imclose(GW,strel('disk',3));
GWf=imfilter(GWc,h);
IGW=~GWf;
IGWa=bwareafilt(IGW,[1 100]);
GWf=GWf | IGWa;
GWfo=imopen(GWf,strel('disk',2));

%% Ki67 segmentation
Rm=medfilt2(R);

RmWf=imbinarize(Rm);
RmWf=bwareaopen(RmWf,10);
RmWfw=Nuclear2(RmWf);
RmWfw=bwareaopen(RmWfw,3);
KiCover=numel(find(RmWf(:)))/numel(find(DWf(:)));

if KiCover>0.2
    Flag=0;
    tb=minminT(Rm);
    Rmb=Rm-1.5*tb;
    RmWf=imbinarize(Rmb);
    RmWf=bwareaopen(RmWf,10);
    RmWfw=Nuclear2(RmWf);
    RmWfw=bwareaopen(RmWfw,3);
    KiCover=numel(find(RmWf(:)))/numel(find(DWf(:)));
end



StatsR=regionprops(RmWfw,'Centroid');
%RC=bwconncomp(RmWf);
nbR=numel(StatsR);





%% Compute metrics

GRW=FRWaf | GWfo; %merge GSII and Cytokeratin masks to get glands

DGRW=GRW & DWf; %retain gland tissue overlapping with DAPI
RmWf2=bwmorph(RmWfw,'thicken',2);
Marker=RmWfw & DGRW;
RmWfR=imreconstruct(Marker,RmWf2); %gather Ki67 nuclei in DAPI positive 
% glands
RC=bwconncomp(RmWfR);

KiGland=RmWfw & GRW;
StatsRglands=regionprops(KiGland,'Centroid');
nbKiG=numel(StatsRglands);




GFrac=numel(find(FRWaf))/numel(find(GRW)); %fraction of GSII per gland tissue
GInt=mean(FR(FRWaf)); % mean intensity of GSII signal
%KI=nbR/numel(find(GRW));
KIG=RC.NumObjects/numel(find(DGRW)); % number of Ki67 nuclei per DAPI positive glands

Stats.GSCov=GFrac;
Stats.GSmInt=GInt;
Stats.NbKITotal=nbR; %total number of Ki67 nuclei
Stats.NbKIGlands=nbKiG; %number of Ki67 nuclei in glands
%Stats.KITotalCov=KI;
%Stats.KIGlandsCov=KIG;
Stats.KIDAPI=nbR/numel(find(DWf)); %relative number of Ki67 nuclei per DAPI area
Stats.KIDAPIGlands=KIG;
Stats.Flag=Flag;


%figure, imshow(R,[])
%figure, imshow(RmWf,[]);figure, imshow(FRWaf);figure, imshow(DWf);
B=bwboundaries(FRWaf);
% 
figure, imshow(GRW);hold on;
for i=1:length(B)
    bb=B{i};
    plot(gca,bb(:,2),bb(:,1),'-','Color',[0.8510    0.3255    0.0980],'Linewidth',3);
end
CC=[StatsR.Centroid];
scatter(gca,CC(1:2:end),CC(2:2:end),'o','MarkerFaceColor',[0.4667    0.6745    0.1882],'MarkerEdgeColor','y');

t1st=[int2str(nbR) ' (' int2str(nbKiG) ')'];
t1=text(100,50,t1st,'Color',[0.4667    0.6745    0.1882],'FontSize',24,'FontWeight','Bold');
t2=text(100,125,num2str(GFrac,1),'Color',[0.8510    0.3255    0.0980],'FontSize',24,'FontWeight','Bold');
t1.BackgroundColor='w';
t2.BackgroundColor='w';
    

end



function thresvalue=minminT(I)

% This subfunction defines a background value for background subtraction of
% image I

thresvalue = max([min(max(I,[],1)) min(max(I,[],2))])	 ;
end


function [NUC, Bn, Ln]= Nuclear2(D,varargin)


% This subfunction segments nuclei using watershed algorithm
% after distance transform of binary image D

Dl=bwdist(~D);
if ~isempty(varargin)
    s=varargin{1};
h=fspecial('disk',s); 
Dl=imfilter(Dl,h);
end
Dl=-Dl;
Dl(~D)=-Inf; % -Inf
%NUC=watershed(Dl,8);
NUC=watershed(Dl);
NUC(~D)=0;
NUC=logical(NUC);





[Bn, Ln, Nn]=bwboundaries(NUC,4,'noholes');
end

