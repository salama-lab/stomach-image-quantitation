function Stats=StomachCytoKnoK(pd,k)

%This function is used to quantify the amount of Ki67 nuclei
%from gastric IF sections labelled with DAPI, CytoK, Ki67 and GSII when
%CytoK signal is dim or absent

%INPUT:
% pd  - a structure array returned by the "dir" function that contains the
% list of all the .lsm files to analyze
% k   - the index of the file/image to analyze

%OUTPUT:
% Stats  - a structure array containing the quantified metrics in each
% field.
% Example:
% Stats=StomachCytoKnoK(dirczi,2);

%% read the .lsm file - requires bioformats readers

IS=bfopen(pd(k).name);

% IS=bfopen('1547_corpus_10x_03.czi');

%Switch the indices depending on the marker's channel, i.e. if DAPI is in
%channel 2, use D=IS{1,1}{2,1} instead.

D=IS{1,1}{1,1}; %DAPI

G=IS{1,1}{2,1}; %CytoK

R=IS{1,1}{4,1}; %Ki67 (or tuft)

FR=IS{1,1}{3,1}; %GSII

% intitialize Flag. 1 if OK, 0 if DAPI area too small (<1/4 of image size)
% or if area of Ki67 signal too large (a hint that it is only background)
Flag=1;

%% DAPI segmentation
h=fspecial('average',4);
T=graythresh(D);
DW=imbinarize(D,0.8*T);
DWf=imfilter(DW,h);
nD=numel(find(DWf));
if nD<numel(D)/4
    Flag=0;
end


%% GSII segmentation
%FRt=imtophat(FR,strel('disk',25));
%FRtm=medfilt2(FRt);
FRtm=medfilt2(FR);
%TF=graythresh(FRtm);

FRW=imbinarize(FRtm);
FRWa=bwareaopen(FRW,10);
FRWaf=imfill(FRWa,'holes');


%% CytoK segmentation
%GW=G>tG;
GW=G>300; %empirical value for thresholding based on current settings and labelling
 
GWc=imclose(GW,strel('disk',5));
GWf=imfilter(GWc,h);

% Rt=imtophat(R,strel('disk',10));
% Rm=medfilt2(Rt);

%% Ki67 segmentation
Rm=medfilt2(R);

RmWf=imbinarize(Rm);
RmWf=bwareaopen(RmWf,10);
RmWff=imfill(RmWf,'holes');
RmWfw=imopen(RmWff,strel('disk',2));
% RmWfw=Nuclear2(RmWf);
% RmWfw=bwareaopen(RmWfw,3);
KiCover=numel(find(RmWf(:)))/numel(find(DWf(:))); %determine fractional coverage of Ki staining

if KiCover>0.2 %probably all background due to absence of specific signal
    Flag=0; 
    tb=minminT(Rm); %use minmin for background subtraction
    Rmb=Rm-1.5*tb;
    RmWf=imbinarize(Rmb);
    RmWf=bwareaopen(RmWf,10);
    RmWff=imfill(RmWf,'holes');
    RmWfw=imopen(RmWff,strel('disk',2));
    %RmWfw=Nuclear2(RmWf);
    %RmWfw=bwareaopen(RmWfw,3);
    KiCover=numel(find(RmWf(:)))/numel(find(DWf(:)));
end

%% Compute metrics

StatsR=regionprops(RmWfw,'Centroid');
%RC=bwconncomp(RmWf);
nbR=numel(StatsR); %total number of Ki67 nuclei





GRW=FRWaf | GWf; %merging Cytokeratin and GSII binary signal to define gland tissue
DGRW=GRW & DWf; % define DAPI positive gland tissue
RmWf2=bwmorph(RmWfw,'thicken',2);
Marker=RmWfw & DGRW;
RmWfR=imreconstruct(Marker,RmWf2); %get Ki67 positive nuclei in the DAPI positive gland tissue
RC=bwconncomp(RmWfR);


GFrac=numel(find(FRWaf))/numel(find(DWf)); %fraction of GSII positive signal over DAPI
GInt=mean(FR(FRWaf)); %mean intensity of GSII signal
KIG=RC.NumObjects/numel(find(DGRW)); %number of Ki67 nuclei in DAPI positive gland tissue

Stats.GSDAPIRatio=GFrac;
Stats.GSmInt=GInt;
Stats.NbKITotal=nbR;
Stats.KIDAPI=nbR/numel(find(DWf)); %number of Ki67 nuclei relative to DAPI
Stats.Flag=Flag; %Flagged if Stats.Flag==0
Stats.KIDAPIGlands=KIG;

%% Plotting
%figure, imshow(R,[]);drawnow
%figure, imshow(RmWf,[]);figure, imshow(FRWaf);figure, imshow(DWf);
B=bwboundaries(FRWaf);

 figure, imshow(DWf);hold on;
for i=1:length(B)
    bb=B{i};
    plot(gca,bb(:,2),bb(:,1),'-','Color',[0.8510    0.3255    0.0980],'Linewidth',3);
end
CC=[StatsR.Centroid];
scatter(gca,CC(1:2:end),CC(2:2:end),'o','MarkerFaceColor',[0.4667    0.6745    0.1882],'MarkerEdgeColor','y');

%t1st=[int2str(nbR) ' (' int2str(RC.NumObjects) ')'];
t1st=int2str(nbR);
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
    

