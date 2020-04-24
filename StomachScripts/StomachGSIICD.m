function Stats=StomachGSIICD(czidir,k)

%This function is used to to quantify the overlap between GSII and CD44v
%from gastric IF sections labelled with DAPI, GSII, CD44v and a 4th marker

%INPUT:
% czidir  -a structure array returned by the "dir" function that contains the
% list of all the .czi files to analyze
% k   -the index of the file/image to analyze

%OUTPUT:
% Stats  -a structure array containing the quantified metrics in each
% field.
% Example:
% Stats=StomachGSIICD(dirczi,1);



%% read the .lsm file - requires bioformats readers
IS=bfopen(czidir(k).name);

% IS=bfopen('1527_corpus_10x_04.czi');

%Switch the indices depending on the marker's channel, i.e. if DAPI is in
%channel 2, use D=IS{1,1}{2,1} instead.

D=IS{1,1}{1,1}; %DAPI

GSII=IS{1,1}{2,1}; %GSII

CD=IS{1,1}{3,1}; %CD44v

Tff=IS{1,1}{4,1}; %Tff3


%% binarize DAPI signal
DW=imbinarize(D);
DWm=medfilt2(DW);

%% Segment and binarize CD44V signal
CDm=medfilt2(CD);
CDim=imadjust(CDm);
%se1=strel('line',5,0);
%se2=strel('line',5,90);
CDmt=imtophat(CDm,strel('disk',2));

CDmta=imadjust(CDmt);

TCD=minminT(CDmta);
CDb=CDmta-0.8*TCD; % background subtraction
CDbmedW=imbinarize(CDb);
CDbmedW=bwareaopen(CDbmedW,10);

%% Segment and binarize GSII signal
GSIIa=imadjust(GSII);
TGS=minminT(GSIIa);
GSIIb=GSIIa-TGS;
GW=imbinarize(GSIIb);
GW=bwareaopen(GW,50);
GWd=imdilate(GW,strel('disk',10)); %Dilate GSII signal to cover plasma membrane

%% Get overlap and metrics

GWCDW=GWd & CDbmedW;

Stats.numDAPI=numel(find(DWm)); %number of DAPI positive pixels
Stats.numCD=numel(find(CDbmedW)); %number of CD44v positive pixels
Stats.numGS=numel(find(GW)); %number of GSII positive pixels
Stats.numoverlap=numel(find(GWCDW)); %number of double positive GSII/CD44v pixels
Stats.PercCD=Stats.numCD/Stats.numDAPI; %fraction of CD44v pixels over DAPI
Stats.PercGS=Stats.numGS/Stats.numDAPI; %Fraction of GSII positive pixels over DAPI
Stats.McoefCD=Stats.numoverlap/Stats.numCD; %Manders coef for CD44v
Stats.McoefGS=Stats.numoverlap/Stats.numGS; %Manders coef for GSII

%% Plotting
%Stack=cat(3,GSII,CD);
A=imoverlay(GSIIa,GW,'r');
B=imoverlay(CDim,CDbmedW,'r');
figure, montage({A,B});
%figure, montage(Stack);
%figure, imshowpair(GWd,CDbmedW);
%figure,imshow(CD,[])
%figure, imshow(RmWf,[]);figure, imshow(FRWaf);figure, imshow(DWf);
% B=bwboundaries(FRWaf);
% 
%  figure, imshow(DWf);hold on;
% for i=1:length(B)
%     bb=B{i};
%     plot(gca,bb(:,2),bb(:,1),'-','Color',[0.8510    0.3255    0.0980],'Linewidth',3);
% end
% CC=[StatsR.Centroid];
% scatter(gca,CC(1:2:end),CC(2:2:end),'o','MarkerFaceColor',[0.4667    0.6745    0.1882],'MarkerEdgeColor','y');
% 
% %t1st=[int2str(nbR) ' (' int2str(RC.NumObjects) ')'];
 t1st=num2str(Stats.McoefCD,3);
 t2st=num2str(Stats.McoefGS,3);
t1=text(100,50,t1st,'Color',[0.4667    0.6745    0.1882],'FontSize',24,'FontWeight','Bold');
t2=text(100,125,t2st,'Color',[0.8510    0.3255    0.0980],'FontSize',24,'FontWeight','Bold');
% t1.BackgroundColor='w';
% t2.BackgroundColor='w';
Stats=struct2table(Stats);
    

