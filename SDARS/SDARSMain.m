function SDARSMain(Macro)
%% SDARSMain
% this is the main script for SDARS
%
% Copyright
% This code is protected by AstraZeneca's copyright
% The code, however, can be freely distributed, used or modified at will,
% as long as the original publication is correctly cited. 
% The citation should say (or similar):
% "Delgado San Martin et al. (2015)
% Tumour stromal morphology impacts nanomedicine cytotoxicity
% in patient-derived xenografts. Nanomedicine: NBM."
%
% Juan A Delgado
% AstraZeneca, UK
% Created 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of Image
I = Macro{1};
if Macro{34} == 1; Figs = zeros(10,1); else; Figs = [];end

%% Aspect Ratio

%% Pre-Processing
[I,Bw,Figs] = PP(I,Macro,Figs);

%% Colour Deconvolution
[II,S,Figs] = CD(I,Figs,Macro);

%% Locate tissue
[III,Figs] = LT(II,Bw,Macro,Figs);

%% Correct by hand

%% Normalisation
[FN,DistN,bb,no] = Norm(III,Macro);

%% Compute Distances
Results = CDist(S,FN,DistN,bb,no,Macro);

%% Export to excell
fname = Macro{2};
filename = [fname(1:end-4) datestr(now,30)];
xlswrite([cd '\' filename],Results)
%% functions
function [I,Bw,Figs] = PP(I,Macro,Figs)
 % Gaussian filter
 H = fspecial('gaussian',[3 3],Macro{3});
 I = imfilter(I,H,'replicate');
 
 % Binarise
 Bw = im2bw(I,Macro{6});
 Bw = 1-Bw;
 figure;imshow(Bw)
 
 % Close morphology
 se = strel('disk',Macro{4});
 Bw = imclose(Bw,se);
 figure;imshow(Bw);

 % Open Area
 Bw = bwareaopen(Bw,Macro{5});
 figure;imshow(Bw)
function [II,S,Figs] = CD(I,Figs,Macro)
%% Rename
TargetColour = reshape(cell2mat(Macro(8:16)),3,3)';
Names = Macro(24:26)';

%% Deinterlace: remove line artefacts from video scan
H_deinterlace = [0 1 0; 0 2 0; 0 1 0] ./4;
[x,y,z] = size(I);
sample_deinterlace = zeros(x,y,z);
for k=1:z
    sample_deinterlace(:,:,k) = filter2(H_deinterlace,double(I(:,:,k)),'same');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert RGB intensity to optical density (absorbance)
RGB_OD = -log((sample_deinterlace+1)./256);

%% Construct color deconvolution matrix
% Take the average around region of interest
H2 = ones(10,10) ./ 45;
RGB_OD_Blur = zeros(x,y,z);

for k=1:z
    RGB_OD_Blur(:,:,k) = filter2(H2,RGB_OD(:,:,k),'same');
end
RGB_OD = RGB_OD_Blur;

%% Create Deconvolution matrix
M = [TargetColour(1,:)/norm(TargetColour(1,:)); TargetColour(2,:)/norm(TargetColour(2,:)); TargetColour(3,:)/norm(TargetColour(3,:))];
D = inv(M);

PP2N_OD = zeros(x,y,z);
for i=1:x
    for j=1:y
        RGB1 = reshape(RGB_OD(i,j,:),z,1);
        PP2N = D * RGB1;

        PP2N_OD(i,j,1) = PP2N(1);
       	PP2N_OD(i,j,2) = PP2N(2);
       	PP2N_OD(i,j,3) = PP2N(3);
    end
end

PP2N_OD = (PP2N_OD - 0.05) .* (PP2N_OD > 0.05);
maxPP2N_OD = max(max(PP2N_OD));
% Enchance contrast
for k = 1:z
PP2N_OD(:,:,k) = PP2N_OD(:,:,k) ./ maxPP2N_OD(k);
end
%% Extract tumor cells that are stained
II = PP2N_OD(:,:,1:3);

% plot it
figure;for i2 =1:3;subplot(3,1,i2);imshow(II(:,:,i2));title(Names{i2},'FontSize',12,'FontWeight','Bold');end


%% Apply thresholds
% identify DAB
n = find(cell2mat(Macro(30:32)));
S = II(:,:,n)>Macro{18};
W = II(:,:,n)>Macro{17}&II(:,:,n)<Macro{18};
M = zeros(size(II));
M(:,:,1) = S;
M(:,:,2) = W;
M(:,:,1) = M(:,:,1)+W;
figure;imshow(M)
function Results = CDist(S,FN,DistN,bb,no,Macro)

for ii = 1:length(no)
noNow = no{ii};
ind = sub2ind(size(S),noNow(:,2),noNow(:,1));
    
% Localise "hotspots"
ix = find(S(ind));
h = [noNow(ix,1),noNow(ix,2)];

% calculate the distances
Dist = zeros(length(bb),length(h(:,1)));
for j = 1:length(h(:,1))%hotspots
Dist(:,j) = sqrt((h(j,1)-bb(:,1)).^2+(h(j,2)-bb(:,2)).^2);
end
MinDist = min(Dist)';

% compute histograms
[D(ii,:),~] = hist(MinDist,DistN);

clear h ix iy ind noNow j 
end

MyF = D./FN;MyF(isnan(MyF)) = 0;

% Compute totals in the image
TotF = sum(MyF);

% Plot it 
figure
ax = subplot(311);set(ax,'FontSize',14,'FontWeight','Bold')
bar(DistN,sum(D),'FaceColor',[0 0 1])
ylabel('# Raw Pix.')
xlim([0 max(DistN)])
ax = subplot(312);set(ax,'FontSize',14,'FontWeight','Bold')
bar(DistN,sum(FN),'FaceColor',[1 0 0])
ylabel('# MC points')
xlim([0 max(DistN)])
ax = subplot(313);set(ax,'FontSize',14,'FontWeight','Bold')
bar(DistN,sum(D)./sum(FN),'FaceColor',[.5 0 1])
ylabel('# Norm. Pix')
xlabel(['Distance to ' Macro{23+find(cell2mat(Macro(27:29)))} ' in ' Macro{20}])
xlim([0 max(DistN)])


% Results
Results = cell(104,4);
Results(1,1) = {'For more details: Delgado San Martin et al (2015) Tumour stromal morphology impact nanomedicine cytotoxicity in patient-derived xenografts. Nanomedicine: NBM.'};
Results(2,1) = {'Created by Juan Delgado. Email: juan.ads.delgado@astrazeneca.com'};
Results(3,:) = {'Distance','Frequency Norm','Frequency Pix','Norm Pix'};
Results(4,:) = {Macro{20},'#','#',' '};
Results(5:end,:) = num2cell([DistN; sum(D); sum(FN); sum(D)./sum(FN)])';
function [III,Figs] = LT(II,Bw,Macro,Figs)
III(:,:,1) = im2bw(II(:,:,1),Macro{21}).*Bw;III(:,:,1) = bwareaopen(III(:,:,1),Macro{5});
III(:,:,2) = im2bw(II(:,:,2),Macro{22}).*Bw;III(:,:,2) = bwareaopen(III(:,:,2),Macro{5});
III(:,:,3) = im2bw(II(:,:,3),Macro{23}).*Bw;III(:,:,3) = bwareaopen(III(:,:,3),Macro{5});
figure;
subplot(321);imshow(II(:,:,1));title(Macro{24})
subplot(323);imshow(II(:,:,2));title(Macro{25})
subplot(325);imshow(II(:,:,3));title(Macro{26})
subplot(322);imshow(III(:,:,1));title(Macro{24})
subplot(324);imshow(III(:,:,2));title(Macro{25})
subplot(326);imshow(III(:,:,3));title(Macro{26})
function [FN,DistN,bb,no] = Norm(III,Macro)
% Select tissue of interest
r = find(cell2mat(Macro(30:32)));
N = logical(III(:,:,r));
no = struct2cell(regionprops(N,'PixelList'));

% Locate "Tissue Boundary"
r = find(cell2mat(Macro(27:29)));
A = logical(III(:,:,r));
[B] = bwboundaries(A);
bb = fliplr(cell2mat(B));lb = length(bb);

% Montecarlo points
for ii = 1:length(no)
noNow = no{ii};

%Select montecarlo points
if length(noNow) > Macro{7};K = Macro{7};else;K = length(noNow);end
MC = noNow(randperm(size(noNow,1),K)',:);

% Compute distances
DN = zeros(lb,length(MC));
for j = 1:length(MC)
DN(:,j) = sqrt((MC(j,1)-bb(:,1)).^2+(MC(j,2)-bb(:,2)).^2);
end

DMin = Macro{19}.*min(DN)';
RawN(ii) = {DMin};
clear ax ay MC DN n F DMin j noNow DN2 lmc lbb
end

% Compute histograms
DistN = linspace(0,250,100); % United vector
for r = 1:length(RawN)
[F,~] = hist(RawN{r},DistN);
FN(r,:) = F;clear F
end
