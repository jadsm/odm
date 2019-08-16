function [beta_H,beta_N,k_S,HIFHM,IHCCode,hours,Output,Profiles,Boundary] = IHCModule(ImName,Par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IHC Module: calculates the IHC profiles and parameters from tissue
% engineered slices for the tumour-stroma model presented in the PhD thesis:
% "Mathematical models for heterogeneous preclinical cancers" by Juan Delgado San Martin
% sumbited for the degree of PhD in physics to the university of Aberdeen.
%
% This piece of work will be submitted to npj: systems biology journal under the name:
% "Tumour-stroma dual relationship can be explained with a multiscalar cellular automaton" in 2016
%
% There is unrestricted license to use this script and modify it as long as the Author is acknowledged
% and either of the above publlications correctly cited.
% 
%       INPUTS: 
%             ImName: cell array of names of images. Images should come in pairs, 
%             named equally, one with a ruler on and an R at the end of the name.
%             Par: structure of parameters
%       OUTPUTS: 
%             beta_H: delay in hypoxia development
%             beta_N: delay in necrosis development
%             HIFHM: Heat maps of hypoxia
%             IHCCode: matrix with the IHC codings
%             k_S: constant of stromal recruitment
%             hours: vector with the times
%             Output: cell array of the total quantities of all the quantifies images
%             Profiles: measured spatial profiles 
%
%
%
% May 2015
% Last modified November 2015
% AstraZeneca, Alderley Park and Cambridge, UK
% Juan Delgado-SanMartin, PhD Student
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocation of code
IHCCode = zeros(4,size(ImName,1));
Output = cell(0); 
NameFile = ['AllProfiles' datestr(now,'yyyymmdd') '.mat'];
if exist(NameFile,'file');load(NameFile);end
        % call for every image
for kk = 1+length(Output):size(ImName,1)
         %% Load Image
            Name = ImName{kk};NameR = [Name(1:end-4) 'R' Name(end-3:end)];
            Images.I = imread(Name);Images.IR = imread(NameR);
            hours(kk) = sscanf(Name,'%*4s %u %*s',[1,Inf]);
            IHCType = Name(1:4);
           
%% --------------- Image processing ----------------------------
        % Conversion to code
        IHCCode(:,kk) = str2code(IHCType);
        
        % Call the image processing
        [Images,Par] = ImageProcessing(Images,Par,IHCCode(:,kk));
        HIFHM(kk) = {Images.HIFHM};
   
    if IHCCode(1,kk)
%% ------------------- Total Quantification ------------------------
        Output(kk) = {TotalQuantification(Images)};
        
    else
%% ------------------- Spatial quantification ----------------------
        [Profiles(:,:,kk),thickness,boundary] = SpatialQuantification(Images,Par);
        thicknessCell(kk) = {thickness};
        
        Boundary(kk).micronpixel = Par.micronpixel;
        Boundary(kk).B = boundary.B;
        Boundary(kk).centroid = boundary.centroid;
        Boundary(kk).Ind = boundary.Ind;
        Boundary(kk).cmlayer = Par.micronpixel*mean(thickness)*1e-4; % cm
 
%% ------------------ Curve integration ----------------------------
        Output(kk) = {CurveIntegration(Profiles(:,:,kk),thickness,Par)};
    end
        %% Save profiles
save(NameFile)
close all
    
end
 
    %% Merge outputs
%     [beta_H,beta_N,k_S] = CalculateParameters(Profiles,Output,SegmentedImages,IHCCode,hours,thicknessCell,boundary);
    [beta_H,beta_N,k_S] = CalculateParameters(Output,hours,thicknessCell,IHCCode,Par);
    
end
function [Images,Parameters] = ImageProcessing(Images,Parameters,IHCCode)
%% Aspect ratio : this section will extract the micron/pixel ratio by comparing two images
    Mask = mean(Images.IR,3)-mean(Images.I,3);
    imshow(Mask)
    h = helpdlg('Now select the extense of the left-right points of the measurement bar. Finish with "double click" or "enter"','Aspect Ratio');uiwait(h);clear h
    [pixX,~,~] = impixel();pixX = sort(pixX);
    input = inputdlg('What is on the label?','Microns');
    Parameters.micronpixel = str2num(input{:})./(pixX(2)-pixX(1));
    clear pixX input Mask IR
    close

%% Image modification
        % Colour deconvolution
        [~] = createfigure(Images.I,[],1);
        
        % first call
        II = CDColour_Callback(0,0,Images.I,transpose(eval(['Parameters.TargetColour' num2str(IHCCode(4))])));
        % wait
        uiwait(gcf) 
               
        % Binarisation      
        h = createfigure(II,Parameters,2);
        % First call
        BW = binarisation_Callback(0,0,sum(II(:,:,2:3),3),h);
        % wait
        uiwait(gcf)
        % Last call
        BW = binarisation_Callback(0,0,sum(II(:,:,2:3),3),h);
                
        % HIF is a 3D logical matrix with the four stainings PosHIF30, PosHIF70,
        % PosHIF100, NegHIF. with different colours.
        [nx,ny,nz] = size(II);
        Colours = [1 0 0;1 .4 0; 1 .8667 .1098; 0 0 1];
        HIFHM = II;Aux2 = HIFHM(:,:,1);HIFHM(:,:,1) = 0;HIFHM = rgb2gray(HIFHM);HIFHM = HIFHM./max(HIFHM(:));
        % Create the mask for the three stainigns: Mask
        Mask(:,:,1) = double(im2bw(HIFHM,.7)); % Red, strong opositive
        Mask(:,:,2) = double(im2bw(HIFHM,.3))-Mask(:,:,1); % Orange: medium Positive [1.0000    0.4000         0]
        Mask(:,:,3) = double(HIFHM>0.05)-Mask(:,:,2)-Mask(:,:,1); % Yellow: weak positive [1.0000    0.8667    0.1098]
        Mask(:,:,4) = Aux2; % Blue: negative
        Mask(Mask<=0) = 0;
        % Transform it into colours
        HIF = zeros(size(II)); % Initialise
        m = zeros(size(II)); % Initialise
        for i = 1:4
            c = Colours(i,:);
            m(:,:,1) = c(1);m(:,:,2) = c(2);m(:,:,3) = c(3);
            HIF = HIF+repmat(Mask(:,:,i),1,1,3).*m;
        end; clear i

        clear Aux2

        % Main modifications
        [~,~,~,~,~,~,f,A] = Main_Modification (BW,1);

        %% rotate the images
        Images.IRot = imrotate(BW,atan(-1/f(1))*180/pi+A,'crop');
        Images.HIFRot = imrotate(HIF,atan(-1/f(1))*180/pi+A,'crop');
        Images.PRot = imrotate(Mask,atan(-1/f(1))*180/pi+A,'crop');
        Images.HIFHM = imrotate(HIFHM,atan(-1/f(1))*180/pi+A,'crop');
       
        %% Crop the image     
        [ny,nx,nz] = size(Images.IRot);% Size of the image
        figure;imshow(Images.IRot);
        txt = text(round(nx/2),100,'Air');set(txt,'FontSize',20,'FontWeight','Bold','Color',[0 .2 1]);
        txt = text(round(nx/2),ny-100,'Filter');set(txt,'FontSize',20,'FontWeight','Bold','Color',[1 .2 0]);
        h = helpdlg('Now, I will ask you to select 2 points to crop the image. First, left, then right','Croping the image');uiwait(h);clear h
        [xi,yi,~] = impixel;
        Images.IRot = Images.IRot(:,xi(1):xi(2));
        Images.HIFRot = Images.HIFRot(:,xi(1):xi(2),:);
        Images.PRot = Images.PRot(:,xi(1):xi(2),:); 
        Images.HIFHM = Images.HIFHM(:,xi(1):xi(2),:); 
        % Main Modifications Again
        [bx,by,bbx,bby,a,b,f,~] = Main_Modification (Images.IRot,0);

        clear PP PPx PPy PPx2 PPy2 PPRot Aux
        
        % Renaming of parameters
        Parameters.bx = bx;
        Parameters.by = by;
        Parameters.bbx = bbx;
        Parameters.bby = bby;
        Parameters.a = a;
        Parameters.b = b;
        Parameters.f = f;
        Parameters.A = A;           
end
function [Profiles,thickness,boundary] = SpatialQuantification(Images,Parameters)

 % Segregate in top and bottom borders
        % Boundaries NON correspondant to border of the frame
        [ny,nx,nz] = size(Images.IRot);
        bbbx = Parameters.bx;bbby = Parameters.by;
        bbbx(Parameters.bx==nx|Parameters.bx==1) = [];
        bbby(Parameters.bx==nx|Parameters.bx==1) = [];

        % Left border
        left = [Parameters.bx(Parameters.bx==1) Parameters.by(Parameters.bx==1)];
        % Right border
        right = [Parameters.bx(Parameters.bx==nx) Parameters.by(Parameters.bx==nx)];
        % Top border
        top = [bbbx(bbby<=Parameters.b) bbby(bbby<=Parameters.b)];
        % Bottom border
        bttm = [bbbx(bbby>=Parameters.b) bbby(bbby>=Parameters.b)];
            figure;imshow(Images.IRot);hold on
            colornow = [1 0 0];
            plot(top(:,1),top(:,2),'Color',colornow,'LineWidth',2)
            aux = top(round(length(top)/2),:);
            txt = text(aux(1),aux(2),'Air Interface');set(txt,'Color',colornow,'BackGround',[1 1 1])
            colornow = [1 1 0];aux = bttm(round(length(bttm)/2),:);
            plot(bttm(:,1),bttm(:,2),'Color',colornow,'LineWidth',2)
            txt = text(aux(1),aux(2),'Filter Interface');set(txt,'Color',colornow,'BackGround',[1 1 1])
            clear aux txt colornow
       
        %Interpolate
        [C,IA,~] = unique(bttm(:,1));
        aux = bttm(:,2);%bttm(:,1) = C;bttm(:,2) = aux(IA);
        vbttm = interp1(C,aux(IA),top(:,1));
        % Redefine bttm
        clear bttm
        bttm = zeros(size(top));
        bttm (:,1) = top(:,1);bttm (:,2) = vbttm;

        % Discretisation
        l = (bttm(:,2)-top(:,2))./Parameters.n; % length of nodes
        GY = repmat(top(:,2),1,Parameters.n)+l*(1:Parameters.n);
        GX = repmat(bttm(:,1),1,Parameters.n);

        % bands
        [bwx,bwy,~] = size(Images.HIFRot);
%         avTh = mean(GY(:,2)-GY(:,1)); % Average Thickness
%         stETh= std(GY(:,2)-GY(:,1))./sqrt(length(GY(:,2))); % Standard Error
%         thickness = median(GY(:,2:end)-GY(:,1:end-1));
        thickness = (1:50)'.*mean(sqrt((bttm(:,1)-top(:,1)).^2+(bttm(:,2)-top(:,2)).^2))./50;
        
        
        for i = 2:Parameters.n
        band (:,:,i) = poly2mask([GX(:,i); flipud(GX(:,i-1));GX(1,i)],round([GY(:,i); flipud(GY(:,i-1));GY(1,i)]),bwx,bwy);
        HIFband(:,:,:,i) = repmat(double(band(:,:,i)),1,1,4).*Images.PRot; 
        % Quantify
        aux1 = HIFband(:,:,1,i);aux2 = HIFband(:,:,2,i);
        aux3 = HIFband(:,:,3,i);aux4 = HIFband(:,:,4,i);
        aux5 = band(:,:,i);
        PosHIF30(i,1) = sum(aux1(:)); 
        PosHIF70(i,1) = sum(aux2(:));
        PosHIF100(i,1) = sum(aux3(:));
        NegHIF(i,1) = sum(aux4(:));
        Total(i,1) = sum(aux5(:));
        end; clear aux1 aux2 aux3 aux4 aux5 i
        
        [Profiles] = [PosHIF30, PosHIF70, PosHIF100, NegHIF,Total];
    
        % Computation of the boundary
        boundary.B = [top;right;flipud(bttm);left]; % cellfun(@(x,y)x,B,'UniformOutput',false);
        boundary.centroid = [Parameters.a Parameters.b]; 
        boundary.Ind = [4*ones(size(top));1*ones(size(right));2*ones(size(flipud(bttm)));3*ones(size(left))];

end
function TotalaSMA = TotalQuantification(Images)
% Total quantification
Tum = sum(Images.IRot(:));TotalaSMA = zeros(1,size(Images.PRot,3));
for l = 1:size(Images.PRot,3);aux = Images.PRot(:,:,l);TotalaSMA(l) = sum(aux(:))./(Tum-sum(aux(:)));end
end
function Integral = CurveIntegration(Profiles,thickness,Parameters)
%% CurveIntegration: sum of all the elements multiplied by the
Integral = thickness'*Profiles.*Parameters.micronpixel;
end
    function [bx,by,bbx,bby,a,b,f,Aout] = Main_Modification (I,Rotation)
% Size image
[ny, nx] = size(I);
% Get boundaries
[boundary] = bwboundaries(I);
bound = boundary{1};
bx = bound(:,2); bbx=sort(unique(bx));
by = bound(:,1); bby=sort(unique(by));
% plot(bx,by,'b')
% Get Centroid
figure;imshow(I)
RP = regionprops(I);
a = RP(1).Centroid(1); % Centroid of the circle 
b = RP(1).Centroid(2);
hold on
plot(a,b,'ro')
plot(a,b,'rx')

% Find the minimum distance
distances = sqrt((bx-a).^2+(by-b).^2);
u = find(distances==min(distances));
M = bound(u,:);
plot(M(2),M(1),'go')
plot(M(2),M(1),'gx')
kk = 0;
%% Questions of rotation
if Rotation == 1
while kk == 0
% Construct a line
po = [1 1];
[f] = fsolve(@(p)FCN(p,[M(2) a],[M(1) b]),po);
y = f(1).*bbx-f(2);
ff2 =  -1./f(1).*a-b;
y2 = -1./f(1).*bbx-ff2;

    % Change the colour of preexisting lines
    if exist('pl')
    oldlinecolor = [.3 .3 .3];% somewhat grey
    set(pl,'Color',oldlinecolor,'LineWidth',.5)
    end
% hold on
pl(1) = plot(bbx,y,'r','LineWidth',1);
pl(2) = plot(bbx,y2,'r','LineWidth',1);
answer = questdlg ('The red lines should go through the midle of the slice. Do you like this result?','Rotation','Yes','No','No');
switch answer
    case 'No'
M2 = impoint;
M = fliplr(M2.getPosition);
    case 'Yes'
  kk = 1;
end

end
% Localise the air
Air = questdlg('Is the selected point on the air?');
switch Air
    case 'Yes'
        A = 1;
    case 'No'
        A = 0;
end; clear Air

Aout = double((M(1)>b & A | M(1)<b & ~A))*180;
else
    f=[];
    Aout = [];
    

end

%% Vertical
plot([a, a],[bby(1) bby(end)],'r--')
close
% Detect points to the periphery
% PeriPoints = bound(bound(:,2)==nx|bound(:,2)==1|bound(:,1)==ny|bound(:,1)==1,:);
% GPP = [PeriPoints(2:end,:); PeriPoints(1,:)]-[PeriPoints(end,:); PeriPoints(1:end-1,:)]; % Gradients of peripheral points
% MGGP = sqrt(GPP(:,1).^2+GPP(:,2).^2); % module
% PP = PeriPoints(find(MGGP>(max(MGGP)./2)),:);

    end
    function save2excell (figs,n,Results,Name)
try
Results = array2table([(1:n)',Results]);
% for i = 1:length(YER(:,1));NameHIF(i) = {['HIFSlide' num2str(i)]};NameER(i) = {['ERSlide' num2str(i)]};end
% Names = [{'Stripe'} NameER NameHIF];
Results.Properties.VariableNames = {'Layer','StrongPositive','MediumPositive','WeakPositive','Negative','Total'};
Results.Properties.VariableUnits = { '#' '#' '#' '#' '#' '#'};
[x,y] = size(Results);
Range = ['A1:' char('A' + y-1) num2str(x+1)];
% Write it in Excell
writetable(Results,Name,'FileType','spreadsheet','Sheet',1,'Range',Range)

% Operate with excell
Excel = actxserver('Excel.Application');
Excel.Visible = true;
wb = Excel.Workbooks.Open([cd '\' Name '.xls']);
Excel.worksheets.Item(1).Name = 'Data';
Excel.worksheets.Item(2).Name = 'Images';
Sheets = Excel.ActiveWorkBook.Sheets;
SheetFigures = get(Sheets, 'Item', 2);

% Activate the figures speichern
SheetFigures.Activate;
Shapes = SheetFigures.Shapes;

for ii = 1:length(figs)
print(figs(ii),'-dtiff', 'Auxiliar');
% Add image
% Image Dymensions
Margin = 20;ImageDim = 300;%mm
if ii > 1; Init = Margin*2+ImageDim;else Init = 0; end
if ii <=4; oo = 1;elseif ii<=7; oo = 2; elseif ii<=10; oo = 3;else oo = 4;end
% Shapes.AddPicture([pwd '\Auxiliar.tif'] ,0,1,Margin*(ii-4*(oo-1))+ImageDim*(ii-4*(oo-1)-1),Margin*oo+ImageDim*(oo-1),ImageDim,ImageDim);
Shapes.AddPicture([pwd '\Auxiliar.tif'] ,0,1,Init+(Margin+ImageDim)*(ii-3*(oo-1)-2),Margin*oo+ImageDim*(oo-1),ImageDim,ImageDim);
end

wb.Save
Excel.Quit

movefile([Name '.xls'],'\\emea.astrazeneca.net\uk\Alderley Park\Users 11\knmg297\Documents\PhD\Programs\Matlab\HIF_ER\Results\Faslodex')
catch
    disp('File was not saved to excell, probably there is a failure in the "cloud"')
end

    end
    function [II] = CDColour_Callback(hObject,callbackdata,I,M) 
        function II = MyCD(I,M)   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Deinterlace: remove line artefacts from video scan
                H_deinterlace = [0 1 0; 0 2 0; 0 1 0] ./4;
                [sx,sy,sz] = size(I);
                sample_deinterlace = zeros(size(I));
                for k=1:sz
                    sample_deinterlace(:,:,k) = filter2(H_deinterlace,double(I(:,:,k)),'same');
                end
                % figure
                % imshow(sample_deinterlace)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Convert RGB intensity to optical density (absorbance)
                RGB_OD = -log((sample_deinterlace+1)./256);

                %% Construct color deconvolution matrix
                % Take the average around region of interest
                H2 = ones(10,10) ./ 45;
                RGB_OD_Blur = zeros(size(I));

                for k=1:sz
                    RGB_OD_Blur(:,:,k) = filter2(H2,RGB_OD(:,:,k),'same');
                end
                % % imshow(sampleRGB_OD_Blur)
                % 
                RGB_OD = RGB_OD_Blur;

                %% Create Deconvolution matrix
                % M = [Val.TargetColour(1,:)/norm(Val.TargetColour(1,:)); Val.TargetColour(2,:)/norm(Val.TargetColour(2,:)); Val.TargetColour(3,:)/norm(Val.TargetColour(3,:))];
                D = inv(M);

                PP2N_OD = zeros(size(I));


                for i=1:sx
                    for j=1:sy
                        RGB1 = reshape(RGB_OD(i,j,:),sz,1);
                        PP2N = D * RGB1;

                        PP2N_OD(i,j,1) = PP2N(1);
                        PP2N_OD(i,j,2) = PP2N(2);
                        PP2N_OD(i,j,3) = PP2N(3);
                    end
                end

                PP2N_OD = (PP2N_OD - 0.05) .* (PP2N_OD > 0.05);
                maxPP2N_OD = max(max(PP2N_OD));
                % Enchance contrast
                for k = 1:sz
                PP2N_OD(:,:,k) = PP2N_OD(:,:,k) ./ maxPP2N_OD(k);
                end

                %% Extract tumor cells that are stained
                II = PP2N_OD(:,:,1:3);
        end
    
    child = get(gcf,'Children');
    
    if M == 0; %calibration
       helpdlg('Select some points of the POSITIVE stain')
       uiwait(gcf);
       figure;imshow(I)
       px1 = impixel;
       helpdlg('Select some points of the CONTRAST stain')
       uiwait(gcf);
       imshow(I)
       px2 = impixel;
       helpdlg('Select some points of the NEGATIVE stain')
       uiwait(gcf);
       imshow(I)
       px3 = impixel;      
       close  
       % Refresh M
       M = [mean(px1,1);mean(px2,1);mean(px3,1)];      
    end
    
    
    II = MyCD(I,M);  
 
    set(gcf,'currentaxes',child(2))
    imshow(II);title('Colour Deconvoluted')
    txt = text(-400,1,'Positive');set(txt,'Color',[0 0 1],'FontWeight','Bold')
    
    set(gcf,'currentaxes',child(1))
    imshow(I);title('Source Image')
    end
    function [F] = FCN(p,x,y)
F(1,1) = y(1)-p(1)*x(1)+p(2);
F(2,1) = y(2)-p(1)*x(2)+p(2);
    end
    function BW = binarisation_Callback(hObject,callbackdata,II,h)     
    % Transformation
    pim2bw = get(h(1),'Value');
    pbwareaopen = round(get(h(2),'Value'));
    pdisk = round(get(h(3),'Value'));

    
    % Calculation of variables
        BW = im2bw(II,pim2bw);
        se = strel('disk',pdisk);
        BW = imclose(BW,se);
        BW = bwareaopen(BW, pbwareaopen);
        clear se
        imshow(BW)
    % update displays
    set(h(4),'String',num2str(pim2bw))
    set(h(5),'String',num2str(pbwareaopen))
    set(h(6),'String',num2str(pdisk))
    
    end 
    function [beta_H,beta_N,k_S] = CalculateParameters(Output,hours,thicknessCell,IHCCode,Parameters)
       
        % HIF
        o24 = cell2mat(Output(IHCCode(2,:)&hours==24)');
        o24 = sum(sum(o24(:,1:2),2)./o24(:,end))*100;
        o48 = cell2mat(Output(IHCCode(2,:)&hours==48)');
        o48 = sum(sum(o48(:,1:2),2)./o48(:,end))*100;       
       
        beta_H = 250*24./(o48-o24)./Parameters.HypThr;

        % H&E
        o24 = cell2mat(Output(IHCCode(4,:)&hours==24)');
        o24 = sum(sum(o24(:,1:2),2)./o24(:,end))*100;
        o48 = cell2mat(Output(IHCCode(4,:)&hours==48)');
        o48 = sum(sum(o48(:,1:2),2)./o48(:,end))*100;   
        beta_N = 350*24./(o48-o24)./Parameters.NecThr;
        
        % Total quantification of the stroma
        aux24 = mean(cell2mat(Output(logical(IHCCode(1,:))&hours==24)'));
        aux48 = mean(cell2mat(Output(logical(IHCCode(1,:))&hours==48)'));
        k_S = (aux24-aux48)/24;%cm3stroma/cm3tumour/hour    
    end
    function IHCCode = str2code(IHCType)
    IHCCode = zeros(4,1);
        switch IHCType
            case 'aSMA';IHCCode(1,1) = 1;
            case 'HIF1';IHCCode(2,1) = 1;
            case 'ER__';IHCCode(3,1) = 1;
            case 'H&E_';IHCCode(4,1) = 1;
        end
    end
    function h = createfigure(II,Par,fl)
    %% createfigure: global support function
        figure
        uicontrol('Style','pushbutton','String','Go on!','CallBack','uiresume(gcbf)');
        
        if fl == 2 % Binarisation
            h(1) = uicontrol('Style','Slider','Position', [100 20 300 20],'Value',Par.im2bw);
            h(2) = uicontrol('Style','Slider','Position', [420 20 300 20],'Value',Par.bwareaopen,'Max',10000);
            h(3) = uicontrol('Style','Slider','Position', [740 20 300 20],'Value',Par.disk,'Max',30);
 
            h(4) = uicontrol('Style','text','Position', [250 40 30 20],'String',num2str(get(h(1),'Value')));
            uicontrol('Style','text','Position', [230 60 80 20],'String','Image to BW')
            h(5) = uicontrol('Style','text','Position', [570 40 30 20],'String',num2str(get(h(2),'Value')));
            uicontrol('Style','text','Position', [550 60 80 20],'String','Area open')
            h(6) = uicontrol('Style','text','Position', [890 40 30 20],'String',num2str(get(h(3),'Value')));
            uicontrol('Style','text','Position', [870 60 80 20],'String','Disk close size')
            set(h(1),'CallBack',{@binarisation_Callback,II,h});
            set(h(2),'CallBack',{@binarisation_Callback,II,h});
            set(h(3),'CallBack',{@binarisation_Callback,II,h});
        elseif fl == 1 % colour deconvolution
        uicontrol('Style','pushbutton','String','Calibrate!','Position', [20 50 60 20],'CallBack',{@CDColour_Callback,II,0});            
        h = [];
        subplot(121);subplot(122);
        end 
    end
