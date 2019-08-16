function SMPar = StromaModel_ParamEstim()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% StromaModel_ParamEstim: This scrip runs the parameter estimation of the 
% tumour-stroma model presented in the PhD thesis:
% "Mathematical models for heterogeneous preclinical cancers" by Juan Delgado San Martin
% sumbited for the degree of PhD in physics to the university of Aberdeen.
%
% This piece of work will be submitted to npj: systems biology journal under the name:
% "Tumour-stroma dual relationship can be explained with a multiscalar cellular automaton" in 2016
%
% There is unrestricted license to use this script and modify it as long as the Author is acknowledged
% and either of the above publlications correctly cited.
% 
%
%        INPUTS: void - none 
%        OUTPUT:
%             SMPar - Parameters of the model
% 
% May 2015
% Modified November 2015
% AstraZeneca, Alderley Park and Cambridge
% Juan Delgado-SanMartin, PhD Student
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add path
% addpath(genpath('\\emea.astrazeneca.net\uk\Alderley Park\Users 11\knmg297\Documents\PhD\Programs\Matlab\Functions'));

%%  ----------------------- WB module ------------------------------------
        % Load an image
        D = {'3x21_3x1.tif'};
        WBImages = cell(size(D));L = length(D);
        for l = 1:L
        WBImages(l,1) = {imread(D{l})};
        end;clear l
        
        % Definition of parameters
        P.colour = [1 0 0;0 .5 .5];
        P.K_HO2 = 769; %L atm/mol 769*1.39e-6; % mmHg(sol)/XO2(gas) - Henry's law coefficient
        P.L = L;
        HypoxiaSetPoint = 20;
        
        % Call module
         [h_H, ft] = WBModule(WBImages, HypoxiaSetPoint,P,D); % output in mmHg(aq)

%%  ----------------------- IHC module ------------------------------------
% Prepare directory
            close all

ImName = {'aSMA24hS3.tif' 'aSMA24hS3R.tif'
          'H&E_24hS3.tif' 'H&E_24hS3R.tif'                                   
          'HIF124h33.tif' 'HIF124h33R.tif'
          'H&E_48hS3.tif' 'H&E_48hS3R.tif'
          'HIF148h3.tif'  'HIF148h3R.tif'
          'aSMA48hS3.tif' 'aSMA48hS3R.tif'};
    
        % Set up
        % Parameters: These parameters have been calibrated according to 
        % expert pathologists
        Par.im2bw = .05;
        Par.disk = 20;
        Par.bwareaopen = 5000;
        Par.TargetColour0 = ([24,18,21;150,50,250;241,243,241]+1)/256;% IHC
        Par.TargetColour1 = ([255,0,102;0,0,0;102,0,204;]+1)/256;% H&E
        Par.HypThr = 20; % 20% of hypoxia expression
        Par.NecThr = 10; % 10% of necrosis expression
        Par.n = 50; % number of layers  
        
 % Call the module
  [SMPar.beta_H,SMPar.beta_N,SMPar.k_S,HIFHM,IHCCode,hours,Output,Profiles,Boundary] = IHCModule(ImName,Par);

%%  ----------------------- PDE module ------------------------------------
% Find kr'
% load rightnow#
ft.a = 100;
ft.b = -24;

%% just for the time being
for i = 1:5;Boundary(i).cmlayer = 3e-3;end
        for i = 2:size(Profiles,3)
        %     try
        close all
        [k_R1D(i), k_R2D(i)] = PDEModule(HIFHM{i},Profiles(:,:,i),Boundary(i),h_H,ft);
        end

%%  ----------------------- Heterogeneity module ------------------------------------

% call the module
SMPar.mus = HetModule(HIFHM,IHCCode,hours,Output);


%% -------------------- GC Module -----------------------------
% Load 2D data
Data = xlsread('\\emea.astrazeneca.net\UK\Alderley Park\Users 11\knmg297\Documents\PhD\Data\GC\HIF MCF7\MCF7_AR_Incucyte growth curve data.xlsx');

% Reorder data
aa = find(Data(:,1)==1);
t = Data(aa-2,3);
for r = 1:length(aa)
    aux3 = Data(aa(r)+2:aa(r)+7,:);
    Cells(r,1) = mean(aux3(~isnan(aux3)));
    StdCells(r,1) = std(aux3(~isnan(aux3)));
end

% Call the module
[SMPar.beta_T,SMPar.alpha_T] = GCModule(t,Cells,StdCells);

end









