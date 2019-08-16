function InitSDARS()
%% InitSDARS: initialisation of the SDARS algorithm
% SDARS - Spatial Distribution of Apoptosis Relative to Stroma
%
% The algorithm is a very simple approach to quantify distance of 
% biomarker relative to other tissues' boundary. The system is then
% normalised for the size of each island containing the marker.
% In the context of the publication, we quantified CC3 sitting on the
% tumour epithelium and measure the distance from the stroma.
% This is just an example of what the technique is capable of.
%
% Example
% See one example by loading Example.tif on the GUI that pops out
% when you run this script(press F5). This example differs from the results
% presented in the paper.
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
% Created by: Juan A Delgado
% juan.ads.delgado@astrazeneca.com
% PhD Student at AstraZeneca, United Kingdom
% Created 2013
% Last modified: 5/Dec/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath([cd '\backupfiles'])
run('SDARSv1.m')

