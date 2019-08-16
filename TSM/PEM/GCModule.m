function [beta_T,alpha_T] = GCModule(t,Cells,StdCells)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GCModule: calculates the growth curve parameters of the 
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
%   INPUTS:
%         t - time
%         Cells - confluence in %
%         StdCells - standard deviation on the replicates
%     OUTPUTS:
%           beta_T - delay on the tumour progression
%           alpha_T - parameter of growth
% 
% May 2015
% AstraZeneca, Alderley Park
% Juan Delgado-SanMartin, PhD Student
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uu = 0;
for e = 2:find(t==50) 
    uu = uu+ (StdCells(e)+StdCells(e-1))./2*(t(e)-t(e-1));
end
uu = uu/50;

% Calculation of alpha
alphaT = 1/(t(1:find(t==50))\log(Cells(1:find(t==50))));

% Calculation of the delay
betaT = SMPar.alphaT*log(uu/1);

% fit 
ft = fit(t(1:find(t==50)),Cells(1:find(t==50)),'exp1');% plot

% plot
figure;hold('all');set(gca,'FontSize',12)
errorbar(t,Cells,StdCells,'ko')
plot(t,ft.a*exp(ft.b*t),'r')
ylim([0,105])
xlim([0 200])
xlabel('hours')
ylabel('% MCF7 Cells')
legend('Data','fit')
print('\\emea.astrazeneca.net\uk\Alderley Park\Users 11\knmg297\Documents\PhD\Cellsalphabeta.tif','-dtiff','-r300')

% Deviation plot
figure;hold('all');set(gca,'FontSize',12)
plot(t,StdCells,'ko-')
ylabel('Std of % MCF7 Cells')
xlabel('hours')
print('\\emea.astrazeneca.net\uk\Alderley Park\Users 11\knmg297\Documents\PhD\Stdalphabeta.tif','-dtiff','-r300')

end