function musOut = HetModule(Images,IHCCode,hours, Output)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HetModule: This module runs a fractal dimension calculation on
% an image as described in the PhD thesis:
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
%         INPUTS:
%             Images - cell array with the segemented images to analyse
%             IHCCode - Code IHC matrix
%         OUTPUT:
%             musOut - parameter motility
% 
% November 2015
% AstraZeneca, Cambridge, UK
% Juan Delgado-SanMartin, PhD Student
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the fractal heterogeneity algorithm
aux = Images(logical(IHCCode(1,:)));
for kk = 1:size(aux,2)
      HF(kk) = HetCalc(aux{kk})./2;
end

% Calculate motility parameter mu_s
    hh = hours(logical(IHCCode(1,:)));
    oo = Output(logical(IHCCode(1,:)))';oo = cell2mat(oo);oo = sum(oo(:,1:3),2);
    hu = unique(hh);
    for i = 1:length(hu)
        HFm(i) = mean(HF(hh==hu(i)));
        stdE(i)=std(HF(hh==hu(i)));
        oom(i) = mean(oo(hh==hu(i)));
    if i ~=1;mus(i-1) = mean((HF(hh==hu(i))-HF(hh==hu(i-1)))'...
            ./(hu(i)-hu(i-1))./oo(hh==hu(i)));end
    end

% Parameter definition
    musOut = mus(1);

% plot it
figure;hold all;set(gca,'FontSize',12)
plot(hh,HF,'ko')
errorbar(hu,HFm,stdE,'ko','MarkerFaceColor','k')
xlim([0 inf])
ylim([.7 .75])
set(gca,'YTick',linspace(.7,.75,3),'XTick',hu)
xlabel('hours')
ylabel('\Psi_{F}')
print('\\emea.astrazeneca.net\uk\Alderley Park\Users 11\knmg297\Documents\PhD\FractalHetaSMA.tif','-dtiff','-r300')


end

function HF = HetCalc(Image)
% Get the fractal dimension
[r,n] = boxcount(Image); 

% Fit
ft = fit(log(r)',log(n'),'poly1');

% Heterogeneity factor
HF = 2+ft.p1;
end
