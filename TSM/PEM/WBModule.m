function [h_H,ft2] = WBModule(WBImages, HypoxiaSetPoint,P,D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WBModule: western Blot module quantifies the western blots and calculates of the 
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
%     INPUTS:
%            WBImages: is a cell array of dimensions Nx1, being N the number of images
%            HypoxiaSetPoint: this is a threshold of hypoxia, default is 20% 
%            P: structure of parameters with fields:
%                 colour - prefered colour for the plor
%                 K_HO2 - Henry's law coefficient
%                 L - length of D
%            D: Images to analyse in a cell array
%     OUTPUT:
%            h_H: level of oxygen at which hypoxia starts to be expressed in mmHg
%            ft2: fit structure
%
%
% May 2015
% Last modified 27/11/2015
% AstraZeneca, Alderley Park
% Juan Delgado-SanMartin, PhD Student
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set up defaults
    if ~iscell(WBImages);error('Input images as a Nx1 cell array');end
    if isempty(HypoxiaSetPoint);HypoxiaSetPoint = 20;end

    % Definition of nested function
    function buttonCallback()
        bp = 1;uiresume;
end


%% Commence of the algorithm
for l = 1:P.L
figure;bp = 0;
button1 = uicontrol('Style','pushbutton','String','Go on!','Position',[30 30 50 25],'Callback',@(src,evnt)buttonCallback);
button2 = uicontrol('Style','pushbutton','String','Make more','Position',[100 30 100 25],'Callback','uiresume');
Iaux = WBImages{l,1};Iaux = double(Iaux(:,:,1))./255;imshow(Iaux);title(D{l,1},'FontSize',14,'FontWeight','Bold')
hd = helpdlg('Select the bands you are interested in. Remember to select one for the background. To create another one, just drop it. When you are ready press: go on!','Select Western Blot');
uiwait(hd)
    % Select the bands of interest (by hand)
        k = 0;
    while bp == 0
        k = k + 1;
        h = imfreehand;Pos = h.getPosition;H(k) = {h};
        txt = text(Pos(1,1),Pos(1,2),num2str(k)); set(txt,'Color',rand(3,1),'FontSize',14,'FontWeight','Bold') 
        uiwait;clear Pos h txt
    end; 
    
        % Label them
        myData = cell(k,1);
        t = uitable('ColumnFormat',({'numeric'}),... 
                    'ColumnEditable', true,'Data',myData,'ColumnName','% of Oxygen');
        [nx,ny,~] = size(Iaux);
        aux = get(t,'Extent');aux(1:2) = [220,aux(4)];
        set(t,'Position',aux)
        uicontrol('Style','text','String','Write the amount of oxygen in the table. Select 0 for background',...
            'Position',[aux(1:2)+aux(3:4),100,50])
        uiwait
        
    % -------------------- Quantify -------------------------------      
    % Define oxygen
        O2 = cell2mat(get(t,'Data'));
        
        % Define the quantification
        for i = 1:length(H)
            V(i,1) = sum(Iaux(H{i}.createMask));
            VTotal(i,1) = sum(sum(H{i}.createMask)); 
        end;clear i
  
        % Re-arrange it
        uO2 = unique(O2);
        BG = sum(V(O2==0))./sum(VTotal(O2==0)); % background average/pixel
        Vcorr = V - VTotal.*BG;
        for r = 1:length(uO2);
        O2M(r,1) = {V(O2==uO2(r)) - VTotal(O2==uO2(r)).*BG};
        end;clear r
        
        % Substract background
        O2M = cell2mat(cellfun(@(x)mean(x),O2M,'UniformOutput',false));
        
        % Results
        ResultMean(l,:) = [{uO2},{O2M}];
        Result(l,:) = [{O2(O2~=0)},{Vcorr(O2~=0)-O2M(uO2==21)}];       
        
% Group/plot
        if l == 1;F = figure;ax = axes;xlabel('% Oxygen');ylabel('Expression HIF1\alpha');end
        set(0,'currentfigure',F);set(F,'currentaxes',ax);
        hold on
        plot(O2,Vcorr,'o','Color',P.colour(l,:))
        plot(ResultMean{l,1},ResultMean{l,2},'--*','Color',P.colour(l,:))
        hold off
        
        if l == 1;G = figure;ax2 = axes;xlabel('% Oxygen');ylabel('Expression HIF1\alpha');end
        set(0,'currentfigure',G);set(G,'currentaxes',ax2);
        hold on
        plot(O2(O2~=0),Vcorr(O2~=0)-O2M(uO2==21),'o','Color',P.colour(l,:))
        plot(ResultMean{l,1},ResultMean{l,2}-O2M(uO2==21),'--*','Color',P.colour(l,:))
        hold off
        
        clear O2M uO2
end;clear l

% Interpolate
set(0,'currentfigure',G);set(G,'currentaxes',ax2);
ft = fit(cell2mat(Result(:,1)),cell2mat(Result(:,2)),'exp1');
hold on;plot((0:1:21),ft.a.*exp(ft.b.*(0:1:21)),'k-')

%% Unit conversion adn refit
    O2UC = 1/100./P.K_HO2*1000;
    MaxHIF = ft.a.*exp(ft.b.*0);
    ft2 = fit(cell2mat(Result(:,1))*O2UC,cell2mat(Result(:,2))./MaxHIF*100,'exp1');


%% Calculate threshold of hypoxia
    % Oxygen concentration to express hypoxia  
    h_H = 1./ft2.b*log(HypoxiaSetPoint./ft2.a);   %mmol/L
    h_HX_O2 = 1./ft.b*log(HypoxiaSetPoint./ft.a);   % in % of oxygen in air
    

    %% Make final plot
%         xv = linspace(0,max(Result{2,1}*O2UC),100); 
%         figure;set(gca,'FontSize',12,'FontWeight','Bold')
%         hold on
%         plot(Result{1,1}.*O2UC,Result{1,2}/MaxHIF*100,'o','Color',P.colour(1,:))
%         plot(Result{2,1}.*O2UC,Result{2,2}/MaxHIF*100,'o','Color',P.colour(2,:))
%         plot(xv,ft2.a.*exp(ft2.b.*xv),'k-')
%         plot(h_H,HypoxiaSetPoint,'ko','MarkerSize',10)
%         plot(h_H,HypoxiaSetPoint,'kx','MarkerSize',10)
%         aux = ResultMean{1,1};aux2 = ResultMean{1,2};aux2 = aux2-aux2(aux==21);plot(aux(aux~=0).*O2UC,aux2(aux~=0)/MaxHIF*100,'*--','Color',P.colour(1,:))
%         aux = ResultMean{2,1};aux2 = ResultMean{2,2};aux2 = aux2-aux2(aux==21);plot(aux(aux~=0).*O2UC,aux2(aux~=0)/MaxHIF*100,'*--','Color',P.colour(2,:))
%         plot([0 h_H],[HypoxiaSetPoint HypoxiaSetPoint],'k--')
%         plot([h_H h_H],[-20 HypoxiaSetPoint],'k--')
%         ylim([-20 100]);xlim([0 .3])
%         xlabel('mmol/L');ylabel('% HIF1\alpha expression');
%         legend('Exp. 1','Exp. 2','Fit','SetPoint')  
%         hold off
end
