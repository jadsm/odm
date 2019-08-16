function [k_R1D, k_R2D] = PDEModule(HIFHM, Profiles, boundary,h_H,ft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDEModule: this module calculates the diffusion parameter of the 
% tumour-stroma model presented in the PhD thesis:
% "Mathematical models for heterogeneous preclinical cancers" by Juan Delgado San Martin
% sumbited for the degree of PhD in physics to the university of Aberdeen.
%
% This piece of work will be submitted to npj: systems biology journal under the name:
% "Dual tumour-stroma relationship explained with a cellular automaton" in 2016
%
% There is unrestricted license to use this script and modify it as long as the Author is acknowledged
% and either of the above publlications correctly cited.
% 
%
% SYNTAX:
% [k_R] = PDEModule(HIFHM, boundary,h_H,ft)
%
% INPUTS:
%      HIFHM - Gray image as a matrix (double)
%      Profiles - nx1 double array - profiles of oxygen 1D
%      boundary - struct containing fields:
%               B: nx2 double array containing coordinates of boundary
%     micronpixel: 1x1 floating point with the aspect ration in um/pixel
%             Ind: nx1 double array containing indices of sides of the
%             slide 1 filter, 2 right, 3 Air, 4 left
%        centroid: 2x1 double array containing the coordinates of centre
%      h_H - 1x1 float containing threshold of hypoxia
%      ft - 1x1 cfit structure (coming from the WBModule)
%      cmlayer - scalar with the average slice thickness in cm/layer 
%
% OUTPUTS
%      k_R1D - kr' apparent oxygen uptake rate by the cells calculated in one
%      dimension
%      k_R2D - kr' apparent oxygen uptake rate by the cells calculated in
%      two
%      dimensions
% Created: May 2015
% AstraZeneca, Alderley Park and Cambridge Institute
% Juan Delgado-SanMartin, PhD Student
%
% Modified: August 2015
% Modified: November 2015
% Modified: February 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Set up
        Par.n = 120;
        Par.ref = 2;
        Par.h_H = h_H;
        Par.cmlayer = boundary.cmlayer;
        [Par.nx,Par.ny] = size(HIFHM);

        % Initial condition
        kro = 50.33; % cm-1
    %% Call the optimiser
           opts = optimoptions(@fmincon,'Algorithm','sqp','DerivativeCheck','off');
              ms = MultiStart('UseParallel','always','Display','iter');
%     % -------------- 1D diffusion --------------------------
         problem = createOptimProblem('fmincon','objective', ...
              @(kr)ObjFcn1D(kr,Profiles,Par.cmlayer,ft),'x0',kro,'lb',[0], ...
              'ub',[inf],'options',opts);
                  [k_R1D,f] = run(ms,problem,50);
%           ProfileCalc = Profilemodel(k_R1D,Par.cmlayer,Profiles,ft);
    % -------------- 2D diffusion --------------------------  
         problem = createOptimProblem('fmincon','objective', ...
              @(kr)ObjFcn(kr,HIFHM,boundary,Par,ft),'x0',kro,'lb',[0], ...
              'ub',[inf],'options',opts);
              [k_R2D,f] = run(ms,problem,50);

end
function OF = ObjFcn (kr,HIFHM,boundary,Par,ft)
% Calculate the PDE profile
tic
Mask = PDEmodel(kr,boundary,Par,HIFHM,ft);
toc

% Calculate the cost function
OF = Mask./100-HIFHM;
OF = sum(OF(:).^2);

% Plot with transparency
% delete(get(gcf,'Children'))
% A = zeros([size(Mask) 3]);
% A(:,:,1) = Mask;
% A(:,:,2) = HIFHM;
% imshow(A)
% title(['(k_R/D)^{1/2} = ' num2str(kr,'%.2f') ' cm^{-1}'])
end
function OF = ObjFcn1D (kr,Profiles,cmlayer,ft)
% Calculate the linear model
% tic
ProfileCalc = Profilemodel(kr,cmlayer,Profiles,ft);
% % toc
% Calculate the cost function
profaux = sum(Profiles(:,1:3),2)./Profiles(:,5);profaux(isnan(profaux)) = 0;
OF = ProfileCalc'/100-profaux;
OF = sum(OF(:).^2);
% figure;semilogy([ProfileCalc'./100,profaux])
end
function ProfileHIF = Profilemodel(kr,cmlayer,Profiles,ft);
% Oxygen levels
O2air = 1.19;
O2filter = 0.28;

% Radius
r = linspace(0,cmlayer*size(Profiles,1)/10,size(Profiles,1));

% Profiles 
ProfileO2 = -(O2air*exp(2*kr) - O2air*exp(2*kr*r) - O2filter*exp(kr) + O2filter*exp(2*kr*r)*exp(kr))./(exp(kr*r) - exp(kr*r)*exp(2*kr));
ProfileHIF = ft.a*exp(ft.b*ProfileO2);
end
function Mask = PDEmodel(kr,boundary,Par,HIFHM,ft)

% Decomposed geometry for PDE
[pdegd,sf,ns] = MyGeometry(boundary,Par.n);
gd = decsg(pdegd,sf,ns);%decsg(pdegd,sf,ns);

% gd = @pdemygeom;
% figure;
% pdegplot(gd,'edgeLabels', 'on','subdomainLabels','on')
% axis equal
%% --------------------------Create Mesh---------------------------------
[p,e,t] = initmesh(gd,'Hmax',100,'Jiggle','off');
% pdeplot(p,e,t)
% pdegplot('lshapeg','edgeLabels','on')
%  Extrafine triangulation 
for ll = 1:Par.ref
[p,e,t] = refinemesh(gd,p,e,t);
end
p = jigglemesh(p,e,t);

%% ---------------------Set PDE conditions-------------------------------
% Set the time steps for the parabolic solver to 50 stepsfrom time 0 to time 1.
% Solve the parabolic PDE.
b = @(p,e,u,time)pdebound(p,e,u,time,boundary); % boundary conditions
% a = @(p,t,u,time)acoeff(p,t,u,time,(kr.*boundary.micronpixel*1e-0).^2);% Oxygen Uptake Rate
a = (kr.*boundary.micronpixel*1e-4).^2;
f = 0;% Independent coefficient
c = @(p,t,u,time)ccoeff(p,t,u,time);%'Proba2(x,y,u,ux,uy,t,sd)';%'(sd.^2)';%@(p,t,u,time)DiffusionCoefficient(p,t,u,time,sdID); % Diffusion coefficient
u = pdenonlin(b,p,e,t,c,a,f); 

%% Create mask
% Interpolation for the centres of the triangles
[x,y,uintrp,~,~,~] = vartransform(p,t,u);

% Convert O2 in HIF expression
muint = ft.a*exp(ft.b*uintrp(end,:)); % At the end (steady state)
mu = ft.a*exp(ft.b*u);

% interpolation on the grid
[ma, mb] = size(HIFHM);
F = scatteredInterpolant(x',y',muint(end,:)');

[xq,yq] = meshgrid(1:mb,1:ma);
%  in = inpolygon(xq,yq,boundary.B(:,1),boundary.B(:,2));
Mask = F(xq,yq);
[~,BW] = roifill(Mask,boundary.B(:,1),boundary.B(:,2));
Mask(~BW) = 0;

% mesh(xq,yq,vq); hold on; plot3(x,y,uintrp,'o'); hold off
%% plot
% plotme(p,e,t,mu,HIFHM)
end
function [pdegd,sf,ns] = MyGeometry(boundary,n)
%% MyGeometry: this script calculates the geometry programatically
% This script has been done generically, for which any sort of geometry would be allowed
% 
% INPUTS: 
%       B is a
%       Bim
%       n
%       A
% OUTPUTS:
%       pdegd is the geometry matrix - here all the objects are defined
%       ns is the labels matrix. The objects will be numbered automatically
%       preceeded by the letters TS (tumour slice). Maximum number of
%       objects is 999.
%       sf is the relationship matrix across all the objects
%

% set up the boundary for the 
[B,ib,~] = unique(round(boundary.B),'rows');
% sort by angle (clockwise)
theta = atan2((B(:,1)-boundary.centroid(1)),(B(:,2)-boundary.centroid(2)));
[~,ii] = sort(theta);B = B(ii,:);
% Ind = boundary.Ind;


% set up geometry matrix
%         pdegd = zeros(n*2+2,length(B));% Allocate geometry matrix
 pdegd = zeros(n*2+2,1);% Allocate geometry matrix
%         hold on
        pdegd(1,:) = 2;
        ns = ['TS1']';

% Calculate boundaries and geometries
    ll = size(B,1); 
    db = downsample(B,ceil(ll/n));
%     dInd = boundary.Ind(1:ceil(ll/n):end);
    l = length(db); 
    pdegd(2,1) = l; 
    pdegd(3:l*2+2,1) = db(:);
     
%     plot(bb(:,2), bb(:,1), 'r','LineWidth',2);

% Check coherence
    ns = ns(:,csgchk(pdegd)==0);% UNCOMMENT
    pdegd = pdegd(:,csgchk(pdegd)==0);

% Set up additive rules
sf = ns(:,1)';
for ii = 1:size(ns,2);sf = [sf '+' ns(:,ii)'];end% UNCOMMENT

end
function [qmatrix,gmatrix,hmatrix,rmatrix] = pdebound(p,e,u,time,boundary)
%% pdebound
% This function computes the boundary conditions of the PDE problem for 
% oxygen diffusion (Fick equation).
% AstraZeneca, Alderley Park
% 17/04/2014 Juan Delgado, PhD Student
% Allocation of the system
N = 1; % Size of the system
ne = size(e,2); % number of edges
qmatrix = zeros(N^2,ne);
gmatrix = zeros(N,ne);
hmatrix = zeros(N^2,2*ne);
rmatrix = zeros(N,2*ne);

   % Find the coordinates
    x1 = p(1,e(1,:)); % x at first point in segment
    x2 = p(1,e(2,:)); % x at second point in segment
    xm = (x1 + x2)/2; % x at segment midpoint
    y1 = p(2,e(1,:)); % y at first point in segment
    y2 = p(2,e(2,:)); % y at second point in segment
    ym = (y1 + y2)/2; % y at segment midpoint

    % Find he top and bottom sides of the slide
    
            M = 1;
      % Matrix h with dimensions (N^2x2ne) - Slope Dirichlet
            hmatrix(:) = 1;

      % Matrix r with dimensions (Nx2ne) - Coefficients Dirichlet
%             rk = [159.6 100 47 100]; %mmHg
%              rk = [1.19 1 0.28 1];
            rk = [0.28 0 1.19 0];
         if ne>1
                rmatrix(repmat(xm==min(xm),1,2)) = rk(2);
                rmatrix(repmat(xm==max(xm),1,2)) = rk(4);
                rmatrix(repmat(ym>boundary.centroid(2)&xm~=min(xm)&xm~=max(xm),1,2)) = rk(1);
                rmatrix(repmat(ym<boundary.centroid(2)&xm~=min(xm)&xm~=max(xm),1,2)) = rk(3);  
            else
              rmatrix(:) = rk(1);
            end




     % Matrix q with dimensions (N^2xne) - Slope Neumann (Set 0 for pure Neumann BC)
            qmatrix(:) = 0;
            
      % Matrix g with dimensions (Nxne) - Coefficients Dirichlet
            gmatrix(:) = 0;

end
function [c] = ccoeff (p,t,u,time)
%Get the coefficients
 [x,y,u,ux,uy,sd] = vartransform(p,t,u);
 
 % Output function
c = ones(1,length(t));
c = 1;
% c = 0.0068; % Epithelium
% c(~sdID) = 1;
end
function [a] = acoeff (p,t,u,time,kr)
%Get the coefficients
 [x,y,u,ux,uy,sd] = vartransform(p,t,u);

% Output of the function
a = zeros(1,length(t));
% a = 1.92; % Epithelium
a = kr;
end
function plotme(p,e,t,u,Mask)

% create axes
    if isempty(get(gcf,'currentaxes'));
            subplot(211);
            title('% HIF1\alpha','FontSize',12,'FontWeight','Bold')
            xlabel('Length mm','FontSize',12,'FontWeight','Bold')
            ylabel('Width mm','FontSize',12,'FontWeight','Bold')
            set(gca,'FontSize',12,'FontWeight','Bold');
            subplot(212);
    end
ax = get(gcf,'Children');

% plots
set(gcf,'currentaxes',ax(end))
pdeplot(p,e,t,'xydata',u(:,end),'colormap','jet','mesh', 'off')


set(gcf,'currentaxes',ax(end-1))
imshow(flip(Mask,1))
drawnow

end
function [x,y,uintrp,ux,uy,sd] = vartransform(p,t,u)

        if nargin == 3;h = 1;
        elseif nargin==2;h = 0;
        elseif nargin<=1;error('Please, introduce at least position and triangulation indexes matrices p and t');
        end

% This fucntion modifies the function form input to the string form input
% Triangle point indices
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);

% Find centroids of triangles
x =(p(1,it1)+p(1,it2)+p(1,it3))/3;
y =(p(2,it1)+p(2,it2)+p(2,it3))/3;
if h == 1
    
        % Interpolate the function value
        uintrp = pdeintrp(p,t,u); % Interpolated values at centroids
    try
        % Calculate gradients
        [ux,uy] = pdegrad(p,t,u); % Approximate derivatives
    catch
    ux = [];uy = [];
    disp('Derivatives not calculated!')
    end
else
    uintrp = [];
    ux = [];uy = [];    
end
% Identify subdomain
sd = t(4,:);
end

