%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODM: this script runs the program for the Oxygen Driven Model (ODM) as described in the paper 
% "Oxygen Driven Tumour Growth Model: A pathophysiology-relevant
% mathematical approach",(2015), JA Delgado-SanMartin, JI Hare, APS Moura, JWT Yates
%
% This script has been used to calculate all the parameter values in the
% manuscript.
%
% Accepted in PLoS Computational Biology
%
% Created by: Juan A. Delgado-SanMartin, M.Eng.
% juan.ads.delgado@astrazeneca.com
% AstraZeneca, Alderley Park, Cheshire, UK, Summer 2013
% Last modified: June 2015
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = ODM ()
%% Transformation of parameters
% Parameters
kp = .0033;
kr = 3.3;
Kh = 25.2;
Kn = 2;
 
% other variables
n = 19; % Number of shells
Po2maxo = 100; % mmHg Maximum oxygen at periphery of tumour
Vto = 0.3; % cm3 intial tumour volume
fh = 0.1; % initial fraction of hypoxia
fn = 0.1; % initial fraction of necrosis
% Volumes
V = ones(n,1).*Vto/n; % layer volume
Vh = V*fh; % initial hypoxia volume
Vn = V*fn; % initial necrosis volume
V = [V;Vh;Vn];

% Specify algebraic ecuations
M = eye(length(V)); M(n+1:end,n+1:end) = 0;
options = odeset('Mass',M);

% Model call
[T,Y] = ode15s(@(t,v)ODM(t,v,kp,kr,Kh,Kn,n,Po2maxo,Vto),[0 30],V,options) ;

% Plot
F = [sum(Y(:,1:n),2) sum(Y(:,n+1:n*2),2) sum(Y(:,2*n+1:n*3),2)];
plot(T,F);legend('Tumour','Hypoxia','Necrotic','Location','Best')
xlabel('days');ylabel('Volume cm^3')

save F2 F

function dYS = ODM(~,v,kp,kr,Kh,Kn,n,Po2maxo,Vto)  
% Volumes
V = v(1:n);
Vh = v(1+n:2*n);
Vn = v(1+2*n:end);
Vt = sum(V(1:n)); % total volume

% Total Radius and Volume
Rend = (Vt*3/4/pi).^(1/3);
 
%% Model Calculation
% Vector Allocation
Po2 = zeros(n,1);dY = zeros(n,1);
 
% Assumptions
%       1: Area/Volume reduction in O2 accessibility
        Po2max = Po2maxo*(Vto/Vt)^(2/3);
%       2: Discretisation. The Oxygen tension drop is constant
        DeltaPo2 = (1-1./cosh(Rend.*kr))./n;
 
% Definition of the oxygen profile
Po2(:,1)= Po2max.*(DeltaPo2.*(1:n)'+1./cosh(Rend.*kr));
 
% Definition of the ODE (tumour growth)
dY (:,1) = kp.*Po2.*V;
 
% Definition of the Algebraic system (Hypoxic/Necrotic count)
 N = Po2./(Po2+Kn).*V - Vn;
 H = Po2./(Po2+Kh).*V - Vh;
 
%% -----------------------Output ---------------------------------------- 
 dYS = [dY;N;H]; % output of the DAE system

