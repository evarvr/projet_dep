function [t,x] = simulate1(p,temps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATE: 
% Simulation des concentrations pour un modele bi-bi ping-pong 
% dans la phase I de la reaction (entre 0 et tsw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global texp  xexp
global tref xref
global tsw Itsw xsw
global E0 Es0 Ed0
global modele

p = abs(p);


B0  = xexp(1,1);
AR0 = xexp(1,2);

%------------------------------------------
% Conditions initiales
%------------------------------------------
%EA0 = p(end)*E0; 	% Initialisation non nulle pour demarrage rapide.

% 	B   AR   R   AB ABA   Ed
x0 = [ B0  AR0   0   0   0    Ed0  ];
if strcmp(modele,'bbpp')
  EA0 = 0;
  x0  = [ x0 EA0];
elseif strcmp(modele,'bbo')
  EARB0  = 0;
  EARAB0 = 0;
  x0 = [ x0 EARB0 EARAB0];
elseif strcmp(modele,'bbppo')
  EA0    = 0;
  EARAB0 = 0;
  x0 = [ x0 EA0 EARAB0];
else
  error('nom de modele incorrect')
end

options = odeset;

if nargin>1 
  %options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
  [t,x]=ode23s(modele,temps ,x0,options, p);
else
  %options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
  [t,x]=ode23s(modele,[0 tsw] ,x0,options, p);
end

% etat x au moment du switch pour reprise modele 2
xsw = x(end,1:6);

