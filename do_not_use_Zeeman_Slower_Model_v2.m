%
%   Zeeman_Slower_Model.m
%
%   Written by: Scott Anderson
%   Written on: June 23, 2014
%
%   This program was written to model the effects of a Zeeman Slower.
%
%   Below are the global variables
% =========================================================================
%P_v = 2 * V^3/V_mp^4 * exp(-(V^2/V_mp^2));
hBar = 1.05457148e-34;  %in units of ...
gamma = 1;%2*pi*3e6/2*pi*3e6;   %***Scale every frequency by this value.***
% =========================================================================

% close all
%Create random number generator uniformly distributed for tau = [0,1).
atoms = 1;
Tau = rand(atoms, 1);

%Calculate the Tau distribution from Tau
T = 300;  %in units of degrees K
K_b = 1.38e-23;    %in units of J/K
mass = 87 * 1.66e-27;	%in units of kg
V_mp = sqrt(2*K_b*T/mass); %Velocity, most probable in units of meters/second
V = V_mp * erfinv(2 * Tau - 1);

options = odeset('RelTol',1e-3,'AbsTol',1e-3);
tspan = [0,100e-6]*2*pi*3e6;
pos = 0;
iCond = [0,0,1,V(1),pos];

norm = 2*pi*3e6;
Delta = 0/norm;
deltaW = 2*pi*100e6/norm;
c = 299762458;
deltaK = deltaW/(c*norm);
phase = 0;
amp = 2;
const = [mass, hBar, gamma, Delta, deltaW, deltaK, phase, amp]

[t,y] = ode45('OBE', tspan, iCond, options, const);
%pause;
fid2 = figure;
plot(t,y)