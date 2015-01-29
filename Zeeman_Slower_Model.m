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
T = 300;  %in units of degrees K
K_b = 1.38e-23;    %in units of J/K
M = 87 * 1.66e-27;	%in units of kg
V_mp = sqrt(2*K_b*T/M); %Velocity, most probable in units of meters/second
%P_v = 2 * V^3/V_mp^4 * exp(-(V^2/V_mp^2));
Del = 1;    %
Omega = 1;  %Rabi frequency in units of ...
H_bar = 1.05457148e-34;  %in units of ...
Lambda = 780.24; %laser wave length in units of meters...
K = 2*pi/Lambda;    %wave vector...
Beta = 2*pi*3e6;    %in units of MHz
B_0 = 1200;    %in units of Gauss
L_0 = 1;  %in units of meters
L = 1; %in units of meters
B = B_0 * sqrt(1 - L/L_0);  %in units of Gauss
Gamma = 2*pi*1e6;   %in units of MHz/Gauss

% =========================================================================

close all
%Create seed for random number generator.
%s = rand(seed);
nbins = 100;  %for testing
%Create random number generator uniformly distributed for tau = [0,1).
Tau = rand(1e4, 1);
%hist(Tau, nbins);

%Calculate the Velocity distribution from Tau
V = V_mp * gammaincinv(Tau,2);%erfinv(2 * Tau - 1);

%Create and view the histogram to verify that the is correct...this is for
%testing
fid1 = figure;
hist(V, nbins);
%std(V); %for testing purposes
%V_mp/sqrt(2)   %for testing purposes

%Apply a force to the velocity distribution (V).
Force = -9e-20; %hbar*k/2Tau
Delta_t = 1e-11;  %the time that the force is supplied
T_step = 1e-8; %step size
Step = 1e4; %number of steps
V_new = zeros(size(V));  %empty array for the new/final velocity distributuition
boundary = std(V);  %velocity value where the force is applied
%set(gca, 'nextplot', 'replacechildren');
for j = 1:Step
    Time = Delta_t + j*T_step;
    for i = 1:size(V)
        V_temp = V(i);
        if ((V_temp >= -boundary)&&(V_temp <= boundary)),
            calculation = V_temp + (Force/M) * Time;
            V_new(i) = calculation;
        else
            V_new(i) = V_temp;
        end;
    end;
    %2
    %hist(V_new, nbins); 
    %F(j) = getframe;
end;
 
%movie(F, 1);
fid2 = figure;
hist(V_new, nbins);

%Gaussian Position Probability distribution...will need this for another
%application...
%