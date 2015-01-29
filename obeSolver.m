%
%   obeSolver.m
%
%   Written by: Scott Anderson
%   Written on: July 31, 2014
%
%   This program was written to solve the equations that comprise the Optical Bloch
%   as written in D. Chieda's thesis.

close all
clear all
%Create seed for random number generator.
nbins = 100;  %for testing
%Create random number generator uniformly distributed for tau = [0,1).
atoms = 1e0;
Tau = rand(atoms, 1);

%Calculate the Velocity distribution from Tau
T = 4.5e2;  %in units of degrees K
K_b = 1.38e-23;    %in units of J/K
mass = 87 * 1.66e-27;	%in units of kg
V_mp = sqrt(2*K_b*T/mass); %Velocity, most probable in units of meters/second
i_Vel = V_mp * gammaincinv(Tau,2); %this is the velocity distribution

norm = 2*pi*3e6;    %normalization factor = to gamma
gamma = 1;  %should be 2pi*3e6, but it is normalized 
Delta = -123;   %detuning...this is a normalize value
deltaW = 2*pi*100e6/norm;   %delta omega...normalized
c = 299762458;  %speed of light
k = 2*pi/780e-9;    %velocity normalization factor
i_Vel =  i_Vel*k/norm;  %normalized
deltaK = norm*deltaW/(k*c);  %normalized
phase = 0;  %Rabi frequency phase
hBar = 1.055e-34;
amp = 8;    %Rabi frequency amplitude
f_Vel = zeros(size(i_Vel)); %initialized final velocity array
pos = 2;%rand(atoms, 1);    %Position of th atom.  Need to change it to a position distribution
pos = pos*k;    %normalized
out = zeros(size(i_Vel));
const = [mass, hBar, pos, Delta, deltaW, deltaK, phase, amp, k, norm, out]; %constants that need to be passed to OBE function
options = odeset('RelTol',1e-3,'AbsTol',1e-3);  %accuracy options for ode45
tspan = [0 100e-6]*norm;    %time span = cooling time of BCF scaled by gamma
%error('stop here');

%structure to calculate the final velocity and position for each atom in
%the velocity and position distributions

flo = floquet_matrix(8.0554e6, 0, 26.6667, 40, 0, 0);

for i = 1:size(i_Vel)
    iCond = [0,0,-1,i_Vel(i),pos];   %pass in the initial conditions 
    [t,y] = ode45('OBE', tspan, iCond, options, const);
    f_Vel(i) = y(length(t),4);
end;

if 0
    subplot(2,1,1);
    plot(t*1e6/norm,y(:,4)*norm/k);
    xlabel('time (us)');
    title('Velocity');
    subplot(2,1,2);
    plot(t*1e6/norm,y(:,5)/k);
    xlabel('time (us)');
    title('Position');
end;  

    A = flo(1,1);
    B = flo(2,1);
    C = flo(3,1);
    D = flo(4,1);
    E = flo(5,1);
    d_1 = 26.6667;
    d_2 = 40;
    t = 100e-6;
    npts = size(y,1);
    rho_12 = zeros(1, npts);
    
    for i = (0:npts),
        time = i*t;
        phi_1 = d_1*time;
        phi_2 = (d_1 + d_2)*time;
        rho = A+B*exp(1i*phi_1)+C*exp(1i*phi_2)+D*exp(-1i*phi_1)+E*exp(-1i*phi_2);
        rho_12(i+1) = real(rho);
    end;

    figure;
    subplot(2,2,1);
    plot(t*1e6/norm, rho_12)   %t*1e6/norm,real(y(:,1)), t*1e6/norm,imag(y(:,1)), 
    xlabel('time (us)');
    title('U');
if 0
    subplot(2,2,2);
    plot(t*1e6/norm,real(y(:,2)),t*1e6/norm,imag(y(:,2)))
    xlabel('time (us)');
    title('V');
    subplot(2,2,3);
    plot(t*1e6/norm,real(y(:,3)),t*1e6/norm,imag(y(:,3)))
    xlabel('time (us)');
    title('W');
end;
    
if 0
    npts = length(t);
    tmax = t(npts)/norm;
    time = (1:npts)*tmax/npts*1e6;
    f = (-npts/2+1:npts/2)/tmax;
    u_new = interp1(t,y(:,1),time, 'spline');

    F = fftshift(fft(u_new));

    figure;
    plot(f, F.*conj(F));
    axis([-3e7,3e7,0,inf]);
    title('FFT Plot of U');
    %xlable('Frequency - Uniformly Spaced');    
end;

if 0
    figure;
    subplot(2,1,1);
    hist(f_Vel*norm/k, 100)  %histogram for the final velocity distribution
    title('Final Velocity Histogram');
    subplot(2,1,2);
    hist(i_Vel*norm/k,100) %histogram for the initial velocity distripution
    title('Initial Velocity Histogram');
end;
