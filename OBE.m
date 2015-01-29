function yDot = OBE(t, y, iflag, const)

%below are the constants needed to solve the system of 5 equations
%pos = const(3);   %already normalized
Delta = const(4);   %already normalized
deltaW = const(5);  %delta omega already normalized
deltaK = const(6);  %delta k already normalized
phase = const(7);   %a dimentionless value
amp = const(8);     %a dimentionless value
k = const(9);   %already normalized.  2*pi/780e-9;
Rabi = amp*(1 + cos(deltaW*t - deltaK*y(5) + phase));   %***this is the Rabi frequency

yDot = zeros(5,1);  %initialize the array
%These equations (yDot1,2,&3) are the Optical Bloch Equations taken from
%D. Cheida's thesis
yDot(1) = -(1/2)*y(1) - Delta*y(2) - imag(Rabi)*y(3); %should be normalized 8/18/14
yDot(2) = Delta*y(1) - (1/2)*y(2) - real(Rabi)*y(3);  %should be normalized 8/18/14 
yDot(3) = imag(Rabi)*y(1) + real(Rabi)*y(2) - (y(3) + 1);   %should be normalized 8/18/14

gradRabi = amp*deltaK*sin(deltaW*t-deltaK*y(5));   %Rabi frequency gradiant...should be normalized 8/18/14
hBar = const(2);
gamma = const(10);  %2*pi*3e6;
F = (-hBar*k*gamma)*(y(1)*gradRabi - y(2)*0);  %Force equation taken from Cheida's Thesis...should be normalized 8/18/14
a = F/const(1); % const(1) = mass.  Newtonian equation for Acceleration...should be normalized 8/18/14
yDot(4) = a;    %assign the acceleration...should be normalized 8/18/14
yDot(5) = y(4) + a*t;   %Newtonian equation to calculate the position.

%const(11) = floquet_matrix(8.0554e6, Delta, 26.6667, 40, yDot(4), yDot(5));
%y(6) = outVals;

close all;