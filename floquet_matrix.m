
function [outVals] = floquet_matrix(k, Delta, d_1, d_2, Z, V)
% clear all;
% close all;
%nargin
if nargin < 6,
    %msgbox('Not enough constants provided');
    prompt = {'Please Enter a Value for "k"', 'Please Enter a Value for "Delta"', ...
             'Please Enter a Value for "delta_1"','Please Enter a Value for "delta_2"', ...
             'Please Enter a Value for "Initial Position"','Please Enter a Value for "Initial Velocity"'};
    dlg_title = 'Not Enough Parameters Provided';
    num_lines = 1;
    defaults = {'8.0554e6', '0', '26.6667', '40', '0', '0'};
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    answer = inputdlg(prompt, dlg_title, num_lines, defaults, options);
    
    k = str2double(answer(1,1)); %2*pi/780e-9;    %velocity normalization factor
    Delta = str2double(answer(2,1)); %0;   %detuning...this is a normalize value
    d_1 = str2double(answer(3,1)); %80e6/3e6;
    d_2 = str2double(answer(4,1)); %120e6/3e6;
    Z = str2double(answer(5,1)); %0;
    V = str2double(answer(6,1)); %0;
end;

W_1 = 2;
W_2 = 2;
W_3 = 2;
W_4 = 2;

flo = zeros(15,15);
flo(1,1) = -1i*(Delta-k*V);
flo(1,9) = -1i*W_1;
flo(1,10) = -1i*W_2;

flo(2,2) = 1i*W_1/2;
flo(2,3) = 1i*W_2/2;
flo(2,12) = -1i*conj(W_1)/2;
flo(2,13) = -1i*conj(W_2)/2;

flo(3,7) = 1i*W_1;
flo(3,8) = 1i*W_2;
flo(3,11) = 1i*(Delta-k*V);

flo(4,2) = -1i*(d_1+Delta-k*V);
flo(4,6) = 1i*conj(W_1)/2;

flo(5,1) = -1i/2*conj(W_3)*exp(-1i*2*k*Z);
flo(5,7) = 1i*d_1;
flo(5,11) = 1i/2*W_1/2;

flo(6,6) = -1i*W_3/2*exp(-1i*2*k*Z);
flo(6,14) = -1i*(d_1-Delta+k*V);

flo(7,4) = 1i*(d_1-Delta+k*V);
flo(7,6) = -1i*conj(W_3)/2*exp(1i*2*k*Z);

flo(8,1) = -1i*conj(W_1)/2;
flo(8,9) = -1i*d_1;
flo(8,11) = 1i*W_3/2*exp(1i*2*k*Z);

flo(9,6) = -1i*conj(W_1)/2;
flo(9,12) = 1i*(d_1+Delta-k*V);

flo(10,3) = -1i*(d_1+d_2+Delta-k*V);
flo(10,6) = -1i*conj(W_2)/2;

flo(11,1) = -1i/2*conj(W_4)*exp(-1i*2*k*Z);
flo(11,8) = -1i*(d_1+d_2);
flo(11,11) = 1i*W_2/2;

flo(12,6) = -1i*W_4/2*exp(-1i*2*k*Z);
flo(12,15) = -1i*(d_1+d_2-Delta+k*V);

flo(13,5) = 1i*(d_1+d_2+Delta-k*V);
flo(13,6) = 1i*conj(W_4)/2*exp(1i*2*k*Z);

flo(14,1) = -1i/2*conj(W_2);
flo(14,10) = 1i*(d_1+d_2);
flo(14,11) = 1i/2*W_4*exp(1i*2*k*Z);

flo(15,6) = -1i*W_2/2;
flo(15,13) = 1i*(d_1+d_2+Delta-k*V);

%      1  2  3   4               5  6                               7                              8  9          10               11 12                        13                             14 15  
ket = [0; 0; 0; -1i*conj(W_1)/2; 0; -1i*conj(W_3)/2*exp(-1i*2*k*Z); -1i*conj(W_3)/2*exp(1i*2*k*Z); 0; -1i*W_1/2; -1i*conj(W_2)/2; 0; -1i*W_4/2*exp(-1i*2*k*Z); -1i*conj(W_4)/2*exp(1i*2*k*Z); 0; -1i*W_2/2;];

outVals = flo\ket;
% 
% A = outVals(1,1);
% B = outVals(2,1);
% C = outVals(3,1);
% D = outVals(4,1);
% E = outVals(5,1);
% clear outVals;
% clear flo;
% clear ket;
% t = 100e-6;
% npts = 1e6;
% rho_12 = zeros(1, npts);
% 
% for i = (0:npts),
%     time = i*t;
%     phi_1 = d_1*time;
%     phi_2 = (d_1 + d_2)*time;
%     rho = A+B*exp(1i*phi_1)+C*exp(1i*phi_2)+D*exp(-1i*phi_1)+E*exp(-1i*phi_2);
%     rho_12(i+1) = real(rho);
% end;
% plot(rho_12(1:npts));
% clear rho_12;




