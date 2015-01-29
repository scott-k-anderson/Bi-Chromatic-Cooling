function yDot = rigid(t, y, iflag, w, A, B);
yDot = zeros(2,1);
yDot(1) = w*y(2);%A*cos(w*t)-B*sin(w*t); 
yDot(2) = -w*y(1);%A*sin(w*t)+B*cos(w*t);
close all