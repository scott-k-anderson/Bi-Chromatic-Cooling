function Omega = Rabi(t,z, deltaW);
c = 299792458;
deltaK = deltaW/c;
Omega = 1 + cos(deltaW * t - deltaK * z + phase);