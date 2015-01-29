function b = gammaincinv(x,a,tail)
%MATLAB Code Generation Library Function

%   Copyright 1984-2011 The MathWorks, Inc.
%#codegen

eml_invariant(nargin >= 2, ...
    eml_message('Coder:toolbox:gammaincinv_1'), ...
    'IfNotConst','Fail');
eml_invariant(isreal(x) && isreal(a) && ...
    isa(x,'float') && isa(a,'float') && ...
    ~issparse(x) && ~issparse(a), ...
    eml_message('Coder:toolbox:gammaincinv_2'), ...
    'IfNotConst','Fail');
if nargin < 3
    upper = false;
else
    eml_prefer_const(tail);
    upper = strcmp(tail, 'upper');
    eml_invariant(upper || strcmp(tail,'lower'), ...
        eml_message('Coder:MATLAB:gammainc_InvalidTailArg'));
end
b = eml_scalexp_alloc(complex(real(eml_scalar_eg(x,a))),x,a);
scalar_a = eml_const(eml_is_const(isscalar(a)) && isscalar(a));
if scalar_a
    ak = double(a(1));
    la = log(ak);
    lga = gammaln(ak);
    lgap1 = gammaln(ak+1);
end
for k = 1:eml_numel(b)
    xk = double(eml_scalexp_subsref(x,k));
    if ~scalar_a
        ak = double(eml_scalexp_subsref(a,k));
        la = log(ak);
        lga = gammaln(ak);
        lgap1 = gammaln(ak+1);
    end
    b(k) = eml_gammaincinv(xk,ak,la,lga,lgap1,upper);
end

%--------------------------------------------------------------------------

function rval = eml_gammaincinv(p,a,la,lga,lgap1,upper)
% Compute the inverse incomplete gamma function, lower or upper tail.
%
% Uses Halley's rational method to find a root of p = gammainc(x,a).  Starting
% values are based on Best and Roberts, Algorithm AS 91, Appl. Statist. (1975)
% Vol.24, p.385-388 and incorporates the suggested changes in AS R85 (vol.40(1),
% pp.233-5, 1991).
eml_prefer_const(upper);
% if upper
%     tail = 1;
% else
%     tail = 0;
% end
LOG2 = 0.693147180559945;
pIn = p;
ftol = 100*eps;
ztol = eps;
c1 = 0.01;
c2 = 0.222222;
c3 = 0.32;
c4 = 0.4;
c5 = 1.24;
c6 = 2.2;
c7 = 4.67;
c8 =  6.66;
c9 = 6.73;
c10 = 13.32;
% Check for trivial cases.
if upper
    if ~(a > 0)  % including a == NaN
        if a == 0
            if (0 <= p) && (p <= 1)
                rval = 0;
                return % gammaincinv(p,0,'upper') = 0, including p == 0;
            else % p == NaN
                eml_error('Coder:MATLAB:gammaincinv_YOutOfRange');
            end
        else
            eml_error('Coder:MATLAB:gammaincinv_NegativeArg');
        end
    elseif ~((0 < p) && (p < 1))  % including p == NaN
        if (p == 0)
            rval = eml_guarded_inf;
            return % gammaincinv(0,a,'upper') = Inf
        elseif p == 1
            rval = 0;
            return % gammaincinv(1,a,'upper') = 0 for a > 0
        else
            eml_error('Coder:MATLAB:gammaincinv_YOutOfRange');
        end
    elseif a == eml_guarded_inf
        if p == 0
            rval = eml_guarded_inf;
        else
            rval = 0;
        end
        return
    end
else
    if ~(a > 0)  % including a == NaN
        if a == 0
            if (0 <= p) && (p <= 1)
                rval = 0;
                return% gammaincinv(p,0,'lower') = 0, including p == 1
            else % p == NaN
                eml_error('Coder:MATLAB:gammaincinv_YOutOfRange');
            end
        else
            eml_error('Coder:MATLAB:gammaincinv_NegativeArg');
        end
    elseif ~((0 < p) && (p < 1))  % including p == NaN
        if p == 0
            rval = 0;
            return % gammaincinv(0,a,'lower') = 0
        elseif p == 1
            rval = eml_guarded_inf;
            return % gammaincinv(1,a,'lower') = Inf for a > 0
        else
            eml_error('Coder:MATLAB:gammaincinv_YOutOfRange');
        end
    elseif a == eml_guarded_inf
        if p == 0
            rval = 0;
        else
            rval = eml_guarded_inf;
        end
        return
    end
end
% If gammaln(a) has overflowed, ...
if lga == eml_guarded_inf
    rval = a;
    return
end
am1 = a - 1;
nu = 2*a;
% Work with the smaller of the two tails, ...
if p > 0.5
    p = 1 - p; % no round-off here, p is > 0.5
    % tail = 1 - tail;
    upper = ~upper;
end
% Use the lower tail for starting value calculation, keeping away from 1
if upper
    pLower = 1 - p;
    if pLower == 1
        pLower = 1 - ftol;
    end
    log1mpLower = log(p); % log(1-pLower)
else
    pLower = p;
    u = 1-pLower;
    if u ~= 1 % log(1-pLower)
        log1mpLower = log(u)*(-pLower/(u-1));
    else
        log1mpLower = -pLower;
    end
end
% Starting approximation for small chi-squared.
if nu < -c5*log(pLower)
    chi2 = power(pLower*exp(lgap1 + a*LOG2), 1/a);
    if chi2 < 100*realmin
        chi2 = 100*realmin;
    end
    % Starting approximation for v less than or equal to 0.32.
elseif nu <= c3
    chi2 = c4;
    maxIter = cast(200,eml_index_class);
    i = 0;
    while i < maxIter
        oldz = chi2;
        p1 = 1 + chi2*(c7+chi2);
        p2 = chi2*(c9 + chi2*(c8 + chi2));
        t = -.5 + (c7 + 2*chi2)/p1 - (c9 + chi2*(c10 + 3*chi2))/p2;
        chi2 = chi2 - (1 - exp(log1mpLower + lga + .5*chi2 + am1*LOG2)*p2/p1)/t;
        if abs(oldz - chi2) < c1*chi2
            break
        end
        i = i + 1;
    end
    % assert(i < maxIter);
else
    z = PHIinv(pLower);
    % Starting approximation using Wilson and Hilferty estimate.
    p1 = c2/nu;
    chi2 = z*sqrt(p1) + 1 - p1;
    chi2 = nu*chi2*chi2*chi2;
    % Starting approximation for p tending to 1.
    if chi2 > c6*nu + 6
        chi2 = -2*(log1mpLower - am1*log(.5*chi2) + lga);
    end
end
z = .5*chi2;
% Newton's method: new x = x - f/f'
% Halley's rational method: new x = x - 1/(f'/f - (f'/f')/2)
% This function: f = gammainc(x,a) - p, f' = x^(a-1)*exp(-x)/gamma(a), f' = ((a-1)/x - 1)*f'
oldf = eml_guarded_inf;
oldz = eml_guarded_inf;
maxIter = cast(1000,eml_index_class);
% Impose a lower limit of realmin on the convergence tolerance for the function value
if p > 1.0021e-294
    ftol1 = ftol*p;
else
    ftol1 = ftol*1.0021e-294;
end
if upper
    sgn = -1;
else
    sgn = 1;
end
zlo = 0;
zhi = realmax;
i = 0;
while i < maxIter
    f = sgn*(eml_gammainc(z,a,la,lgap1,upper) - p);
    % Halve the step if it overshot the root badly.
    if (f*oldf < 0) && (abs(oldf) <= abs(f))
        z = 0.5*z + 0.5*oldz;
        f = sgn*(eml_gammainc(z,a,la,lgap1,upper) - p);
    end
    % Update the brackets around the root.
    if f > 0
        zhi = z;
    else
        zlo = z;
    end
    % If (gammainc(z,a) - p) is small enough, exit. ...
    if abs(f) < ftol1
        break
    end
    if abs(z-oldz) < ztol*z+realmin
        break
    end
    oldz = z; oldf = f;
    % Use Halley's method.
    if i < 500
        df = exp((a-1)*log(z) - z - lga);
        z = z*(1 - f/(z*df + f*(z+1-a)/2)); 
        % Make sure we don't step out of the brackets.
        if z <= zlo
            % The lower bracket might be at z==0.  ...
            if zlo == 0
                if abs(double(upper)-p) < ...
                        abs(eml_gammainc(realmin,a,la,lgap1,upper)-p)
                    z = 0;
                    break
                else
                    zlo = realmin;
                end
            end
            z = 0.99*zlo + 0.01*zhi; % zhi always finite
        elseif z >= zhi
            % The upper bracket is always < Inf, ...
            z = 0.01*zlo + 0.99*zhi;
        end
        % Reasons for Halley's method not converging: ...
    else
        % If zlo is still small, ...
        if 1e8*zlo < zhi % typically when zhi is near realmax
            oldz = 1e8*zlo;
            oldf = sgn*(eml_gammainc(oldz,a,la,lgap1,upper) - p);
            if oldf > 0
                zhi = oldz;
            else
                zlo = oldz;
            end
            % Resetting oldz and oldf prevents the step halving this time ...
        end
        z = 0.5*zlo + 0.5*zhi; % zhi always finite
    end    
    i = i + 1;
end
if i >= maxIter
    if ~coder.setwarningflag
        eml_warning('Coder:MATLAB:gammaincinv_FailedToConverge', ...
            eml_flt2str(pIn),eml_flt2str(a));
    end
end
rval = z;

%--------------------------------------------------------------------------

function z = PHIinv(p)
% Based on ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
q = p - 0.5;
% Coefficients for the central region
A0 = 3.3871328727963666080e+0;
A1 = 1.3314166789178437745e+2;
A2 = 1.9715909503065514427e+3;
A3 = 1.3731693765509461125e+4;
A4 = 4.5921953931549871457e+4;
A5 = 6.7265770927008700853e+4;
A6 = 3.3430575583588128105e+4;
A7 = 2.5090809287301226727e+3;
B1 = 4.2313330701600911252e+1;
B2 = 6.8718700749205790830e+2;
B3 = 5.3941960214247511077e+3;
B4 = 2.1213794301586595867e+4;
B5 = 3.9307895800092710610e+4;
B6 = 2.8729085735721942674e+4;
B7 = 5.2264952788528545610e+3;
% Coefficients for the near tails
C0 = 1.42343711074968357734e+0;
C1 = 4.63033784615654529590e+0;
C2 = 5.76949722146069140550e+0;
C3 = 3.64784832476320460504e+0;
C4 = 1.27045825245236838258e+0;
C5 = 2.41780725177450611770e-1;
C6 = 2.27238449892691845833e-2;
C7 = 7.74545014278341407640e-4;
D1 = 2.05319162663775882187e+0;
D2 = 1.67638483018380384940e+0;
D3 = 6.89767334985100004550e-1;
D4 = 1.48103976427480074590e-1;
D5 = 1.51986665636164571966e-2;
D6 = 5.47593808499534494600e-4;
D7 = 1.05075007164441684324e-9;
% Coefficients for the extreme tails
E0 = 6.65790464350110377720e+0;
E1 = 5.46378491116411436990e+0;
E2 = 1.78482653991729133580e+0;
E3 = 2.96560571828504891230e-1;
E4 = 2.65321895265761230930e-2;
E5 = 1.24266094738807843860e-3;
E6 = 2.71155556874348757815e-5;
E7 = 2.01033439929228813265e-7;
F1 = 5.99832206555887937690e-1;
F2 = 1.36929880922735805310e-1;
F3 = 1.48753612908506148525e-2;
F4 = 7.86869131145613259100e-4;
F5 = 1.84631831751005468180e-5;
F6 = 1.42151175831644588870e-7;
F7 = 2.04426310338993978564e-15;
if p <= 0
    if p == 0
        z = -eml_guarded_inf;
    else
        z = eml_guarded_nan;
    end
    return
elseif p >= 1
    if p == 1
        z = eml_guarded_inf;
    else
        z = eml_guarded_nan;
    end
    return
end
% Rational approximation for central region
if abs(q) <= 0.425  % .075 <= p <= .925
    r = 0.180625 - q*q;
    z = q*(((((((A7*r + A6)*r + A5)*r + A4)*r + A3)*r + A2)*r + A1)*r + A0)/...
        (((((((B7*r + B6)*r + B5)*r + B4)*r + B3)*r + B2)*r + B1)*r + 1);
else
    if q < 0
        r = sqrt(-log(p));
    else
        r = sqrt(-log(1-p));
    end
    % Rational approximation for near tails
    if r <= 5  % exp(-25) < p < .075 or .925 < p < 1-exp(-25)
        r = r - 1.6;
        z = (((((((C7*r + C6)*r + C5)*r + C4)*r + C3)*r + C2)*r + C1)*r + C0)/...
            (((((((D7*r + D6)*r + D5)*r + D4)*r + D3)*r + D2)*r + D1)*r + 1);
        
        % Rational approximation for extreme tails
    else % p < exp(-25) or 1-exp(-25) < p
        r = r - 5;
        z = (((((((E7*r + E6)*r + E5)*r + E4)*r + E3)*r + E2)*r + E1)*r + E0)/...
            (((((((F7*r + F6)*r + F5)*r + F4)*r + F3)*r + F2)*r + F1)*r + 1);
    end
    if q < 0
        z = -z;
    end
end

%--------------------------------------------------------------------------
