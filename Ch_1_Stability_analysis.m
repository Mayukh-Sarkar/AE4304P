% Filename:  cit2a.m
%
% Calculation of state matrix and input matrix for calculation
% of asymmetric aircraft response to atmospheric turbulence.
% The system model is in the form
%
%       .
%       x = Ax + Bu
%       -    -    -
% with x = [beta phi pb/2V rb/2V ug_ u_g* alpha_g alpha_g* beta_g betag*]'
%
% and
%
% u = [delta_a delta_r w1 w2 w3]'.
%
%
% The turbulence filters are derived using the approximated
% effective one-dimensional power spectral densities for u_g, alpha_g and
% beta_g.
% 
% Deig

% AIRCRAFT- AND FLIGHT CONDITION 'cruise'.
V   = 121.3; %check
W = 53361; %check
rho = 0.4587;%check
S   = 24.2;%check
b   = 13.36;%check
mub = 32;%check
KX2 = 0.013; %check
KZ2 = 0.037; %check
KXZ = 0.002; %check
CL  = W/(0.5*rho*V^2*S); %check

% TURBULENCE PARAMETERS APPROXIMATED POWER SPECTRAL DENSITIES
Lg        = 150; %check
B         = b/(2*Lg);
sigma     = 2; %check
sigmaug_V = sigma/V;
sigmavg   = 1; %check
sigmabg   = sigmavg/V;
sigmaag   = sigma/V;

Iug0 = 0.0249*sigmaug_V^2;
Iag0 = 0.0182*sigmaag^2;
tau1 = 0.0991;     tau2 = 0.5545;     tau3 = 0.4159;
tau4 = 0.0600;     tau5 = 0.3294;     tau6 = 0.2243;

% AIRCRAFT ASYMMETRIC AERODYNAMIC DERIVATIVES 
CYb  =-1.3250;     Clb  =-0.1070;     Cnb  = 0.1835;
CYp  =-0.1320;     Clp  =-0.3684;     Cnp  =0.0035;
CYr  = 0.4300;     Clr  = 0.1750;     Cnr  =-0.1930;
CYda = 0.0000;     Clda =-0.2292;     Cnda = 0.0071;
CYdr = 0.3037;     Cldr = 0.0446;     Cndr =-0.1261;
 
                   Clpw = 0.8*Clp;    Cnpw = 0.9*Cnp;
                   Clrw = 0.7*Clr;    Cnrw = 0.2*Cnr;
CYfb = 0;
Clfb = 0;
Cnfb = 0;

%CYfbg = CYfb+0.5*CYr;
%Clfbg = Clfb+0.5*Clr;
%Cnfbg = Cnfb+0.5*Cnr;

% CALCULATION OF AIRCRAFT ASYMMETRIC STABILITY DERIVATIVES
yb   = (V/b)*CYb/(2*mub);
yphi = (V/b)*CL/(2*mub);
yp   = (V/b)*CYp/(2*mub);
yr   = (V/b)*(CYr-4*mub)/(2*mub);
ybg  = yb;
ydr  = (V/b)*CYdr/(2*mub);
den  = b*4*mub*(KX2*KZ2-KXZ^2)/V;
lb   = (Clb*KZ2+Cnb*KXZ)/den;
lp   = (Clp*KZ2+Cnp*KXZ)/den;
lr   = (Clr*KZ2+Cnr*KXZ)/den;
lda  = (Clda*KZ2+Cnda*KXZ)/den;
ldr  = (Cldr*KZ2+Cndr*KXZ)/den;
lug  = (-Clrw*KZ2-Cnrw*KXZ)/den;
lbg  = lb;
lag  = (Clpw*KZ2+Cnpw*KXZ)/den;
nb   = (Clb*KXZ+Cnb*KX2)/den;
np   = (Clp*KXZ+Cnp*KX2)/den;
nr   = (Clr*KXZ+Cnr*KX2)/den;
nda  = (Clda*KXZ+Cnda*KX2)/den;
ndr  = (Cldr*KXZ+Cndr*KX2)/den;
nug  = (-Clrw*KXZ-Cnrw*KX2)/den;
nbg  = nb;
nag  = (Clpw*KXZ+Cnpw*KX2)/den;
aug1 =-(V/Lg)^2*(1/(tau1*tau2));
aug2 =-(tau1+tau2)*(V/Lg)/(tau1*tau2);
aag1 =-(V/Lg)^2*(1/(tau4*tau5));
aag2 =-(tau4+tau5)*(V/Lg)/(tau4*tau5);
abg1 =-(V/Lg)^2;
abg2 =-2*(V/Lg);
bug1 = tau3*sqrt(Iug0*V/Lg)/(tau1*tau2);
bug2 = (1-tau3*(tau1+tau2)/(tau1*tau2))*sqrt(Iug0*(V/Lg)^3)/(tau1*tau2);
bag1 = tau6*sqrt(Iag0*V/Lg)/(tau4*tau5);
bag2 = (1-tau6*(tau4+tau5)/(tau4*tau5))*sqrt(Iag0*(V/Lg)^3)/(tau4*tau5);
bbg1 = sigmabg*sqrt(3*V/Lg);
bbg2 = (1-2*sqrt(3))*sigmabg*sqrt((V/Lg)^3);


A = [yb yphi yp    yr 0    0    0    0    ybg  0;
     0  0    2*V/b 0  0    0    0    0    0    0;
     lb 0    lp    lr lug  0    lag  0    lbg  0;
     nb 0    np    nr nug  0    nag  0    nbg  0;
     0  0    0     0  0    1    0    0    0    0;
     0  0    0     0  aug1 aug2 0    0    0    0;
     0  0    0     0  0    0    0    1    0    0;
     0  0    0     0  0    0    aag1 aag2 0    0;
     0  0    0     0  0    0    0    0    0    1;
     0  0    0     0  0    0    0    0    abg1 abg2];

B = [0   ydr 0    0    0;
     0   0   0    0    0;
     lda ldr 0    0    0;
     nda ndr 0    0    0;
     0   0   bug1 0    0;
     0   0   bug2 0    0;
     0   0   0    bag1 0;
     0   0   0    bag2 0;
     0   0   0    0    bbg1;
     0   0   0    0    bbg2];
 C = eye(10,10);
 D = zeros(10,5);

eig_A = eig(A)
sys_uncontrolled = ss(A,B,C,D);
figure(1)
pzmap(sys_uncontrolled)
grid on

Ar = [yb     yr 0    0    0    0    ybg  0;
     nb     nr nug  0    nag  0    nbg  0;
     0       0  0    1    0    0    0    0;
     0       0  aug1 aug2 0    0    0    0;
     0       0  0    0    0    1    0    0;
     0       0  0    0    aag1 aag2 0    0;
     0       0  0    0    0    0    0    1;
     0       0  0    0    0    0    abg1 abg2];

Br = [0   ydr 0    0    0;
     nda ndr 0    0    0;
     0   0   bug1 0    0;
     0   0   bug2 0    0;
     0   0   0    bag1 0;
     0   0   0    bag2 0;
     0   0   0    0    bbg1;
     0   0   0    0    bbg2];
Cr = eye (8,8);
Dr = zeros(8,5); 
eig_Ar = eig(Ar)
sys_reduced = ss(Ar,Br,Cr,Dr);
figure(2)
pzmap(sys_reduced)
grid on
Kphi = -0.2;
K    = [0 Kphi 0 0  0 0  0 0  0 0];
A2   = A-B(:,1)*K;

eig_A2 = eig(A2)
sys_controlled = ss(A2,B,C,D);
figure(3)
pzmap(sys_controlled)
grid on
