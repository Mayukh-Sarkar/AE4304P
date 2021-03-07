%% Analytical
dw = diff(w);
dw(length(dw)+1) = 0;
var__beta_c = sum(Sxx(:,1)'.*dw)/pi
var__phi_c = sum(Sxx(:,2)'.*dw)/pi
var__pbv_c = sum(Sxx(:,3)'.*dw)/pi
var__rbv_c = sum(Sxx(:,4)'.*dw)/pi
var__ay_c = sum(Sxx(:,5)'.*dw)/pi

ana = [var__beta_c var__phi_c var__pbv_c var__rbv_c var__ay_c];

var__beta_r = sum(Sxx_r(:,1)'.*dw)/pi
var__rbv_r = sum(Sxx_r(:,2)'.*dw)/pi

ana_r = [var__beta_r var__rbv_r];
%var__ay_r = sum(Sxx_r(:,3)'.*dw)/pi
%% var.m
vbeta_c  = var(beta_c)
vphi_c  = var(yt1(:,2))
vpbv_c  = var(yt1(:,3))
vrbv_c  = var(yt1(:,4))
vay_c  = var(a_y);
var_f = [vbeta_c vphi_c vpbv_c vrbv_c vay_c ];
%Reduced model
vbeta_r  = var(ytr1(:,1))
vrbv_r  = var(ytr1(:,2))
%vay_r  = var(a_y_r)
var_r = [vbeta_r vrbv_r ]
%% lyapunov


W = eye(2,2)

%Full Model
Bl=B(:,4:5);
L1   = lyap(A2,Bl*W*Bl');
L   = L1(1:4,1:4);
var_L = diag(L)';

var_ay = C1(5,:)*L1*C1(5,:)';
lya = [var_L(1) var_L(2) var_L(3) var_L(4) var_ay ];

Full_var = [ana' var_f' lya']
%Reduced model
Bl1=Br(:,4:5);
Lr   = lyap(Ar,Bl1*W*Bl1');
Lr   = Lr(1:2,1:2);
var_L_r = diag(Lr)';

lya1 = [var_L_r(1) var_L_r(2)];

red_var = [ana_r' var_r' lya1']




