%% Analytical
dw = diff(w);
dw(length(dw)+1) = 0;
var__beta_c = sum(Sxx(:,1)'.*dw)/pi;
var__phi_c = sum(Sxx(:,2)'.*dw)/pi;
var__pbv_c = sum(Sxx(:,3)'.*dw)/pi;
var__rbv_c = sum(Sxx(:,4)'.*dw)/pi;
var__ay_c = sum(Sxx(:,5)'.*dw)/pi;

var__beta_r = sum(Sxx_r(:,1)'.*dw)/pi;
var__rbv_r = sum(Sxx_r(:,2)'.*dw)/pi;
%var__ay_r = sum(Sxx_r(:,3)'.*dw)/pi
%% var.m
vbeta_c  = var(yt1(:,1));
vphi_c  = var(yt1(:,2));
vpbv_c  = var(yt1(:,3));
vrbv_c  = var(yt1(:,4));


vay_c  = var(a_y);
%Reduced model
vbeta_r  = var(ytr1(:,1))
vrbv_r  = var(ytr1(:,2))
%vay_r  = var(a_y_r)
    
%% lyapunov


W = 1

%Full Model
Bl=B(:,4:5);
L1   = lyap(A2,Bl*W*Bl');
L   = L1(1:4,1:4);
var_L = diag(L1)'

var_ay = C(5,:)*L1*C(5,:)'

%Reduced model
Bl1=Br(:,4:5);
Lr   = lyap(Ar,Bl1*W*Bl1');
Lr   = Lr(1:2,1:2);
var_L_r = diag(Lr)';






