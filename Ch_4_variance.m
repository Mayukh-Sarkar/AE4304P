%% Analytical
dw = diff(w);
dw(length(dw)+1) = 0;
var__beta_c = sum(Sxx(:,1)'.*dw)/pi
var__phi_c = sum(Sxx(:,2)'.*dw)/pi
var__pbv_c = sum(Sxx(:,3)'.*dw)/pi
var__rbv_c = sum(Sxx(:,4)'.*dw)/pi
var__ay_c = sum(Sxx(:,5)'.*dw)/pi

var__beta_r = sum(Sxx_r(:,1)'.*dw)/pi
var__rbv_r = sum(Sxx_r(:,2)'.*dw)/pi
var__ay_r = sum(Sxx_r(:,3)'.*dw)/pi
%% var.m
vbeta_c  = var(yt(:,1))
vphi_c  = var(yt(:,2))
vpbv_c  = var(yt(:,3))
vrbv_c  = var(yt(:,4))
vay_c  = var(a_y)
%Reduced model
vbeta_r  = var(ytr(:,1))
vrbv_r  = var(ytr(:,2))
vay_r  = var(a_y_r)
    
%% lyapunov


W = eye(2,2);
%Full Model
B=B(:,4:5);
L   = lyap(A2,B*W*B');
L   = L(1:4,1:4);
var_L = diag(L)'

%Reduced model
Br=Br(:,4:5);
Lr   = lyap(Ar,Br*W*Br');
Lr   = Lr(1:2,1:2);
var_L_r = diag(Lr)'


