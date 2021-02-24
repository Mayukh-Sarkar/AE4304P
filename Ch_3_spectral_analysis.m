
w = logspace(-2,2,300);
%% analystical
% controlled model
temp = bode(A2,B,C(1,:),D(1,:),5,w); Sbeta_c  = temp.*temp;
temp = bode(A2,B,C(2,:),D(2,:),5,w); Sphi_c   = temp.*temp;
temp = bode(A2,B,C(3,:),D(3,:),5,w); Spp_c    = temp.*temp;
temp = bode(A2,B,C(4,:),D(4,:),5,w); Srb_c    = temp.*temp;
temp = bode(A2,B,V*(A2(1,:)+[0 0 0 2*V/b 0    0    0    0    0  0]),V*B(1,:),4,w); Say_c = temp.*temp;

Sxx  = [Sbeta_c Sphi_c Spp_c Srb_c Say_c];

%reduced model
temp = bode(Ar,Br,Cr(1,:),Dr(1,:),5,w); Sbeta_r  = temp.*temp;
temp = temp + bode(Ar,Br,Cr(2,:),Dr(2,:),5,w); Srb_r    = temp.*temp;
temp = bode(Ar,Br,V*(Ar(1,:)+[0 2*V/b 0    0    0    0    0  0]),V*Br(1,:),5,w);Say_r = temp.*temp;

Sxx_r = [Sbeta_r,Srb_r,Say_r];

v_g = randn(N,1)/sqrt(dt);    % sqrt(dt) because of lsim
w_g = randn(N,1)/sqrt(dt);
nn = zeros(N,1);
u  = [nn nn nn w_g v_g];

Y_c = lsim(A2,B,C,D,u,t);
beta_c = Y_c(:,1);
phi_c = Y_c(:,2);
pbV_c = Y_c(:,3);
rbV_c= Y_c(:,4);

Y_r = lsim(Ar,Br,Cr,Dr,u,t);
beta_r = Y_r(:,1);
rbV_r= Y_r(:,2);
% COMPUTE PERIODOGRAM AND ESTIMATE PSD
%% using fft
BETA_c  = dt*fft(beta_c);
Phi_c   = dt*fft(phi_c);
Pbv_c     = dt*fft(pbV_c);
Rbv_c     = dt*fft(rbV_c);
ay_c = dt*fft(a_y);

%  PSD ESTIMATE
beta_psd_c  = (1/T)*( BETA_c.*conj(BETA_c));
phi_psd_c  = (1/T)*(  Phi_c.*conj(Phi_c));
pbv_psd_c     = (1/T)*(    Pbv_c.*conj(Pbv_c));
rbv_psd_c= (1/T)*(    Rbv_c.*conj(Rbv_c));
ay_psd_c = (1/T)*(    ay_c.*conj(ay_c));

%reduced model
BETA_r  = dt*fft(beta_r);
Rbv_r     = dt*fft(rbV_r);
ay_r = dt*fft(a_y_r);

beta_psd_r  = (1/T)*( BETA_r.*conj(BETA_r));
rbv_psd_r= (1/T)*(    Rbv_r.*conj(Rbv_r));
ay_psd_r = (1/T)*(    ay_r.*conj(ay_r));


Periodograms = [beta_psd_c phi_psd_c pbv_psd_c rbv_psd_c ay_psd_c];
Periodograms_r = [beta_psd_r rbv_psd_r ay_psd_r];
% DEFINE FREQUENCY VECTOR
fs = 1/dt;     % sample frequency
omega = 2*pi*fs*(0:(N/2)-1)/N;
%% using pwelch

pBETA_c  = pwelch(beta_c);
pPhi_c   = pwelch(phi_c);
pPbv_c     = pwelch(pbV_c);
pRbv_c     = pwelch(rbV_c);
pay_c = pwelch(a_y);

%  PSD ESTIMATE
% pbeta_psd_c  = ( pBETA_c.*conj(pBETA_c));
% pphi_psd_c  = (  pPhi_c.*conj(pPhi_c));
% ppbv_psd_c     = (    pPbv_c.*conj(pPbv_c));
% prbv_psd_c= (    pRbv_c.*conj(pRbv_c));
% pay_psd_c = (    pay_c.*conj(pay_c));
figure(1)
subplot(2,1,1); loglog(w,Sxx(:,1),'--',omega,beta_psd_c (1:N/2),'-',omega(1:257),pBETA_c) 
axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('\beta PSD_c');
% subplot(2,1,2); loglog(w,Sxx_r(:,1),'--',omega,beta_psd_r(1:N/2),beta_r(1:N/2))
% axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('\beta PSD_c')

% figure(2)
% subplot(2,1,1); loglog(w,Sxx(:,2),'--',omega,phi_psd_c(1:N/2));
% axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('\phi PSD_c')
% subplot(2,1,2); loglog(w,Sxx(:,3),'--',omega,pbv_psd_c(1:N/2));
% axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('pb/V PSD_c')
% 
% figure(3)
% subplot(2,1,1); loglog(w,Sxx(:,4),'--',omega,rbv_psd_c(1:N/2));
% axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('rb/V PSD_c')
% subplot(2,1,2); loglog(w,Sxx_r(:,2),'--',omega,rbv_psd_r(1:N/2));
% axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('rb/V PSD_r')
% 
