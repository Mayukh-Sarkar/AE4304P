
w = logspace(-2,2,300);
%% analystical
% controlled model
temp = bode(A2,B,C(1,:),D(1,:),4,w)+ bode(A2,B,C(1,:),D(1,:),5,w); Sbeta_c  = temp.*temp;
temp = bode(A2,B,C(2,:),D(2,:),4,w)+ bode(A2,B,C(2,:),D(2,:),5,w); Sphi_c   = temp.*temp;
temp = bode(A2,B,C(3,:),D(3,:),4,w)+ bode(A2,B,C(3,:),D(3,:),5,w); Sbp_c    = temp.*temp;
temp = bode(A2,B,C(4,:),D(4,:),4,w)+ bode(A2,B,C(4,:),D(4,:),5,w); Srb_c    = temp.*temp;
temp = bode(A2,B,C1(5,:),D(5,:),4,w)+ bode(A2,B,C1(5,:),D(5,:),5,w); Say_c = temp.*temp;

Sxx  = [Sbeta_c Sphi_c Sbp_c Srb_c Say_c];

%reduced model
temp = bode(Ar,Br,Cr1(1,:),Dr(1,:),4,w) + bode(Ar,Br,Cr1(1,:),Dr(1,:),5,w); Sbeta_r  = temp.*temp;
temp = bode(Ar,Br,Cr1(2,:),Dr(2,:),4,w)+ bode(Ar,Br,Cr1(2,:),Dr(2,:),5,w); Srb_r    = temp.*temp;
temp = bode(Ar,Br,Cr1(3,:),D(3,:),4,w) + bode(Ar,Br,Cr1(3,:),D(3,:),5,w);Say_r = temp.*temp;

Sxx_r = [Sbeta_r Srb_r Say_r];
u = [nn' nn' nn'  w_g'  v_g'];

yt1 = lsim(A2,B,C1,D,u,t);
ytr1= lsim(Ar,Br,Cr,Dr,u,t);


beta_c = yt1(:,1);
phi_c = yt1(:,2);
pbV_c = yt1(:,3);
rbV_c= yt1(:,4);


beta_r = ytr1(:,1);
rbV_r= ytr1(:,2);
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


% DEFINE FREQUENCY VECT
fs = 1/dt;     % sample frequency
omega = 2*pi*fs*(0:(N/2)-1)/N;
%% using pwelch

pBETA_c  = pwelch(beta_c,omega,[],N,fs);
pBETA_c = pBETA_c/2;
pPhi_c   = pwelch(phi_c,omega,[],N,fs);
pPhi_c = pPhi_c/2;
pPbv_c   = pwelch(pbV_c,omega,[],N,fs);
pPbv_c = pPbv_c/2;
pRbv_c     = pwelch(rbV_c,omega,[],N,fs);
pRbv_c = pRbv_c/2;
pay_c = pwelch(a_y,omega,[],N,fs);
pay_c = pay_c/2;
pBETA_r  = pwelch(beta_r,omega,[],N,fs);
pBETA_r = pBETA_r/2;
pRbv_r     = pwelch(rbV_r,omega,[],N,fs);
pRbv_r = pRbv_r/2;
pay_r = pwelch(a_y_r,omega,[],N,fs);
pay_r = pay_r/2;

 %PSD ESTIMATE
figure(9)
subplot(2,1,1); loglog(omega,beta_psd_c (1:N/2),'-',omega,pBETA_c(2:N/2+1),w,Sxx(:,1),'k') 
axis(10.^[-1 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('S_\beta_\beta_c[rad^2/rad/s]');
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
subplot(2,1,2); loglog(omega,beta_psd_r(1:N/2),'-',omega,pBETA_r(2:N/2+1),w,Sxx_r(:,1),'k')
axis(10.^[-1 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('S_\beta_\beta_r[rad^2/rad/s] ')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on

figure(10)
subplot(2,1,1); loglog(omega,phi_psd_c(1:N/2),'-',omega,pPhi_c(2:N/2+1),w,Sxx(:,2),'k');
axis(10.^[-1 2 -12 0]); xlabel('omega [rad/s]'); ylabel('S_\phi_\phi_c [rad^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
subplot(2,1,2); loglog(omega,pbv_psd_c(1:N/2),'-',omega,pPbv_c(2:N/2+1),w,Sxx(:,3),'k');
axis(10.^[-1 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('S_p_p [rad^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on

figure(11)
subplot(2,1,1)
loglog(omega,rbv_psd_c(1:N/2),'-',omega,pRbv_c(2:N/2+1),w,Sxx(:,4),'k');
axis(10.^[-1 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('S_r_r_c [rad^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
subplot(2,1,2); loglog(omega,rbv_psd_r(1:N/2),'-',omega,pRbv_r(2:N/2+1),w,Sxx_r(:,2),'k');
axis(10.^[-1 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('S_r_r_r [rad^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on

figure(12)
subplot(2,1,1)
loglog(omega,ay_psd_c(1:N/2),'-',omega,pay_c(2:N/2+1),w,Sxx(:,5),'k');
axis(10.^[-1 2 -10 5]); xlabel('omega [rad/s]'); ylabel('S_a_a [(m/s^2)^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
subplot(2,1,2); 
loglog(omega,ay_psd_r(1:N/2),'-',omega,pay_r(2:N/2+1),w,Sxx_r(:,3),'k');
axis(10.^[-1 2 -10 5]); xlabel('omega [rad/s]'); ylabel('S_a_a [(m/s^2)^2/rad/s]')
legend('Experimental Periodogram','Smoothed Periodogram','Analytical PSD')
grid on
figure(13) % for beta comptete model
bode(A2,B,C(1,:),D(1,:),4,w)
hold on
bode(A2,B,C(1,:),D(1,:),5,w)
legend('lateral turbulence','vertical turbulence')
% for rb/2v full model
figure(14)
bode(A2,B,C(4,:),D(4,:),4,w)
hold on
bode(A2,B,C(4,:),D(4,:),5,w)
legend('lateral turbulence','vertical turbulence')