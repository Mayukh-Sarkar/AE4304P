dt = 0.005; T  = 600; 
t = [0:dt:T]; N = length(t);
nn = zeros(1,N);

% TURBULENCE INPUTS
u_g = randn(1,N)/sqrt(dt);    % sqrt(dt) because of lsim characteristics
v_g = randn(1,N)/sqrt(dt);
w_g = randn(1,N)/sqrt(dt);

% INPUT VECTORS
u1 = [nn' nn' u_g' nn'  nn'];     
u2 = [nn' nn' nn'  v_g' nn'];
u3 = [nn' nn' nn'  nn'  w_g'];

% DEFINE OUTPUT MATRICES
C1 = [1 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0 0;
     V*yb V*yphi (b/(2*Lg))*yp V*(yr+(2*V/b)) 0 0 0 0 V*ybg 0];
 
   
D = [0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0];

Cr1= [1 0 0 0 0 0 0 0;
      0 1 0 0 0 0 0 0;
      V*yb  V*(yr+(2*V/b)) 0 0 0 0 V*ybg 0];
Dr1 = [0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0];
% RESPONSE to u_g
y1 = lsim(A2,B,C1,D,u1,t);
% RESPONSE to v_g
y2 = lsim(A2,B,C1,D,u2,t);
% RESPONSE to w_g
y3 = lsim(A2,B,C1,D,u3,t);
% RESPONSE to all together (linear system!)
yt = y1+y2+y3;
%xt = x1+x2+x3;
a_y =  yt(:,5);
 
y1r = lsim(Ar,Br,Cr1,Dr1,u1,t);
% RESPONSE to v_g
y2r = lsim(Ar,Br,Cr1,Dr1,u2,t);
% RESPONSE to w_g
y3r = lsim(Ar,Br,Cr1,Dr1,u3,t);
ytr = y1r+y2r+y3r;
a_y_r = ytr(:,3);

% PLOT RESULTS
figure(4)
subplot(2,1,1);
plot(t,yt(:,1),'g')
xlabel('Time (s) ')
ylabel('\beta_c (rad/s)')
grid on
subplot(2,1,2); 
plot(t,ytr(:,1),'b')
xlabel('Time (s)')
ylabel('\beta_r (rad/s')
grid on

figure(5)
subplot(2,1,1);
plot(t,yt(:,4),'g')
xlabel('Time (s) ')
ylabel('rb/2V_c (rad/s)')
grid on
subplot(2,1,2); 
plot(t,ytr(:,2),'b')
xlabel('Time (s)')
ylabel('rb/2V_r (rad/s)')
grid on
figure(6)
subplot(2,1,1);
plot(t,yt(:,2),'g')
xlabel('Time (s) ')
ylabel('\phi (rad/s)')
grid on
subplot(2,1,2); 
plot(t,yt(:,3),'b')
xlabel('Time (s)')
ylabel('pb/2V (rad/s)')
grid on
figure(7)
subplot(2,1,1)
plot(t,a_y,'g')
xlabel('Time (s) ')
ylabel('a_y_c (m/s)')
grid on
subplot(2,1,2); 
plot(t,a_y_r,'b','LineWidth',1)
xlabel('Time (s)')
ylabel('a_y_r (m/s)')
grid on
figure(8)
subplot(2,1,1);
plot(t,yt,'g')
xlabel('Time (s) ')
ylabel('y_c')
grid on
subplot(2,1,2); 
plot(t,ytr,'b')
xlabel('Time (s)')
ylabel('y_r ')
grid on

