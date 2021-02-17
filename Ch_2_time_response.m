dt = 0.05; T  = 60; t = [0:dt:T]; N = length(t);
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
C = [1 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0 0];
   
D = [0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0];

% RESPONSE to u_g
y1 = lsim(A2,B,C,D,u1,t);
% RESPONSE to v_g
y2 = lsim(A2,B,C,D,u2,t);
% RESPONSE to w_g
y3 = lsim(A2,B,C,D,u3,t);
% RESPONSE to all together (linear system!)
yt = y1+y2+y3;
a_y = V*(yt(:,1)+yt(:,2));
y1r = lsim(Ar,Br,Cr,Dr,u1,t);
% RESPONSE to v_g
y2r = lsim(Ar,Br,Cr,Dr,u2,t);
% RESPONSE to w_g
y3r = lsim(Ar,Br,Cr,Dr,u3,t);
ytr = y1r+y2r+y3r;
a_y_r = V*(ytr(:,1)+ytr(:,2));
% PLOT RESULTS
figure(1)
subplot(2,1,1);
plot(t,yt(:,1),'g','LineWidth',2)
xlabel('Time (s) ')
ylabel('\beta_c (rad/s)')
grid on
subplot(2,1,2); 
plot(t,ytr(:,1),'b','LineWidth',2)
xlabel('Time (s)')
ylabel('\beta_r (rad/s')
grid on

figure(2)
subplot(2,1,1);
plot(t,yt(:,4),'g','LineWidth',2)
xlabel('Time (s) ')
ylabel('rb/2V (rad/s)')
grid on
subplot(2,1,2); 
plot(t,ytr(:,2),'b','LineWidth',2)
xlabel('Time (s)')
ylabel('rb/2V (rad/s)')
grid on
figure(3)
subplot(2,1,1);
plot(t,yt(:,2),'g','LineWidth',2)
xlabel('Time (s) ')
ylabel('\phi (rad/s)')
grid on
subplot(2,1,2); 
plot(t,yt(:,3),'b','LineWidth',2)
xlabel('Time (s)')
ylabel('pb/2V (rad/s)')
grid on
figure(4)
subplot(2,1,1);
plot(t,a_y,'g','LineWidth',2)
xlabel('Time (s) ')
ylabel('a_y_c (m/s)')
grid on
subplot(2,1,2); 
plot(t,a_y_r,'b','LineWidth',2)
xlabel('Time (s)')
ylabel('a_y_r (m/s)')
grid on
figure(5)
subplot(2,1,1);
plot(t,yt,'g','LineWidth',2)
xlabel('Time (s) ')
ylabel('y_c')
grid on
subplot(2,1,2); 
plot(t,ytr,'b','LineWidth',2)
xlabel('Time (s)')
ylabel('y_r ')
grid on