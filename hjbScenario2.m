close all; clear all;
tf = 5;dt=0.001;

% A= [0 1 ; 0 1]; B=[0; 1];H=[0 0; 0 4]; Q=[ 50 0; 0 0.01]; R=3;
% Q=[1 0; 0 0];

A= [-1 2 ; 0 1]; B=[-1; 1];
H=[2 0; 0 1]; 
R=1;
Q=[2 1; 1 2];
% Q=[ 50 0; 0 0.01];

params.A= A;
params.B= B;
params.H= H;
params.R= R;
params.Q= Q;

% Q=[ 50 0; 0 0.01];
t=0:dt:tf;
tspan = flip(t);
state_initial =H;

[tspan,PP] = ode45(@(t,yy) dynamics(t,yy,params), tspan, state_initial);
PP = flip(PP);
% K=0.5*yy;
% % A=1;B=1;Q=0;H=0.5;R=0.5;
% A=3; B=11; Q=7;R=400;
len = length(t);

for i=1:length(PP)-1   
    P = reshape(PP(i,:),2,2);
    K(:,i) = inv(R)*B'*P;
    
end
X = zeros(2,len);
Xs = zeros(2,len);
U = zeros(1,len-1);
Uss = zeros(1,len-1);

% A= 1+A*dt;
% B= B*dt;
X(:,1)=[-6 8];
Xs(:,1) = [-6 8];
% xcurr = [1 1];
for i= 1:len-1
    U(i) = -K(:,i)'*X(:,i);
    Uss(i) = -K(:,1)'*Xs(:,i);
    %     X(:,i+1) = A*X(:,i)+B*U(i);
    X(:,i+1)= (A-B*K(:,i)')*X(:,i)*dt+ X(:,i);
    Xs(:,i+1)= (A-B*K(:,1)')*Xs(:,i)*dt+ Xs(:,i);

end

% subplot(2,1,1);
figure
plot(t,PP(:,1),'r', 'LineWidth',2); 
hold on;grid on
plot(t,PP(:,4),'g', 'LineWidth',2); 
title(['R=', num2str(R)])
xlabel('Time (sec)');ylabel('P')
legend('P1(k)', 'P2(k)')
saveas(gcf,[pwd '/imgQ12/PplotQ12'],'epsc')


figure
plot(t(1:end-1),K(1,:), 'r','LineWidth',2); hold on; grid;
plot(t(1:end-1),K(1,1)*ones(1,len-1), 'b','LineWidth',1.5);
xlabel('Time (sec)');ylabel('Gains')
legend('K1(t)', 'K1(ss)','Location','southwest')
title(['Plot of feedback gain K1 with R=', num2str(R)])
c= num2str(K(1,1));
text(t(10), K(1,1)+0.09, c, 'Fontsize', 10);
saveas(gcf,[pwd '/imgQ12/GainplotQ12'],'epsc')


figure
plot(t(1:end-1),K(2,:), 'r','LineWidth',2); hold on; grid;
plot(t(1:end-1),K(2,1)*ones(1,len-1), 'b','LineWidth',2);
title(['Plot of feedback gain K2 with R=', num2str(R)])
legend('K2(t)', 'K2ss)','Location','southwest')
c= num2str(K(2,1));
text(t(10), K(2,1)+0.09, c, 'Fontsize', 10);
xlabel('Time (sec)');ylabel('Gains')
saveas(gcf,[pwd '/imgQ12/Gain2plotQ12'],'epsc')


figure
% subplot(2,1,1);
plot(t,X,'Linewidth',2); grid; 
title(['Plot of states(X1 and X2) with R=', num2str(R)])
legend('X_1(t)', 'X_2(t)')
xlabel('Time (sec)');ylabel('States')
saveas(gcf,[pwd '/imgQ12/statesplotQ12'],'epsc')



figure
plot(t,Xs,'Linewidth',2); grid; 
legend('X_1(ss)', 'X_2(ss)')
xlabel('Time (sec)');ylabel('States')
% title(['Plot of states(X1 and X2) with R=', num2str(R)], 'and K_ss')
saveas(gcf,[pwd '/imgQ12/XsplotQ12'],'epsc2')


figure
plot(t(1:end-1),U,'r','Linewidth',2); grid; 
xlabel('Time (sec)');ylabel('input');
% figure
hold on;
plot(t(1:end-1),Uss,'b','Linewidth',1); 
legend('U(t)',"U(ss)")
saveas(gcf,[pwd '/imgQ12/inputplotQ12'],'epsc')


function [dxstate] = dynamics(t,state,params)
    A=params.A;
    B=params.B;
    H=params.H;
    R=params.R;
    Q=params.Q;

    P= reshape(state(:,:),2,2);
%     K = R\B'*P;
    dxstate = -(Q+P*A+A'*P-P*B*inv(R)*B'*P);
    dxstate = reshape(dxstate,2*2,1);

end


