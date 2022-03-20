close all; clear;
tf = 10;dt=0.001;

t=0:dt:tf;
tspan = flip(t);


A=1.5;B=0.5;H=1;Q=1;R=4;
params.A= A;
params.B= B;
params.H= H;
params.R= R;
params.Q= Q;
state_initial =H;
[tspan,PP] = ode45(@(t,yy) dynamics(t,yy,params), tspan, state_initial);
PP = flip(PP);
len = length(t);

K= zeros(1,len-1);
for i=1:length(PP)-1   
    
    K(i) = R\B'*PP(i);
    
end
X = zeros(1,len);
U = zeros(1,len-1);
A= 1+A*dt;
B= B*dt;
X(1)=10;

for i= 1:len-1
    U(i) = -K(i)*X(i);
    X(i+1) = A*X(i)+B*U(i);
end


plot(t,PP,'r', 'LineWidth',1.5); 
grid on
% ylim([0.4 max(max(PP),max(K))])
hold on;
% figure
plot(t(1:end-1),K,'b', 'LineWidth',1.5); grid on;
title(['Plot of P(t) and K(t) using ricaati eqn numerical'])
legend('PP(k)', 'K(k)')
saveas(gcf,[pwd '/imgQ2/gainplotQ2'],'epsc')


figure
plot(t,X(1:len),'r','Linewidth',2);hold on; grid
plot(t(1:end-1),U,'b','Linewidth',2),legend('X(k)','U(k)')
title(['Plot of X(t) and U(t) using ricaati eqn numerical'])
saveas(gcf,[pwd '/imgQ2/XandUplotQ2'],'epsc')


function [dxstate] = dynamics(t,state,params)
    A=params.A;
    B=params.B;
    H=params.H;
    R=params.R;
    Q=params.Q;
    P= state;
    dxstate = -1-3*P+1/16*(P^2);

% dxstate = -(Q+P*A+A'*P-P*B*inv(R)*B'*P);
end