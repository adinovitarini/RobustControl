clear all;clc
%% Nominal Plant Dynamics
s=tf('s');
N = {[1 -100],[10 10];[-10 -10],[1 -100]};
D = [1 0 100];
Pnom = tf(N,D);
Pnom = ss(Pnom,'min')
SystemPole = pole(Pnom);
SystemZero = tzero(Pnom);
%% Plant with uncertainty
for i = 1:2
    k(i) = ureal('k',1,'Range',[0.8,1.2]);
    theta(i) = ureal('theta',0.01,'Range',[0,.02]);
    f(i) = k(i)*((-.5*theta(i)*s)+1)/((.5*theta(i)*s)+1);
end
P = [f(1) 0;0 f(2)]*Pnom;
WDelta = (P/Pnom)-1;
r0 = 0.2; rinf1 = 2.3; tau = 0.021;
rinf2 = 100;
wM = tf([tau r0], [tau/rinf1 1]);
WM = eye(2)*wM;
WM = ss(WM);
wM2= tf([tau r0], [tau/rinf2 1]);
WM2 = eye(2)*wM2;
WM2 = ss(WM2) ;
%% Multiplicative Output Uncertainty 
time = [0:0.01:3];
step_ref = ones(1,length(time));
ref = [0.1*step_ref' 0*step_ref'];
ref2 = [0.1*sin(10*time)' 0.1*cos(10*time)'];
Ms = 2; A = 1e-2; wb = 11.5; 
Ms2 = 8;wb2 = 20;
wP = tf([1/Ms wb], [1 wb*A]);
WP = eye(2)*wP; WP = ss(WP) ;
wP2 = tf([1/Ms2 wb2], [1 wb2*A]);
WP2 = eye(2)*wP2; WP2 = ss(WP2) ;
systemnames = 'Pnom WP WM';
inputvar = '[w(2);u(2)]';
outputvar = '[WP;WM;-w-Pnom]';
input_to_Pnom= '[u]';
input_to_WP = '[w+Pnom]';
input_to_WM = ' [Pnom]';
G = sysic
nmeas = 2; ncon = 2;
[Khi,CLhi,ghi,hiinfo] = hinfsyn(G,nmeas,ncon);
Fhi=loopsens(Pnom,Khi);
%% Trade off gamma based on Gmax and Gmin
[Khi,CLHi,ghi,hiinfo] = hinfsyn(G,nmeas,ncon,'Gmax',1,'Gmin',0);
lft_norm = hinfnorm(lft(G,Khi))
[SV,w] = sigma(CLHi)
%% Plotting Figure
figure(1);clf
pzmap(tf(N,D));
title('Pole-Zero map of $P(s)$','Interpreter', 'LaTeX');
figure(2);clf
sigmaplot(Pnom);
text(13,20,'$$\bar{\sigma}P(s)$$','Interpreter','latex');
hold on;
text(10,-12,'$$\underline{\sigma}P(s)$$','Interpreter','latex')
title('sigmaplot of P(s)');
figure(3);bodemag(WDelta(1,1))
hold on;bodemag(WM(1,1),'r');
grid on
title('Weights for Multiplicative Uncertainty')
legend('$$Wm(s)\Delta(s)$$','$$Wm(s)=\frac{0.021s+0.2}{0.0091s+1}$$','Interpreter','latex')
figure(4);bodemag(wM,'b');hold on;bodemag(wM2,'r');
legend('r(\infty) =2.3','r(\infty) = 100')
figure(5)
lsim(Pnom(1,:),ref,time','r')
hold on
lsim(P(1,:),ref,time','b')
hold on
plot(time',ref(:,1),'-.r')
lsim(Pnom(2,:),ref,time','m')
hold on
lsim(P(2,:),ref,time','c')
hold on
plot(time',ref(:,2),'-.m')
hl = title('Time response with $\omega=[0.1step;0]^{T}$');
set(hl ,'Interpreter','latex');
legend('$P_{1}(s)$','$$\tilde{P}_{1}(s)$$','$\omega_{1}(s)$','$P_{2}(s)$','$$\tilde{P}_{2}(s)$$','$\omega_{2}(s)$','Interpreter', 'LaTeX');
grid on
figure(6)
lsim(Pnom(1,:),ref2,time','r')
hold on
lsim(P(1,:),ref2,time','b')
hold on
plot(time',ref2(:,1),'-.r')
lsim(Pnom(2,:),ref2,time','m')
hold on
lsim(P(2,:),ref2,time','c')
hold on
plot(time',ref2(:,2),'-.m')
hl = title('Time response with $\omega=[0.1sin(\omega t);0.1cos(\omega t)]^{T}$');
set(hl ,'Interpreter','latex');
legend('$P_{1}(s)$','$$\tilde{P}_{1}(s)$$','$\omega_{1}(s)$','$P_{2}(s)$','$$\tilde{P}_{2}(s)$$','$\omega_{2}(s)$','Interpreter', 'LaTeX');
grid on
figure(7);clf
bodemag(WP(1,1),'r');hold on;bodemag(WP2(1,1),'b')
hold on
bodemag(1/WM(1,1),'g')
hl = legend('$\omega_b(s)=11.5$','$\omega_b(s)=20$','$\underline{\sigma}(WM^{-1})$')
set(hl ,'Interpreter','latex')
figure(8);semilogx(w,SV)
figure(9)
lsim(CLHi(1,:),ref2,time','b')
hold on
plot(time',ref2(:,1),'r')
hold on
lsim(CLHi(2,:),ref2,time','g')
hold on
plot(time',ref2(:,2),'m')
hl = title('Time response with $\omega=[0.1sin(\omega t);0.1cos(\omega t)]^{T}$');
set(hl ,'Interpreter','latex');
legend('$$\tilde{P}_{1}(s)$$','$\omega_{1}(s)$','$$\tilde{P}_{2}(s)$$','$\omega_{2}(s)$','Interpreter', 'LaTeX');
grid on
