clear all;clc
%% Nominal Plant Dynamics
s=tf('s');
N = {[1 -100],[10 10];[-10 -10],[1 -100]};
D = [1 0 100];
Pnom = tf(N,D);
Pnom = ss(Pnom,'min')
SystemPole = pole(Pnom);
SystemZero = tzero(Pnom);
nmeas = 2; ncon = 2;
Ms = 2; A = 1e-2; wb = 11.5; 
Ms2 = 8;wb2 = 20;
wP2 = tf([1/Ms2 wb2], [1 wb2*A]);
WP2 = eye(2)*wP2; WP2 = ss(WP2) ;
%% Plant with uncertainty
for i = 1:2
    k(i) = ureal('k',1,'Range',[0.8,1.2]);
    theta(i) = ureal('theta',0.01,'Range',[0,.02]);
    f(i) = k(i)*((-.5*theta(i)*s)+1)/((.5*theta(i)*s)+1);
end
P = [f(1) 0;0 f(2)]*Pnom;
bb = (P-Pnom)*inv(P);
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
wP = tf([1/Ms wb], [1 wb*A]);
WP = eye(2)*wP; WP = ss(WP) ;
systemnames = 'Pnom WP2 WM';
inputvar = '[w(2);u(2)]';
outputvar = '[WP2;WM;-w-Pnom]';
input_to_Pnom= '[u]';
input_to_WP = '[w+Pnom]';
input_to_WM = ' [Pnom]';
G = sysic
[Khiold,CLhi,ghiold,hiinfo] = hinfsyn(G,nmeas,ncon);
Fhiold=loopsens(Pnom,Khiold);
[SVold,wold] = sigma(CLhi);
%% %% %% Trade off gamma based on Gmax and Gmin
%Ms = 8 dan wb = 20
systemnames = 'Pnom WP WM';
inputvar = '[w(2);u(2)]';
outputvar = '[WP;WM;-w-Pnom]';
input_to_Pnom= '[u]';
input_to_WP = '[w+Pnom]';
input_to_WM = ' [Pnom]';
Gnew = sysic
[Khi,CLHi,ghi,hiinfo] = hinfsyn(Gnew,nmeas,ncon,0.6719)
[KhiF,CLHiF,ghiF,hiinfoF] = hinfsyn(Gnew,nmeas,ncon,1.5)
lft_norm = hinfnorm(lft(Gnew,Khi));
[SV,w] = sigma(CLHi);
[SVF,wF] = sigma(CLHiF);
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
figure(3);sigmaplot(bb)
hold on;sigmaplot(WM,'r');
grid on
title('Weights for Multiplicative Uncertainty')
legend('$$\bar{\sigma}(\tilde{P}-P)P^{-1}$$','$$Wm(s)=\frac{0.021s+0.2}{0.0091s+1}$$','Interpreter','latex')
figure(3);sigmaplot(bb,'b')
hold on;sigmaplot(WM,'r');
grid on
title('Weights for Multiplicative Uncertainty')
legend('$$\bar{\sigma}(\tilde{P}-P)P^{-1}$$','$$Wm(s)=\frac{0.021s+0.2}{0.0091s+1}$$','Interpreter','latex')
figure(4);bodemag(wM,'b');hold on;bodemag(wM2,'r');
legend('r(\infty) =2.3','r(\infty) = 100')
figure(5);clf
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
figure(6);clf;
lsim(Pnom(1,:),ref2,time','r');hold on
lsim(P(1,:),ref2,time','b');hold on
plot(time',ref2(:,1),'-.r');hold on
lsim(Pnom(2,:),ref2,time','m');hold on
lsim(P(2,:),ref2,time','c');hold on
plot(time',ref2(:,2),'-.m')'
hl = title('Time response with $\omega=[0.1sin(\omega t);0.1cos(\omega t)]^{T}$');
set(hl ,'Interpreter','latex');
legend('$P_{1}(s)$','$$\tilde{P}_{1}(s)$$','$\omega_{1}(s)$','$P_{2}(s)$','$$\tilde{P}_{2}(s)$$','$\omega_{2}(s)$','Interpreter', 'LaTeX');
grid on
figure(7);clf
bodemag(WP(1,1),'r');hold on;bodemag(WP2(1,1),'b')
hold on/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*/*
bodemag(1/WM(1,1),'g')
hl = legend('$\omega_b(s)=11.5$','$\omega_b(s)=20$','$\underline{\sigma}(WM^{-1})$')
set(hl ,'Interpreter','latex')
figure(8);clf;semilogx(w,SV)
figure(9);clf
lsim(CLHi(1,:),ref,time','b')
hold on
plot(time',ref(:,1),'r')
hold on
lsim(CLHi(2,:),ref,time','g')
hold on
plot(time',ref(:,2),'m')
hl = title('Time response with $\omega=[0.1step;0]^{T}$');
set(hl ,'Interpreter','latex');
legend('$$\tilde{P}_{1}(s)$$','$\omega_{1}(s)$','$$\tilde{P}_{2}(s)$$','$\omega_{2}(s)$','Interpreter', 'LaTeX');
grid on
figure(10);clf
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
figure(11);semilogx(w,SV(1,:),'r')
hold on
figure(11);semilogx(wF,SVF(1,:),'b')
grid on
ylabel('Magnitude')
xlabel('Freq[rad/s]')
title('Close loop TF $\bar{\sigma}(L_f(G,K))$','Interpreter','latex')
legend('$\gamma=0.6719$','$\gamma=1.5$','Interpreter','latex')
figure(12);clf;sigmaplot(Khi,'b')
grid on
text(1,30,'$$\bar{\sigma}K_{hi}$$','Interpreter','latex');
hold on;
text(15,-12,'$$\underline{\sigma}K_{hi}$$','Interpreter','latex')
