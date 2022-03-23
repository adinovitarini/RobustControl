clear all;clc
theta = ureal('theta',.05,'range',[0,.1]);
s=tf('s');
G = tf(1,[1 1],'InputDelay',0.05)
%Nominal 
Go = tf(1,[1 1])
%Weight 
WDelta = (G/Go)-1
W = makeweight(.2,.25,2.5)
%Transform ss ke TF 
[W_num,W_denum] = ss2tf(W.A,W.B,W.C,W.D);
W_tf = tf(W_num,W_denum)
Delta = ultidyn('Delta',[1 1])
G_test = Go*(1+W*Delta)
figure(1);clf
bodemag(WDelta)
hold on
bodemag(W_tf,'r')
title('Weights for Multiplicative Uncertainty')
legend('$W(s)\Delta(s)$','$W(s)=\frac{2.5s+0.1169}{s+0.5846}$','Interpreter','latex')
figure(2);clf
subplot(121)
bode(G,'g')
title('G(s) without model uncertainty')
subplot(122)
bode(G_test,'k')
title('G(s) with model uncertainty')
hold on
figure(3)
bodemag(G_test,{0.5,5},'k')
title('All stability margins')
