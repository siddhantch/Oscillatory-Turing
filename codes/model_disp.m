%Auxillary Function
DispersionRelation1 = @(k,J,D)maxk(real(eig(J-k^2*D)),3);
DispersionRelation2 = @(k,J,D)maxk((eig(J-k^2*D)),1,'ComparisonMethod','real');

%Parameters
global eps eps1 f f1 q del del1;
eps = 0.215;
f = 1.1;
eps1 = 0.5;
f1 = 0.65;
q = 0.01;
del = 2*eps;
del1 = 2*eps1;
Dx = 0.1;
Dz = 0.1;
Dr = 0.1;
Du = 3;
Dw = 100;

%Steady State
SS = fsolve(@(Y) SS_fun(Y), [0.1;0.1;0;0.1;0.1]);
x = SS(1);
z = SS(2);
r = SS(3);
u = SS(4);
w = SS(5);

% x = 0.147;
% z = 0.147;
% r = 0.18;
% u = 0.253;
% w = 0.253;

%Jacobian and Diffusion matrices
J = zeros(5);
J(1,1) = 1/eps*(1-2*x-f*z*(1/(x+q)-(x-q)/(x+q)^2))- 1/del;
J(1,2) = 1/eps*(-f*(x-q)/(x+q));
J(1,3) = 1/del;
J(2,1) = 1;
J(2,2) = -1;
J(3,1) = 1/del;
J(3,3) = -1/del - 1/del1;
J(3,4) = 1/del1;
J(4,3) = 1/del1;
J(4,4) = 1/eps1*(1-2*u-f1*w*(1/(u+q)-(u-q)/(u+q)^2))- 1/del1;
J(4,5) = 1/eps1*(-f1*(u-q)/(u+q));
J(5,4) = 1;
J(5,5) = -1;


D = diag([Dx Dz Dr Du Dw]);

%DispersionRelation(1.9,J,D)
%plot(k,arrayfun(@(k)DispersionRelation(k,J,D),k),'LineWidth',1);
A0 = [];
B0 = [];
for k = 0:0.01:2
    A = DispersionRelation1(k,J,D);
    B = imag(DispersionRelation2(k,J,D));
    A0 = [A0 A];
    B0 = [B0 B];
end
k = 0:0.01:2;
figure(1);
plot(k,A0);
hold on
plot(k,zeros(1,201),'b--'); 
xlabel('k');
ylabel('Re(\lambda)');
% ylim([-0.3 0.3]);

% xlim([0 0.8]);
figure(2);
plot(k,B0);
hold on
xlabel('k');
ylabel('img(\lambda)');

[m, i] = max(max(A0));
imag_m = B0(i);