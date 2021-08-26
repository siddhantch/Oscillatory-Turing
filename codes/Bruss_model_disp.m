%Auxillary Function
DispersionRelation1 = @(k,J,D)maxk(real(eig(J-k^2*D)),3);
DispersionRelation2 = @(k,J,D)maxk((eig(J-k^2*D)),1,'ComparisonMethod','real');

%Parameters
global a b del del1
a = 2;
b = 4;
del = 0.5;
del1 = 0.5;
Dx = 1;
Dz = 1;
Dr = 1;
Du = 20;
Dw = 20;

%Steady State
SS = fsolve(@(Y) SS_fun_bruss(Y), [0.1;0.1;0;0.1;0.1]);
x = SS(1);
z = SS(2);
r = SS(3);
u = SS(4);
w = SS(5);

%Jacobian and Diffusion matrices
J = zeros(5);
J(1,1) = -(1+b) + 2*x*z - 1/del;
J(1,2) = x^2;
J(1,3) = 1/del;
J(2,1) = b-2*x*z;
J(2,2) = -x^2;
J(3,1) = 1/del;
J(3,3) = -1/del - 1/del1;
J(3,4) = 1/del1;
J(4,3) = 1/del1;
J(4,4) = -(1+b) + 2*u*w - 1/del;
J(4,5) = u^2;
J(5,4) = b-2*u*w;
J(5,5) = -u^2;


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