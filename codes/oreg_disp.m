DispersionRelation = @(k,J,D)maxk(real(eig(J-k^2*D)),1);
DispersionRelation2 = @(k,J,D)maxk((eig(J-k^2*D)),1,'ComparisonMethod','real');

%Parameters
global eps f q Da Db
eps=0.23;f=1.4;q = 0.01; Da=0.17;Db=0.17;

%Equilibria
SS = fzero(@(Y) oreg_SS(Y),0.5);
Ass=SS;Bss=SS; 

% Classical oregonator
J = [1/eps*(1 - 2*Ass - f*Bss*(1/(Ass+q) - (Ass-q)/(Ass+q)^2)), 1/eps*(-f*(Ass-q)/(Ass+q));
    1,-1];
D = diag([Da Db]);

plot(0:0.01:2,arrayfun(@(k)DispersionRelation(k,J,D),0:0.01:2));
hold on;
plot(0:0.01:2,zeros(1,201),'b--');
ylabel('Re(\lambda)');
xlabel('k');
%ylim([-0.45 0.45]);
%xlim([0 1.7]);
figure(2);
plot(0:0.01:2,imag(arrayfun(@(k)DispersionRelation2(k,J,D),0:0.01:2)));