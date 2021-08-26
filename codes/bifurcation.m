global eps f q
DispersionRelation = @(k,J,D)max(real(eig(J-k^2*D)));

%Turing Pattern
Da = 3;
Db = 100;
q = 0.01;
Tur = zeros(1,200);
i = 1;
for f = linspace(0.4,2.4,200)
    for eps = logspace(-3,1,300)
        SS = fzero(@(Y) oreg_SS(Y),0.5);
        Ass=SS;Bss=SS; 
        J = [1/eps*(1 - 2*Ass - f*Bss*(1/(Ass+q) - (Ass-q)/(Ass+q)^2)), 1/eps*(-f*(Ass-q)/(Ass+q));1,-1];
        D = diag([Da Db]);
        if max(arrayfun(@(k)DispersionRelation(k,J,D),0:0.01:2)) < 0
            Tur(i) = eps;
            break;
        end
    end
    i = i+1;
end

%Hopf Instability
Da = 0.1;
Db = 0.1;
q = 0.01;
Hop = zeros(1,200);
i = 1;
for f = linspace(0.4,2.4,200)
    for eps = logspace(-3,1,300)
        SS = fzero(@(Y) oreg_SS(Y),0.5);
        Ass=SS;Bss=SS; 
        J = [1/eps*(1 - 2*Ass - f*Bss*(1/(Ass+q) - (Ass-q)/(Ass+q)^2)), 1/eps*(-f*(Ass-q)/(Ass+q));1,-1];
        D = diag([Da Db]);
        if max(arrayfun(@(k)DispersionRelation(k,J,D),0:0.01:2)) < 0
            Hop(i) = eps;
            break;
        end
    end
    i = i+1;
end

plot(linspace(0.4,2.4,200),Tur);
hold on;
plot(linspace(0.4,2.4,200),Hop);
xlim([0.4 2.3]);
ylim([5*10^-2 10])
set(gca, 'YScale', 'log');