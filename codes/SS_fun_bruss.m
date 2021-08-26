function Z = SS_fun_bruss(Y)
global a b del del1
x = Y(1);
y = Y(2);
r = Y(3);
u = Y(4);
w = Y(5);
Z(1) = a - (1+b)*x + x^2*y - 1/del*(x-r);
Z(2) = b*x - x^2*y;
Z(3) = 1/del*(x-r) + 1/del1*(u-r);
Z(4) = a - (1+b)*u + u^2*w - 1/del*(u-r);
Z(5) = b*u - u^2*w;
