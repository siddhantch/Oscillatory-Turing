function OF = oreg_F(x,z,flag)
global eps eps1 f f1 q;
e = eps;
ef = f;
if flag == 1
    e = eps1;
    ef = f1;
end
OF = (1/e).*(x - x.^2 - ef*z.*((x-q)./(x+q)));
end