function F = SS_fun(Y)
format long
global del del1 %f f1 q eps eps1 ;
F = zeros(3,1);
F(1) = oreg_F(Y(1),Y(2),0) - (1/del)*(Y(1) - Y(3));
F(2) = oreg_G(Y(1),Y(2));
F(3) = (1/del)*(Y(1)-Y(3)) + (1/del1)*(Y(4)-Y(3));
F(4) = oreg_F(Y(4),Y(5),1) - (1/del1)*(Y(4) - Y(3));
F(5) = oreg_G(Y(4),Y(5));
% F(1) = 2*Y(1)^2 - Y(1) + 2*f*Y(1)*((Y(1)-q)/(Y(1)+q)) - Y(2);
% F(2) = (Y(1) - Y(2))/eps + (Y(3) - Y(2))/eps1;
% F(3) = 2*Y(3)^2 - Y(3) + 2*f1*Y(3)*((Y(3)-q)/(Y(3)+q)) - Y(2);
end