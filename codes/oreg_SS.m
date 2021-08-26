function op = oreg_SS(y)
global eps q f;
op = 1/eps*(y - y^2 - f*y*(y-q)/(y+q));
end