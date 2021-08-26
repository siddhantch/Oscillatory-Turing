function Ad = del_m(A0)
n = size(A0);
dx = 0.4;
n = n(1);
Ad = zeros(n);
for i = 2:n-1
    for j = 2:n-1
        Ad(i,j) = (-20*A0(i,j) + A0(i+1,j+1) + A0(i-1,j+1) + A0(i+1,j-1) + A0(i-1,j-1) + 4 * (A0(i+1,j) + A0(i,j+1) + A0(i-1,j) + A0(i,j-1))) / (6 * dx * dx);
    end
end
