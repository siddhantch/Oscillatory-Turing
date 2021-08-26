function [x,A] = oregsim(T,plotting)
% BrussPDE
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

if nargin<2
  plotting=0;
end

%initialize parameters
global eps f q;
eps = 0.5;
f = 1.1;
q = 0.01;
Da = 0.5;
Db = 12;

%Obtain Steady state for given set of parameters
SS = fzero(@(Y) oreg_SS(Y),0.05);

%Mesh initialization
dx=0.2;
dt=0.2*dx^2/Db;
L=51.2;
s=0.01;
x=0:dx:L;
n=length(x);

%Perturbations to Steady State
A = SS*ones(n)+s*rand(n,n);
B = SS*ones(n)+s*rand(n,n);
t=0;
i=0;

if plotting
  ts = [];
  Amin = [];
  Amax = [];
  fh = figure;
  fh.Position = [200 300 1000 400];
end


while t<T
    
  A0 = zeros(n+2,n+2);
  A0(2:n+1,2:n+1) = A;
  A0(1,2:n+1) = A(n,:);
  A0(n+2,2:n+1) = A(1,:);
  A0(2:n+1,1) = A(:,n);
  A0(2:n+1,n+2) = A(:,1);
  B0 = zeros(n+2,n+2);
  B0(2:n+1,2:n+1) = B;
  B0(1,2:n+1) = B(n,:);
  B0(n+2,2:n+1) = B(1,:);
  B0(2:n+1,1) = B(:,n);
  B0(2:n+1,n+2) = B(:,1);
  
  
  d2A=4*del2(A0,dx); d2A = d2A(2:n+1,2:n+1);
  d2B=4*del2(B0,dx); d2B = d2B(2:n+1,2:n+1);

  
  dA = dt*(oreg_F(A,B,0) + Da*d2A);
  dB = dt*(oreg_G(A,B) + Db*d2B);
  
  
  A = A + dA;
  B = B + dB;
 
  
  
%    A(1,:)=A(2,:);A(n,:)=A(n-1,:);A(:,1)=A(:,2);A(:,n)=A(:,n-1);
%    B(1,:)=B(2,:);B(n,:)=B(n-1,:);B(:,1)=B(:,2);B(:,n)=B(:,n-1);

  
  
  t=t+dt;
  i=i+1;
  
  if plotting 
    if mod(i,100)==0
      t/T;
      subplot(1,2,1)
      surf(x,x,A);
      set(gca,'layer','top','tickdir','out')
      shading flat
      grid off
      view([0 90])
      %imagesc(X)
      colorbar
      title(sprintf('t = %1.3f',t))
      drawnow;
      ts = [ts t];
      Amin = [Amin min(min(A))];
      Amax = [Amax max(max(A))];
      subplot(1,2,2)
      plot(ts,Amax,ts,Amin)
      ylabel('Range of A')
    end
  else
    if mod(i,1000)==0
      fprintf('.')
    end
  end  
end

end