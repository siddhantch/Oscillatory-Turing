function [x,X] = oscillating_Turing(T,plotting)
% BrussPDE
% Copyright (c) Microsoft Corporation. All rights reserved.
% Licensed under the MIT License.

if nargin<2
  plotting=0;
end

%initialize parameters
%global Xinit Zinit Rinit Uinit Winit

global eps eps1 f f1 q del del1;
eps = 0.215;
f = 1.1;
eps1 = 0.5;
f1 = 0.65;
q = 0.01;
del = 2*eps;%eps/0.365;
del1 = 2*eps1;%eps1/0.365;
Dx = 0.1;
Dz = 0.1;
Dr = 0.1;
Du = 3;
Dw = 100;

%Obtain Steady state for given set of parameters
%OPTIONS = optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
SS = fsolve(@(Y) SS_fun(Y), [0.1;0.1;0.1;0.1;0.1]);

%Mesh initialization
dx = 0.4;
dt=0.2*dx^2/Dw;
L = 102.4;
s=0.0001;
x=0:dx:L;
n=length(x);

%Perturbations to Steady State
X = SS(1)*ones(n)+s*rand(n,n);
Z = SS(2)*ones(n)+s*rand(n,n);
R = SS(3)*ones(n)+s*rand(n,n);
U = SS(4)*ones(n)+s*rand(n,n);
W = SS(5)*ones(n)+s*rand(n,n);

% Xinit = X;
% Zinit = Z;
% Rinit = R;
% Uinit = U;
% Winit = W;


t=0;
i=0;

if plotting
  %ts = [];
  ts = zeros(1,7000);
%   Xmin = [];
%   Xmax = [];
  %Xl = [];
  X1 = zeros(1,7000);
  XA = zeros(1,7000);
  XB = zeros(1,7000);
  XC = zeros(1,7000);
  fh = figure;
  fh.Position = [200 300 1000 400];
end


while t<T
    
  X0 = zeros(n+2,n+2);
  X0(2:n+1,2:n+1) = X;
  X0(1,2:n+1) = X(n,:);
  X0(n+2,2:n+1) = X(1,:);
  X0(2:n+1,1) = X(:,n);
  X0(2:n+1,n+2) = X(:,1);
  Z0 = zeros(n+2,n+2);
  Z0(2:n+1,2:n+1) = Z;
  Z0(1,2:n+1) = Z(n,:);
  Z0(n+2,2:n+1) = Z(1,:);
  Z0(2:n+1,1) = Z(:,n);
  Z0(2:n+1,n+2) = Z(:,1);
  R0 = zeros(n+2,n+2);
  R0(2:n+1,2:n+1) = R;
  R0(1,2:n+1) = R(n,:);
  R0(n+2,2:n+1) = R(1,:);
  R0(2:n+1,1) = R(:,n);
  R0(2:n+1,n+2) = R(:,1);
  U0 = zeros(n+2,n+2);
  U0(2:n+1,2:n+1) = U;
  U0(1,2:n+1) = U(n,:);
  U0(n+2,2:n+1) = U(1,:);
  U0(2:n+1,1) = U(:,n);
  U0(2:n+1,n+2) = U(:,1);
  W0 = zeros(n+2,n+2);
  W0(2:n+1,2:n+1) = W;
  W0(1,2:n+1) = W(n,:);
  W0(n+2,2:n+1) = W(1,:);
  W0(2:n+1,1) = W(:,n);
  W0(2:n+1,n+2) = W(:,1);
  
  
  d2X=4*del2(X0,dx); d2X = d2X(2:n+1,2:n+1);
  d2Z=4*del2(Z0,dx); d2Z = d2Z(2:n+1,2:n+1);
  d2R=4*del2(R0,dx); d2R = d2R(2:n+1,2:n+1);
  d2U=4*del2(U0,dx); d2U = d2U(2:n+1,2:n+1);
  d2W=4*del2(W0,dx); d2W = d2W(2:n+1,2:n+1);

  
  dX = dt*(oreg_F(X,Z,0) - (1/del)*(X-R) + Dx*d2X);
  dZ = dt*(oreg_G(X,Z) + Dz*d2Z);
  dR = dt*((1/del)*(X-R) + (1/del1)*(U-R) + Dr*d2R);
  dU = dt*(oreg_F(U,W,1) - (1/del1)*(U-R) + Du*d2U);
  dW = dt*(oreg_G(U,W) + Dw*d2W);
  
  X = X + dX;
  Z = Z + dZ;
  R = R + dR;
  U = U + dU;
  W = W + dW;
  
  
%   X(1,:)=X(2,:);X(n,:)=X(n-1,:);X(:,1)=X(:,2);X(:,n)=X(:,n-1);
%   Z(1,:)=Z(2,:);Z(n,:)=Z(n-1,:);Z(:,1)=Z(:,2);Z(:,n)=Z(:,n-1);
%   R(1,:)=R(2,:);R(n,:)=R(n-1,:);R(:,1)=R(:,2);R(:,n)=R(:,n-1);
%   U(1,:)=U(2,:);U(n,:)=U(n-1,:);U(:,1)=U(:,2);U(:,n)=U(:,n-1);
%   W(1,:)=W(2,:);W(n,:)=W(n-1,:);W(:,1)=W(:,2);W(:,n)=W(:,n-1);
  
  
  t=t+dt;
  i=i+1;
  
  if plotting 
    if mod(i,100)==0
      %t/T
      subplot(1,2,1)
      figure(1);
      surf(x,x,X);
      set(gca,'layer','top','tickdir','out')
      shading flat
      grid off
      view([0 90])
      %imagesc(X)
      colorbar
      title(sprintf('t = %1.3f',t))
      drawnow;
      
%       figure(2);
%       surf(x,x,U);
%       set(gca,'layer','top','tickdir','out')
%       shading flat
%       grid off
%       view([0 90])
%       %imagesc(X)
%       colorbar
%       title(sprintf('t = %1.3f',t))
%       drawnow;
      
%       ts = [ts t];
        ts(i/100) = t;
%       Xmin = [Xmin min(min(X))];
%       Xmax = [Xmax max(max(X))];
%       Xl = [Xl X(128,128)];
%         XC(i/100) = X(49,229);
%         XA(i/100) = X(64,198);
%         XB(i/100) = X(78,226);
%         X1(i/100) = X(128,128);
%       subplot(1,2,2)
%       figure(3)
%       plot(ts,Xl);
%       figure(4)
%       plot(ts,XA);
%       hold on
%       plot(ts,XB);
%       plot(ts,XC);
      %plot(ts,Xmax,ts,Xmin)
%       ylabel('Range of X')
        t
    end
  else
    if mod(i,1000)==0
      fprintf('.')
    end
  end  
end

      %t/T
%       subplot(1,2,1)
%       figure(1);
%       surf(x,x,X);
%       set(gca,'layer','top','tickdir','out')
%       shading flat
%       grid off
%       view([0 90])
%       %imagesc(X)
%       colorbar
%       title(sprintf('t = %1.3f',t))
%       drawnow;
%       
%       figure(2);
%       surf(x,x,U);
%       set(gca,'layer','top','tickdir','out')
%       shading flat
%       grid off
%       view([0 90])
%      %imagesc(X)
%       colorbar
%       title(sprintf('t = %1.3f',t))
%       drawnow;
%       
%       figure(3)
%       plot(ts,X1);
%       figure(4)
%       plot(ts,XA);
%       hold on
%       plot(ts,XB);
%       plot(ts,XC);

end