% test 1D EM 

clear all;
close all;

%% setup the analytic example
omega = 1;  % w = 1
sigma = @(x) x.^2+1;  % 电导率随位置变化
mu    = @(x) cos(x) + 2;  % 磁导率随位置变化

b     = @(x) x.*(x-1).*mu(x);  % 磁场 note b(0) = b(1) = 0;
e     = @(x) cos(2*pi*x);  % 电场

% The system
% 1i*w*b + e'           = s1
% (1/mu * b)' - sigma*e = s2

% fictitious source
s1 = @(x) 1i*omega*b(x) - 2*pi*sin(2*pi*x);
s2 = @(x) 2*x-1 - sigma(x).*e(x);


for n = [8  16  32  64  128  256 512  1024   2048]  
  % setup a random nodal mesh
  h = rand(n,1)*0+1; L = sum(h); h = h/L;  % 这一行等效于h = h / n
  x = [0; cumsum(h)];  % x为节点坐标
  % cell-centered mesh
  xc = x(1:end-1) + diff(x)/2; % xc 为每个cell的中心点坐标
 
  % the sources
  sh = s1(xc);  % 磁场源
  se = s2(x);  % 电场源
  
  % the linear system
  nc = length(sigma(xc));  % [0,1]区间的份数

  G    = ddx(nc);  % n*(n+1)矩阵，对应文章中的G
  Av   = ave(nc);  % n*(n+1)矩阵，对应文章中的Av
  Linv = sdiag(1./h);  % 相当于文章中的L-1或Mf
  Mmu  = sdiag(h./mu(xc));  % 相当于文章中的M(u,f)
  Msig = sdiag(Av'*(sigma(xc).*h));  % 相当于文章中的M(e,sigma)
  M    = sdiag(Av'*h);  % 相当于文章中的Me

  % set the matrix system
  % [1i*w         d/dz] [b] = [s1]
  % [1/mu d/dz  -sigma] [e] = [s2]

  A = [1i*omega*speye(n)  Linv*G; ...
        -G'*Linv'*Mmu      -Msig];  % 按照书上公式，此处应该是Linv,不应该是Linv'
  bb = [sh; M*se]; % M*se相当于书上的Me*s

  
  eh = A\bb;
  subplot(2,1,1)
  plot(xc,real(eh(1:n)),xc,real(b(xc)),'r');
  subplot(2,1,2)
  plot(x,real(eh(1+n:end)),x,real(e(x)),'r');
  
  fprintf('L2 error  %3.2e %3.2e  Linf error %3.2e  %3.2e\n', ...
           norm(Linv\(real(eh(1:n) -b(xc)))), ...
           norm(Linv\(0.5*real(eh(1+n:end-1)+eh(2+n:end))-e(xc))), ...
           norm(real(eh(1:n) -b(xc)),'inf'), ...
           norm(real(eh(1+n:end)-e(x)),'inf'));
           
  pause()
  
end