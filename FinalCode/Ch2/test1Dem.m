% test 1D EM 

clear all;
close all;

%% setup the analytic example
omega = 1;  % w = 1
sigma = @(x) x.^2+1;  % �絼����λ�ñ仯
mu    = @(x) cos(x) + 2;  % �ŵ�����λ�ñ仯

b     = @(x) x.*(x-1).*mu(x);  % �ų� note b(0) = b(1) = 0;
e     = @(x) cos(2*pi*x);  % �糡

% The system
% 1i*w*b + e'           = s1
% (1/mu * b)' - sigma*e = s2

% fictitious source
s1 = @(x) 1i*omega*b(x) - 2*pi*sin(2*pi*x);
s2 = @(x) 2*x-1 - sigma(x).*e(x);


for n = [8  16  32  64  128  256 512  1024   2048]  
  % setup a random nodal mesh
  h = rand(n,1)*0+1; L = sum(h); h = h/L;  % ��һ�е�Ч��h = h / n
  x = [0; cumsum(h)];  % xΪ�ڵ�����
  % cell-centered mesh
  xc = x(1:end-1) + diff(x)/2; % xc Ϊÿ��cell�����ĵ�����
 
  % the sources
  sh = s1(xc);  % �ų�Դ
  se = s2(x);  % �糡Դ
  
  % the linear system
  nc = length(sigma(xc));  % [0,1]����ķ���

  G    = ddx(nc);  % n*(n+1)���󣬶�Ӧ�����е�G
  Av   = ave(nc);  % n*(n+1)���󣬶�Ӧ�����е�Av
  Linv = sdiag(1./h);  % �൱�������е�L-1��Mf
  Mmu  = sdiag(h./mu(xc));  % �൱�������е�M(u,f)
  Msig = sdiag(Av'*(sigma(xc).*h));  % �൱�������е�M(e,sigma)
  M    = sdiag(Av'*h);  % �൱�������е�Me

  % set the matrix system
  % [1i*w         d/dz] [b] = [s1]
  % [1/mu d/dz  -sigma] [e] = [s2]

  A = [1i*omega*speye(n)  Linv*G; ...
        -G'*Linv'*Mmu      -Msig];  % �������Ϲ�ʽ���˴�Ӧ����Linv,��Ӧ����Linv'
  bb = [sh; M*se]; % M*se�൱�����ϵ�Me*s

  
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