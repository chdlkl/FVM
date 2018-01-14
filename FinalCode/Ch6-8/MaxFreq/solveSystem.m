function en = solveSystem(Ke,Msig,Mmuinvn,Grad,w,rhs)
% Solve the maxwell system in the freq domain
% using a-phi preconditioner
%

% Ke = CTC + 1i*w*Msig
% Ke*a + STBa*a + 1i*w*Msig*Grad*phi
% 

pcgtol = 1e-7;
% setup preconditioner using Aphi system
STBa = Grad*Mmuinvn*Grad';
Aap = [Ke + STBa,         1i*w*Msig*Grad; ...
       1i*w*Grad'*Msig,   1i*w*Grad'*Msig*Grad];
 
Map = @(x) tril(Aap)\(sdiag(diag(Aap))*(triu(Aap)\x));
%v = randn(size(Aap,1),1); v = v/norm(v);   
%lmax = norm(Aap*v)*2; lmin = 0+1e-12;   
%Map = @(x) SolChebyshev(Aap,x,x*0,5,lmax,lmin);
P1  = [speye(size(Ke,2)); Grad'];
P2  = [speye(size(Ke,2)), Grad];


MM = @(x) P2*(Map(P1*x));
%en = Ke\rhs;
en = rhs*0;
for i=1:size(rhs,2)
     en(:,i) = bicgstab(Ke,rhs(:,i), pcgtol,100,MM);
end
