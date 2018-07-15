Subroutine sdiag(flag)
  use EM_1D
  implicit none
  Integer :: i, j, flag
  real(kind=8) :: temp_sigma(n,1), tempMsig(n+1,1), tempM(n+1,1)
  
  if ( flag == 1 ) then
    Linv = 0.d0
    forall( i = 1:n, j = 1:n, i == j ) Linv(i,j) = 1.d0 / h(i,1)
  else if ( flag == 2 ) then
    Mmu = 0.d0
    forall( i = 1:n, j = 1:n, i == j ) Mmu(i,j) = h(i,1) / mu
  else if ( flag == 3 ) then
    Msig = 0.d0
    tempMsig = matmul( transpose(Av), (sig*h) )
    forall( i = 1:n+1, j = 1:n+1, i == j ) Msig(i,j) = tempMsig(i,1)
  else if ( flag == 4 ) then
    M = 0.d0
    tempM = matmul( transpose(Av), h )
    forall( i = 1:n+1, j = 1:n+1, i == j ) M(i,j) = tempM(i,1)
  else
    write( *,'(1x,a)' ) 'error, please check the argument of sub sdiag!'
  end if
    
End subroutine
  