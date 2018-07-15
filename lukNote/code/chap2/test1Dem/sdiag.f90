Subroutine sdiag(flag)  !.. ∂‘Ω«æÿ’Û
  use EM_1D
  implicit none
  Integer :: i, j, flag
  real(kind=8) :: temp, tempMmu(n), temp_sigma(n,1), tempMsig(n+1,1)
  real(kind=8) :: temp_h(n,1), tempM(n+1,1)
  
  if ( flag == 1 ) then
    temp = 1.d0 / h
    Linv = 0.d0
    forall( i = 1:n, j = 1:n, i == j ) Linv(i,j) = temp
  else if ( flag == 2 ) then
    Mmu = 0.d0
    tempMmu = mu()
    forall( i = 1:n, j = 1:n, i == j ) Mmu(i,j) = h / tempMmu(i)
  else if ( flag == 3 ) then
    Msig = 0.d0
    temp_sigma(:,1) = sigma()
    temp_sigma = h * temp_sigma
    tempMsig = matmul( transpose(Av), temp_sigma )
    forall( i = 1:n+1, j = 1:n+1, i == j ) Msig(i,j) = tempMsig(i,1)
  else if ( flag == 4 ) then
    M = 0.d0
    temp_h(:,1) = h
    tempM = matmul( transpose(Av), temp_h )
    forall( i = 1:n+1, j = 1:n+1, i == j ) M(i,j) = tempM(i,1)
  else
    write( *,'(1x,a)' ) 'error, please check the argument of sub sdiag!'
    write( *,'(1x,a)' ) 'please input Enter and proram stop!'
    read( *,* )
  end if

End subroutine sdiag