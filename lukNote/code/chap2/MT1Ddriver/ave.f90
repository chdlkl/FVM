Subroutine ave()
  use EM_1D
  Implicit none
  Integer :: i, j
  
  Av = 0.d0
  forall( i = 1:n, j = 1:n+1, i == j ) Av(i,j) = 5.d-1
  forall( i = 1:n, j = 1:n+1, j-i == 1 ) Av(i,j) = 5.d-1

End subroutine ave