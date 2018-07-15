Subroutine ddx()
  use EM_1D
  implicit none
  Integer :: i, j
  
  G = 0.d0
  forall( i = 1:n, j = 1:n+1, i == j ) G(i,j) = -1.d0
  forall( i = 1:n, j = 1:n+1, j-i == 1 ) G(i,j) = 1.d0
  
End subroutine ddx