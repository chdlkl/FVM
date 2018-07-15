Module EM_1D
  implicit none
  Integer :: n  !.. n为单元数 
  real(kind=8), allocatable :: G(:,:), Av(:,:)  !.. G为差分矩阵, Av是平均矩阵
  real(kind=8), allocatable :: x(:), xc(:)  !.. x为长度为n+1的节点坐标; xc为每个cell的中心点坐标，长度为n
  real(kind=8), allocatable :: Linv(:,:), Mmu(:,:) !.. Linv相当于文章中的L-1或Mf，生成1./h为主对角线元素的对称矩阵; Mmu相当于文章中的M(u,f)
  real(kind=8), allocatable :: Msig(:,:), M(:,:)  !.. Msig相当于文章中的M(e,sigma); M相当于文章中的Me
  real(kind=8), allocatable :: se(:,:)   !.. se为电场源，长度为n+1
  real(kind=8), allocatable :: tempbb(:,:)
  complex(kind=8), allocatable :: sh(:,:)  !.. sh为磁场源，长度为n
  complex(kind=8), allocatable :: A(:,:), bb(:,:)  !.. A为左边系数矩阵; bb为右端项，大小为2*n+1
  real(kind=8), parameter :: startPoint = 0.d0, endPoint = 1.d0  !.. 计算区间为[0,1]
  real(kind=8), parameter :: pi = acos(-1.d0), omiga = 1.d0 
  real(kind=8) :: h  !.. h为单元长度。
  
Contains
  !.. set up the analytic example
  function sigma() 
    implicit none
    real(kind=8), allocatable :: sigma(:)
    
    If ( allocated( sigma ) ) then
      write( *,'(1x,a)' ) 'sigma has already been allocated!'
    else
      allocate( sigma(size(xc)) )
    End if
    
    sigma = xc**2 + 1.d0
    
  end function sigma
  
  function sigma_e()  !.. set up sigma
    implicit none
    real(kind=8), allocatable :: sigma_e(:)
    
    If ( allocated( sigma_e ) ) then
      write( *,'(1x,a)' ) 'sigma_e has already been allocated!'
    else
      allocate( sigma_e(size(x)) )
    End if
    
    sigma_e = x**2 + 1.d0
    
  end function sigma_e
  
  function mu()  !.. set up mu
    implicit none
    real(kind=8), allocatable :: mu(:)
    
    If ( allocated( mu ) ) then
      write( *,'(1x,a)' ) 'mu has already been allocated!'
    else
      allocate( mu(size(xc)) )
    End if
    
    mu = cos(xc) + 2.d0
    
  end function mu
  
  !.. set up b and e
  function b()  
    implicit none
    real(kind=8), allocatable :: b(:)
    
    If ( allocated( b ) ) then
      write( *,'(1x,a)' ) 'b has already been allocated!'
    else
      allocate( b(size(xc)) )
    End if
    
    b = xc * (xc-1.d0) * mu()  !.. note b(startPoint) = b(endPoint) = 0  要满足边界条件
    
  end function b
  
  function e()
    implicit none
    real(kind=8), allocatable :: e(:)
    
    If ( allocated( e ) ) then
      write( *,'(1x,a)' ) 'e has already been allocated!'
    else
      allocate( e(size(x)) )
    End if
    
    e = cos( 2.d0 * pi * x )
    
  end function e
    
  !.. the system 
  !.. i*w*b + e'            = s1
  !.. (1/mu * b)' - sigma*e = s2
  function s1(mm)  !.. 磁场源
    implicit none
    Integer, intent(in) :: mm
    complex(kind=8), allocatable :: s1(:)
  
    If ( allocated( s1 ) ) then
      write( *,'(1x,a)' ) 's1 has already been allocated!'
    else
      allocate( s1(mm) )
    End if
    
    s1 = cmplx( -2.d0*pi*sin(2.d0*pi*xc), omiga*b() )
  
  end function s1
  
  function s2(mm)  !.. 电场源
    implicit none
    Integer, intent(in) :: mm
    real(kind=8), allocatable :: s2(:)
  
    If ( allocated( s2 ) ) then
      write( *,'(1x,a)' ) 's2 has already been allocated!'
    else
      allocate( s2(mm) )
    End if
  
    s2 = 2.d0*x-1.d0 - sigma_e()*e()
  
  end function s2  
  
  Subroutine output() !.. output data
    implicit none
    Integer :: i, fileid
    real(kind=8) :: tempB(n), tempE(n+1)
  
    tempB(:) = b()
    open( newunit = fileid, file = 'b.dat' )
    do i = 1, n
      write( fileid,* ) xc(i), real( bb(i,1) ), real( tempB(i) )
    end do
    close( fileid )
  
    tempE(:) = e()
    open( newunit = fileid, file = 'e.dat' )
    do i = n+1, 2*n+1
      write( fileid,* ) x(i-n), real( bb(i,1) ), real( tempE(i-n) )
    end do
    close( fileid )
  End subroutine output
  
End module EM_1D

Program test1Dem
  use EM_1D
  use lapack95
  implicit none
  Integer :: i, j
  Integer, parameter :: start = 1, finish = 9
  Character(len=20) date,time
  
  write( *,'(1x,a)' ) 'start time as flow:'
  call date_and_time(date,time)  
  write( *,"('date:',a8,/,'time:',a2,':',a2,':',a2)" ) date, time(1:2), time(3:4), time(5:6)
  Do i = start, finish  !.. 计算7批，每次的单元数目为2^(i-1)*8
    n = 2**(i-1) * 8
    write( *,'(1x,a,g0,a,g0)' ) '正在计算第 ', i, ' 次，步长为 ', n
    allocate( x(n+1), xc(n) )
    allocate( se(n+1,1), sh(n,1) )
    allocate( G(n,n+1), Av(n,n+1), Linv(n,n), Mmu(n,n), Msig(n+1,n+1), M(n+1,n+1) )
    allocate( A(2*n+1,2*n+1), bb(2*n+1,1) )
    allocate( tempbb(n+1,1) )
    
    !.. set up a random nodal mesh
    h = (endPoint-startPoint) / n
    do j = 1, n+1  !.. 获取节点坐标
      x(j) = startPoint + dble(j-1) * h
    end do
    
    !.. cell-centered mesh
    do j = 1, n
      xc(j) = ( x(j) + x(j+1) ) / 2.d0  !.. 获取中心点坐标
    end do
    
    !.. the sources
    sh(:,1) = s1(n)
    se(:,1) = s2(n+1)
    call ddx
    call ave
    call sdiag(1)
    call sdiag(2)
    call sdiag(3)
    call sdiag(4)
    
    !.. set the matrix system
    !.. [i*w           d/dz] [b] - s1
    !.. [1/mu  d/dz  -sigma] [e] - s2
    !.. 构建系数矩阵，大小为(2*n+1,2*n+1)
    A = 0.d0
    !.. 先构造左上角大小为n*n的子矩阵
    forall( i = 1:n, j = 1:n, i == j ) A(i,j) = cmplx( 0.d0, omiga )
    !.. 构造右上角大小为n*(n+1)的子矩阵
    A(1:n,n+1:2*n+1) = cmplx( matmul( Linv, G ), 0.d0 )
    !.. 构造左下角大小为(n+1)*n的子矩阵
    A(n+1:2*n+1,1:n) = cmplx( -1.d0 * matmul( matmul( transpose(G),transpose(Linv) ), Mmu ), 0.d0 )
    !.. 构造右下角大小为(n+1)*(n+1)的子矩阵
    A(n+1:2*n+1,n+1:2*n+1) = cmplx( -1.d0*Msig, 0.d0 )

    !.. 右端项
    bb(1:n,1) = sh(:,1)
    tempbb = matmul( M, se )
    bb(n+1:2*n+1,1) = tempbb(:,1)
    
    !.. 求解
    call gesv( A, bb )  !.. 调用gesv之后，bb为求出的解

    if ( i == finish )  call output()
    deallocate( x, xc )
    deallocate( se, sh )
    deallocate( G, Av, Linv, Mmu, Msig, M )
    deallocate( A, bb ) 
    deallocate( tempbb )
  End do
  
  write( *,'(1x,a)' ) 'stop time as flow:'
  call date_and_time(date,time)  
  write( *,"('date:',a8,/,'time:',a2,':',a2,':',a2)" ) date, time(1:2), time(3:4), time(5:6)
  
End program test1Dem  