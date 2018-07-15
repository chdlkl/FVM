Module EM_1D
  implicit none
  Integer :: n  !.. nΪ��Ԫ�� 
  real(kind=8), allocatable :: G(:,:), Av(:,:)  !.. GΪ��־���, Av��ƽ������
  real(kind=8), allocatable :: x(:), xc(:)  !.. xΪ����Ϊn+1�Ľڵ�����; xcΪÿ��cell�����ĵ����꣬����Ϊn
  real(kind=8), allocatable :: Linv(:,:), Mmu(:,:) !.. Linv�൱�������е�L-1��Mf������1./hΪ���Խ���Ԫ�صĶԳƾ���; Mmu�൱�������е�M(u,f)
  real(kind=8), allocatable :: Msig(:,:), M(:,:)  !.. Msig�൱�������е�M(e,sigma); M�൱�������е�Me
  real(kind=8), allocatable :: se(:,:)   !.. seΪ�糡Դ������Ϊn+1
  real(kind=8), allocatable :: tempbb(:,:)
  complex(kind=8), allocatable :: sh(:,:)  !.. shΪ�ų�Դ������Ϊn
  complex(kind=8), allocatable :: A(:,:), bb(:,:)  !.. AΪ���ϵ������; bbΪ�Ҷ����СΪ2*n+1
  real(kind=8), parameter :: startPoint = 0.d0, endPoint = 1.d0  !.. ��������Ϊ[0,1]
  real(kind=8), parameter :: pi = acos(-1.d0), omiga = 1.d0 
  real(kind=8) :: h  !.. hΪ��Ԫ���ȡ�
  
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
    
    b = xc * (xc-1.d0) * mu()  !.. note b(startPoint) = b(endPoint) = 0  Ҫ����߽�����
    
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
  function s1(mm)  !.. �ų�Դ
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
  
  function s2(mm)  !.. �糡Դ
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
  Do i = start, finish  !.. ����7����ÿ�εĵ�Ԫ��ĿΪ2^(i-1)*8
    n = 2**(i-1) * 8
    write( *,'(1x,a,g0,a,g0)' ) '���ڼ���� ', i, ' �Σ�����Ϊ ', n
    allocate( x(n+1), xc(n) )
    allocate( se(n+1,1), sh(n,1) )
    allocate( G(n,n+1), Av(n,n+1), Linv(n,n), Mmu(n,n), Msig(n+1,n+1), M(n+1,n+1) )
    allocate( A(2*n+1,2*n+1), bb(2*n+1,1) )
    allocate( tempbb(n+1,1) )
    
    !.. set up a random nodal mesh
    h = (endPoint-startPoint) / n
    do j = 1, n+1  !.. ��ȡ�ڵ�����
      x(j) = startPoint + dble(j-1) * h
    end do
    
    !.. cell-centered mesh
    do j = 1, n
      xc(j) = ( x(j) + x(j+1) ) / 2.d0  !.. ��ȡ���ĵ�����
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
    !.. ����ϵ�����󣬴�СΪ(2*n+1,2*n+1)
    A = 0.d0
    !.. �ȹ������ϽǴ�СΪn*n���Ӿ���
    forall( i = 1:n, j = 1:n, i == j ) A(i,j) = cmplx( 0.d0, omiga )
    !.. �������ϽǴ�СΪn*(n+1)���Ӿ���
    A(1:n,n+1:2*n+1) = cmplx( matmul( Linv, G ), 0.d0 )
    !.. �������½Ǵ�СΪ(n+1)*n���Ӿ���
    A(n+1:2*n+1,1:n) = cmplx( -1.d0 * matmul( matmul( transpose(G),transpose(Linv) ), Mmu ), 0.d0 )
    !.. �������½Ǵ�СΪ(n+1)*(n+1)���Ӿ���
    A(n+1:2*n+1,n+1:2*n+1) = cmplx( -1.d0*Msig, 0.d0 )

    !.. �Ҷ���
    bb(1:n,1) = sh(:,1)
    tempbb = matmul( M, se )
    bb(n+1:2*n+1,1) = tempbb(:,1)
    
    !.. ���
    call gesv( A, bb )  !.. ����gesv֮��bbΪ����Ľ�

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