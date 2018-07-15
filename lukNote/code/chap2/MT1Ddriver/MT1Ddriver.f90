Module EM_1D
  implicit none
  Integer, parameter :: n = 1025  !.. nΪ��Ԫ�� 
  Integer, parameter :: minFre = 1, maxFre = 4, numFre = 32  !.. Ƶ����
  real(kind=8) :: omiga(numFre), skin(numFre), h(n,1), sig(n,1)  !.. hΪ��Ԫ����, sigΪ�絼�ʾ���
  real(kind=8), allocatable :: G(:,:), Av(:,:)  !.. GΪ��־���, Av��ƽ������
  real(kind=8), allocatable :: z(:), zc(:)  !.. zΪ����Ϊn+1�Ľڵ�����; zcΪÿ��cell�����ĵ����꣬����Ϊn
  real(kind=8), allocatable :: Linv(:,:), Mmu(:,:) !.. Linv�൱�������е�L-1��Mf������1./hΪ���Խ���Ԫ�صĶԳƾ���; Mmu�൱�������е�M(u,f)
  real(kind=8), allocatable :: Msig(:,:), M(:,:)  !.. Msig�൱�������е�M(e,sigma); M�൱�������е�Me
  real(kind=8), parameter :: pi = acos(-1.d0), mu = 4.d0*pi*1.d-7, sig0 = 1.d-2  !.. sig0Ϊ�����絼�� 
  complex(kind=8), allocatable :: A(:,:), AA(:,:), bb(:,:)   !.. AΪ���ϵ������;bbΪ�Ҷ����СΪ2*n
  complex(kind=8) :: b(n), e(n), d(numFre)  !.. d��ʾÿ��Ƶ���ڵر�����ֵ
  real(kind=8) :: L  !.. LΪ��������
  
Contains
  Subroutine get_para()  !.. �õ�omiga��skin
    implicit none
    Integer :: i
    real(kind=8) :: temp
    
    do i = 1, numFre
      temp = minFre + dble(i-1)*( dble(maxFre-minFre)/(numFre-1) )
      omiga(i) = 10**temp
      skin(i) = sqrt( 2.d0/omiga(i)/mu/sig0 )
    end do
    
  End subroutine get_para
  
  Subroutine output()  !.. output data
    implicit none
    Integer :: i, fileid
    
    open( newunit = fileid, file = 'Pa_and_angle.dat' )
    do i = 1, numFre
      write( fileid,* ) omiga(i), mu/omiga(i)*abs(d(i))**2, atan( imag(d(i))/real(d(i)) ) + pi  !.. �˴���piֻ�ǵ�����ת������
    end do
    close( fileid )

  End subroutine output
  
End module EM_1D

Program MT1Ddriver
  use EM_1D
  use lapack95
  implicit none
  Integer :: i, j, k
  real(kind=8) :: temp_h
  
  allocate( z(n+1), zc(n) )
  allocate( G(n,n+1), Av(n,n+1), Linv(n,n), Mmu(n,n), Msig(n+1,n+1), M(n+1,n+1) )
  call get_para
  
  Do j = 1, size(omiga)
    write( *,'(1x,a,g0,a)' )  '���ڼ���� ', j, ' ��Ƶ��!'
    !.. set up a deep enough mesh
    L = 3.d0 * skin(j)
    
    temp_h = (L-0.d0) / dble(n-1)
    h(1:n-1,1) = temp_h
    h(n,1) = h(n-1,1)*1.1d0
    
    sig(:,1) = sig0
    
    !.. mesh
    z(1) = 0.d0
    do i = 2, n+1
      z(i) = z(i-1) + h(i-1,1)
    end do
    
    do i = 1, n
      zc(i) = ( z(i)+z(i+1) ) / 2.d0
    end do
    
    !.. the linear system
    call ddx
    call ave
    call sdiag(1)
    call sdiag(2)
    call sdiag(3)
    call sdiag(4)
    
    !.. set the matrix system
    !.. [i*w           d/dz] [b] - s1  where s1 = 0
    !.. [1/mu  d/dz  -sigma] [e] - s2  where s2 = 0
    allocate( A(2*n,2*n), AA(2*n+1,2*n+1) )
    allocate( bb(2*n,1) )
    AA = 0.d0
    !.. �ȹ������ϽǴ�СΪn*n���Ӿ���
    forall( i = 1:n, k = 1:n, i == k ) AA(i,k) = cmplx( 0.d0, omiga(j) )
    !.. �������ϽǴ�СΪn*(n+1)���Ӿ���
    AA(1:n,n+1:2*n+1) = cmplx( matmul( Linv, G ), 0.d0 )
    !.. �������½Ǵ�СΪ(n+1)*n���Ӿ���
    AA(n+1:2*n+1,1:n) = cmplx( -1.d0 * matmul( matmul( transpose(G),transpose(Linv) ), Mmu ), 0.d0 )
    !.. �������½Ǵ�СΪ(n+1)*(n+1)���Ӿ���
    AA(n+1:2*n+1,n+1:2*n+1) = cmplx( -1.d0*Msig, 0.d0 )

    !.. ���任
    A(:,:) = AA(2:,2:)
    bb(:,1) = -1.d0 * AA(2:,1) 
    
    call gesv( A, bb )  !.. solve
    
    !if ( j == size(omiga) ) print*, bb
    b(1) = 1.d0  !.. �߽�����
    b(2:) = bb(1:n-1,1)
    e(:) = bb(n+1:,1)
    
    !.. extract data
    d(j)  = e(1)
    if ( j == size(omiga) ) call output
    deallocate( A, AA, bb )
  End do
  deallocate( z, zc )
  deallocate( G, Av, Linv, Mmu, Msig, M )
  
End program MT1Ddriver