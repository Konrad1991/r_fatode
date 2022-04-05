module lapack
      implicit none
      save
      integer :: nvar
      double precision,allocatable :: fjac(:,:),e1(:,:)
      integer, allocatable :: ip1(:), ip2(:)
      complex(kind=selected_real_kind(14,300)),allocatable :: e2(:,:)
contains
      subroutine lapack_decomp(ising,hgamma)
      integer :: ising,i,j
      double precision :: hgamma
      do j=1,nvar
        do i=1,nvar
          e1(i,j) = -fjac(i,j)
        end do
        e1(j,j) = e1(j,j) +hgamma
      end do
      call dgetrf( nvar, nvar, e1, nvar, ip1, ising )
      end subroutine lapack_decomp

      subroutine lapack_decomp_cmp(ising,alpha,beta)
      integer :: ising,i,j
      double precision :: alpha,beta
      do j=1,nvar
         do i=1,nvar
            e2(i,j) = dcmplx( -fjac(i,j), 0d0)
         end do
         e2(j,j) =e2(j,j) + cmplx( alpha, beta )
      end do
      call zgetrf( nvar, nvar, e2, nvar, ip2, ising )
      end subroutine lapack_decomp_cmp

      subroutine lapack_solve(rhs)
      double precision ::rhs(nvar)
      integer :: ising
      call dgetrs( 'n', nvar, 1, e1, nvar, ip1, rhs, nvar, ising)
      end subroutine lapack_solve

      subroutine lapack_solve_cmp(b,bz)
      double precision :: b(nvar), bz(nvar)
      complex(kind=selected_real_kind(14,300)) :: rhs(nvar)
      integer :: ising,i
      do i=1,nvar
         rhs(i) = dcmplx(b(i),bz(i))
      end do
      call zgetrs( 'n', nvar, 1, e2, nvar, ip2, rhs, nvar, ising)
      do i = 1,nvar
        b(i) = dble(rhs(i))
        bz(i) = aimag(rhs(i))
      end do
      end subroutine lapack_solve_cmp


      subroutine lapack_init(n)
      integer :: n,state
      nvar = n
      allocate(ip1(nvar),ip2(nvar),fjac(nvar,nvar),e1(nvar,nvar),&
                                       e2(nvar,nvar),STAT = state)
      if(state .ne. 0)stop 'Allocation error in lapack_init'
      end subroutine lapack_init

      subroutine lapack_free
      integer :: state
      deallocate(ip1,ip2,fjac,e1,e2,STAT=state)
      if(state .ne.0) stop 'Deallocation error in lapack_free'
      end subroutine lapack_free

end module lapack



module ls_solver
      use lapack
      implicit none
      public
contains
      subroutine lss_jac(t,y,jac)
      double precision :: t,y(nvar)
      external :: jac
      call jac(nvar,t,y,fjac)
      end subroutine lss_jac

      subroutine lss_decomp(hgamma, ising)
      implicit none
      integer :: ising
      double precision :: hgamma
      call lapack_decomp(ising,hgamma)
      end subroutine lss_decomp

      subroutine lss_decomp_cmp(alpha, beta, ising)
      integer :: ising
      double precision :: alpha,beta

      call lapack_decomp_cmp(ising,alpha,beta)
      end subroutine lss_decomp_cmp

      subroutine lss_solve(rhs)
      double precision ::rhs(nvar)

      call lapack_solve(rhs)

!      print *, rhs
      end subroutine lss_solve

      subroutine lss_solve_cmp(b,bz)
      double precision  :: b(nvar), bz(nvar)

      call lapack_solve_cmp(b,bz)

!      print *, rhs
      end subroutine lss_solve_cmp

      subroutine lss_init(n,nn)
      integer :: n,nn

      call lapack_init(n)

      end subroutine lss_init

      subroutine lss_free

      call lapack_free

      end subroutine

end module ls_solver
