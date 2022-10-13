module linpack

use mpi_f08

implicit none

integer, parameter  :: sp = selected_real_kind(6,37)      !single precision
integer, parameter  :: dp = selected_real_kind(15,307)    !double precisio

contains

subroutine linpack_real_sp(n)
  integer, intent(in)   :: n
  !local variables
  integer, parameter    :: niter = 100
  integer               :: iter
  real(dp)              :: tstart, tend
  real(sp), allocatable :: a(:,:)
  real(sp), allocatable :: b(:,:)
  real(sp), allocatable :: c(:,:)

  allocate(a(n,n), b(n,n), c(n,n))
  call random_seed()
  call random_number(a)
  call random_number(b)
  
  tstart = MPI_Wtime()
  do iter = 1, niter
    call sgemm("n","n",n,n,n,1.0_sp,a,n,b,n,0.0_sp,c,n)
  end do
  tend = MPI_Wtime()

  write(*,*)
  write(*,*)           "**************************************************"
  write(*,"(a,i6)")    " linpack for real(sp) n = ", n
  write(*,"(a,f12.2)") " Elasped time in s      = ", tend-tstart
  write(*,"(a,f8.2)")  " Perfromance in GFLOPS  = ", 2.0_dp*n**3/(tend-tstart)
  write(*,*)           "**************************************************"
  write(*,*)
end subroutine linpack_real_sp

subroutine linpack_real_dp(n)
  integer, intent(in)   :: n
  !local variables
  integer, parameter    :: niter = 100
  integer               :: iter
  real(dp)              :: tstart, tend
  real(dp), allocatable :: a(:,:)
  real(dp), allocatable :: b(:,:)
  real(dp), allocatable :: c(:,:)


  allocate(a(n,n), b(n,n), c(n,n))
  call random_seed()
  call random_number(a)
  call random_number(b)
  
  tstart = MPI_Wtime()
  do iter = 1, niter
    call dgemm("n","n",n,n,n,1.0_dp,a,n,b,n,0.0_dp,c,n)
  end do
  tend = MPI_Wtime()

  write(*,*)
  write(*,*)           "**************************************************"
  write(*,"(a,i6)")    " linpack for real(sp) n = ", n
  write(*,"(a,f12.2)") " Elasped time in s      = ", tend-tstart
  write(*,"(a,f8.2)")  " Perfromance in GFLOPS  = ", 2.0_dp*n**3/(tend-tstart)
  write(*,*)           "**************************************************"
  write(*,*)
end subroutine linpack_real_dp

subroutine linpack_cmplx_sp(n)
  integer, intent(in)      :: n
  !local variables
  integer, parameter       :: niter = 100
  integer                  :: iter
  real(dp)                 :: tstart, tend
  complex(sp), allocatable :: a(:,:)
  complex(sp), allocatable :: b(:,:)
  complex(sp), allocatable :: c(:,:)
  real(sp), allocatable    :: a_r(:,:), a_i(:,:)
  real(sp), allocatable    :: b_r(:,:), b_i(:,:)
  real(sp), allocatable    :: c_r(:,:), c_i(:,:)

  allocate(a(n,n), b(n,n), c(n,n))
  allocate(a_r(n,n), a_i(n,n), b_r(n,n), b_i(n,n), c_r(n,n), c_i(n,n))
  call random_seed()
  call random_number(a_r)
  call random_number(a_i)
  call random_number(b_r)
  call random_number(b_i)
  a = cmplx(a_r, a_i, kind=sp)
  b = cmplx(b_r, b_i, kind=sp)

  tstart = MPI_Wtime()
  do iter = 1, niter
    call cgemm("n","n",n,n,n,(1.0_sp,0.0_sp),a,n,b,n,(0.0_sp,0.0_sp),c,n)
  end do
  tend = MPI_Wtime()

  write(*,*)
  write(*,*)           "**************************************************"
  write(*,"(a,i6)")    " linpack for real(sp) n = ", n
  write(*,"(a,f12.2)") " Elasped time in s      = ", tend-tstart
  write(*,"(a,f8.2)")  " Perfromance in GFLOPS  = ", 6.0_dp*n**3/(tend-tstart)
  write(*,*)           "**************************************************"
  write(*,*)
end subroutine linpack_cmplx_sp

subroutine linpack_cmplx_dp(n)
  integer, intent(in)      :: n
  !local variables
  integer, parameter       :: niter = 100
  integer                  :: iter
  real(dp)                 :: tstart, tend
  complex(dp), allocatable :: a(:,:)
  complex(dp), allocatable :: b(:,:)
  complex(dp), allocatable :: c(:,:)
  real(dp), allocatable    :: a_r(:,:), a_i(:,:)
  real(dp), allocatable    :: b_r(:,:), b_i(:,:)
  real(dp), allocatable    :: c_r(:,:), c_i(:,:)

  allocate(a(n,n), b(n,n), c(n,n))
  allocate(a_r(n,n), a_i(n,n), b_r(n,n), b_i(n,n), c_r(n,n), c_i(n,n))
  call random_seed()
  call random_number(a_r)
  call random_number(a_i)
  call random_number(b_r)
  call random_number(b_i)
  a = cmplx(a_r, a_i, kind=dp)
  b = cmplx(b_r, b_i, kind=dp)

  tstart = MPI_Wtime()
  do iter = 1, niter
    call zgemm("n","n",n,n,n,(1.0_dp,0.0_dp),a,n,b,n,(0.0_dp,0.0_dp),c,n)
  end do
  tend = MPI_Wtime()

  write(*,*)
  write(*,*)           "**************************************************"
  write(*,"(a,i6)")    " linpack for real(sp) n = ", n
  write(*,"(a,f12.2)") " Elasped time in s      = ", tend-tstart
  write(*,"(a,f8.2)")  " Perfromance in GFLOPS  = ", 6.0_dp*n**3/(tend-tstart)
  write(*,*)           "**************************************************"
  write(*,*)
end subroutine linpack_cmplx_dp

end module linpack