program linpack_test

  use mpi_f08
  use linpack

  implicit none

  integer :: n, niter
  integer :: ierr


  call MPI_Init(ierr)

  n = 100
  niter = 1000
  call linpack_real_sp(n, niter)
  call linpack_real_dp(n, niter)
  call linpack_cmplx_sp(n, niter)
  call linpack_cmplx_dp(n, niter)

  n = 500
  niter = 100
  call linpack_real_sp(n, niter)
  call linpack_real_dp(n, niter)
  call linpack_cmplx_sp(n, niter)
  call linpack_cmplx_dp(n, niter)

  n = 2000
  niter = 10
  call linpack_real_sp(n, niter)
  call linpack_real_dp(n, niter)
  call linpack_cmplx_sp(n, niter)
  call linpack_cmplx_dp(n, niter)

  call MPI_Finalize(ierr)

end program linpack_test