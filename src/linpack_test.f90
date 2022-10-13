program linpack_test

  use mpi_f08
  use linpack

  implicit none

  integer :: n
  integer :: ierr


  call MPI_Init(ierr)

  n = 500

  call linpack_real_sp(n)
  call linpack_real_dp(n)
  call linpack_cmplx_sp(n)
  call linpack_cmplx_dp(n)

  call MPI_Finalize(ierr)

end program linpack_test