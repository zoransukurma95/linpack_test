program linpack_test

  use mpi_f08
  use linpack

  implicit none

  integer :: n, niter
  integer :: rank, ierr


  call MPI_Init(ierr)
  call MPI_Comm_Rank(mpi_comm_world, rank)

  n = 500
  niter = 200
  call linpack_real_sp(n, niter, rank)
  call linpack_real_dp(n, niter, rank)
  call linpack_cmplx_sp(n, niter, rank)
  call linpack_cmplx_dp(n, niter, rank)
  if (rank == 0) then
    write(*,*)           "**************************************************"
    write(*,*)
    write(*,*)
    write(*,*)
  end if

  n = 1000
  niter = 100
  call linpack_real_sp(n, niter, rank)
  call linpack_real_dp(n, niter, rank)
  call linpack_cmplx_sp(n, niter, rank)
  call linpack_cmplx_dp(n, niter, rank)
  if (rank == 0) then
    write(*,*)           "**************************************************"
    write(*,*)
    write(*,*)
    write(*,*)
  end if

  n = 2000
  niter = 10
  call linpack_real_sp(n, niter, rank)
  call linpack_real_dp(n, niter, rank)
  call linpack_cmplx_sp(n, niter, rank)
  call linpack_cmplx_dp(n, niter, rank)
  if (rank == 0) then
    write(*,*)           "**************************************************"
    write(*,*)
    write(*,*)
    write(*,*)
  end if

  n = 3500
  niter = 5
  call linpack_real_sp(n, niter, rank)
  call linpack_real_dp(n, niter, rank)
  call linpack_cmplx_sp(n, niter, rank)
  call linpack_cmplx_dp(n, niter, rank)
  if (rank == 0) then
    write(*,*)           "**************************************************"
  end if

  call MPI_Finalize(ierr)

end program linpack_test