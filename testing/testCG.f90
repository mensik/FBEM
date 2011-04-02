program testCG

  use la
  use timer
  use files

  double precision, allocatable, dimension(:,:) :: A
  double precision, allocatable, dimension(:) :: x,b,temp

  integer :: m, n, i, j

  call get_matrix_size(m,n)

  allocate(A(m,n), x(n), b(m),temp(m))

  call load_matrix(A)
  call load_vector(b)

  call start_timer()
  call cg(A,b,x)
  call update_timer()

  write (*,*) 'Time : ', time
end program testCG
