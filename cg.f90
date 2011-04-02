!******************************************************************************
      module la
!******************************************************************************
  real, parameter :: PREC = 1e-10


contains
!******************************************************************************

!------------------------------------------------------------------------------
  subroutine cg(A,b,x)
!------------------------------------------------------------------------------

    implicit none

    real, dimension (:,:), intent(in) :: A
    real, dimension (:), intent(in) :: b
    real, dimension (:), intent(out) :: x
  
    real, dimension (:), allocatable :: r, p, Ap
    real :: alpha, r_sold, r_sold_new
    integer counter, i, n

!------------------------------------------------------------------------------
    n = size(b)
    allocate(r(n), p(n), Ap(n))

    x = b
    call mv_mult(A,x,Ap)
    r = b - Ap
    p = r
    r_sold = sum(r(:)**2)

    counter = 0
    do
      counter = counter + 1
      call mv_mult(A,p,Ap)
      alpha = r_sold / sum(p(:) * Ap(:))
      x = x + alpha * p
      r = r - alpha * Ap
      r_sold_new = sum(r(:)**2)

      if (r_sold_new < PREC) then
        exit
      end if
      
      p = r + (r_sold_new / r_sold) * p
      r_sold = r_sold_new
    end do
!------------------------------------------------------------------------------ 
  end subroutine cg
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine mv_mult(A,x,y)
!------------------------------------------------------------------------------
  
    implicit none

    real, dimension (:,:), intent(in) :: A
    real, dimension (:), intent(in) :: x
    real, dimension (:), intent(out) :: y
  
    real :: row_sum
    integer :: m , n, i, j

!------------------------------------------------------------------------------
    m = size(A,1)
    n = size(A,2)

    if (n /= size(x)) then
      write (*,*) 'ERROR: input vector has incorrect size'
      stop
    else if (m /= size(y)) then
      write (*,*) 'ERROR: output vector has incorrect size', PREC
      stop
    end if

    !$OMP PARALLEL DO SHARED (A,x,y) PRIVATE (i,j,row_sum)
    do i = 1, m
      row_sum = 0
      do j = 1, n
        row_sum = row_sum + A(i,j)*x(j)
      end do
      y(i) = row_sum
    end do
    !$OMP END PARALLEL DO
!------------------------------------------------------------------------------
  end subroutine mv_mult
!------------------------------------------------------------------------------

end module la
!******************************************************************************
