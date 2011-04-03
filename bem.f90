module bem

    use la

contains
!------------------------------------------------------------------------------
  subroutine assembleV(i_start, i_end, j_start, j_end, nodes, edges, V)
!------------------------------------------------------------------------------ 
    implicit none

    integer, intent(in) :: i_start, i_end, j_start, j_end
    real, intent(in) :: nodes(:,:)
    integer, intent(in) :: edges(:,:)
    real, intent(out) :: V(:,:)

    integer :: i, j
    real :: a(2),b(2),x(2)

    do i = i_start, i_end
      do j = j_start, j_end
        x = 0.5 * (nodes(edges(i,1),:) + nodes(edges(i,2),:))
        a = nodes(edges(j,1),:)
        b = nodes(edges(j,2),:)
        V(i,j) = V_loc(a, b, x, i == j)
      end do
    end do

  end subroutine assembleV
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  real function matA(a,b,x)
!------------------------------------------------------------------------------
    implicit none
    
    real, intent(in) :: a(2), b(2), x(2)
    real :: matA

    matA = atan(in_prod(b - a, x - a) / (cross_prod(a, b - x) + cross_prod(b, x))) &
         - atan(in_prod(b - a, x - b) / (cross_prod(a, b - x) + cross_prod(b, x)))

  end function matA
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  real function V_loc(a, b, x, same_element)
!------------------------------------------------------------------------------ 
    implicit none
    
    real, intent(in) :: x(2), a(2), b(2)
    logical, intent(in) :: same_element
    real :: V_loc
    real :: el_length
    
    el_length = norm_2(b - a)

    V_loc = 0

    if (same_element) then
      V_loc = el_length * (log(el_length) - log(2.0) - 1.0)
    else
      V_loc = 1 / el_length * ((cross_prod(x, a - b) + cross_prod(a,b)) * matA(a, b, x) &
           + in_prod(b - a, (x - a) * log(norm_2(x - a)) - (x - b) * log(norm_2(x - b))))
    end if
  
  end function V_loc
!------------------------------------------------------------------------------
end module bem
