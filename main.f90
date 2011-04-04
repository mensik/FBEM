program fbem

  use discretization
  use bem
  use files

  real , parameter :: h = 3
  real, parameter :: er = 5
  real, allocatable :: nodes(:,:)
  real, allocatable :: A(:,:), temp(:,:)
  integer, allocatable :: edges(:,:)
  integer :: b(3,2)
  integer :: node_count
  integer :: v_size, r_size,  i, j, bl_s(4), bl_e(4)

  node_count = nodeCount(2.0,6.0,2.0,2.0,h)

  write (*,*) node_count

  allocate(nodes(node_count, 2), edges(node_count,2))

  call make_borders(2.0,6.0,2.0,2.0,2.0,h, nodes, edges, b)

  v_size = b(1,2) - b(1,1) + 1
  r_size = b(2,2) - b(2,1) + 1

  bl_s = (/ 1, r_size + 1, 2*r_size + 1, 2*r_size + v_size + 1/)
  bl_e = (/ r_size, 2*r_size, 2*r_size + v_size, 2*r_size + 2*v_size/)

  allocate(A((r_size + v_size)* 2, (r_size + v_size) * 2))
  
  A = 0

! 1st line of matrix A
!------------------------------------------------------------------------------
  call assembleL(b(2,1), b(2,2), b(2,1), b(2,2), nodes, edges, A(bl_s(1):bl_e(1), bl_s(1):bl_e(1)))
  A(bl_s(1):bl_e(1) , bl_s(2):bl_e(2)) =  A(bl_s(1):bl_e(1), bl_s(1):bl_e(1))
  call assembleL(b(2,1), b(2,2), b(3,1), b(3,2), nodes, edges, A(bl_s(1):bl_e(1), bl_s(3):bl_e(3)))
  call assembleL(b(2,1), b(2,2), b(1,1), b(1,2), nodes, edges, A(bl_s(1):bl_e(1), bl_s(4):bl_e(4)))

  do i = 1,r_size
    A(i,i) = A(i,i) + 0.5
    A(i, i + r_size) = A(i, i + r_size) - 0.5 
    
    do j = 1, r_size
      A(i, j) = -er * A(i, j)
    end do
  end do
!------------------------------------------------------------------------------
! 2nd line of matrix A
!------------------------------------------------------------------------------
  call assembleV(b(2,1), b(2,2), b(2,1), b(2,2), nodes, edges, A(bl_s(2):bl_e(2), bl_s(1):bl_e(1)))
  A(bl_s(2):bl_e(2), bl_s(1):bl_e(1)) = -A(bl_s(2):bl_e(2), bl_s(1):bl_e(1))
  A(bl_s(2):bl_e(2), bl_s(2):bl_e(2)) = -A(bl_s(2):bl_e(2), bl_s(1):bl_e(1))
  call assembleV(b(2,1), b(2,2), b(3,1), b(3,2), nodes, edges, A(bl_s(2):bl_e(2), bl_s(3):bl_e(3)))
  call assembleV(b(2,1), b(2,2), b(1,1), b(1,2), nodes, edges, A(bl_s(2):bl_e(2), bl_s(4):bl_e(4)))
!------------------------------------------------------------------------------
! 3rd line fo matrix A
!------------------------------------------------------------------------------
  call assembleV(b(3,1), b(3,2), b(2,1), b(2,2), nodes, edges, A(bl_s(3):bl_e(3), bl_s(2):bl_e(2)))
  call assembleV(b(3,1), b(3,2), b(3,1), b(3,2), nodes, edges, A(bl_s(3):bl_e(3), bl_s(3):bl_e(3)))
  call assembleV(b(3,1), b(3,2), b(1,1), b(1,2), nodes, edges, A(bl_s(3):bl_e(3), bl_s(4):bl_e(4)))
!------------------------------------------------------------------------------
! 4th line of matrix A
!------------------------------------------------------------------------------
  call assembleV(b(1,1), b(1,2), b(2,1), b(2,2), nodes, edges, A(bl_s(4):bl_e(4), bl_s(2):bl_e(2)))
  call assembleV(b(1,1), b(1,2), b(3,1), b(3,2), nodes, edges, A(bl_s(4):bl_e(4), bl_s(3):bl_e(3)))
  call assembleV(b(1,1), b(1,2), b(1,1), b(1,2), nodes, edges, A(bl_s(4):bl_e(4), bl_s(4):bl_e(4)))
!------------------------------------------------------------------------------

  call save_matrix(A)
  
  deallocate(nodes, edges, A)
end program fbem
