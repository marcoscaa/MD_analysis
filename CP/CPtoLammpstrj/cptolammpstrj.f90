!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE      :: pos
  DOUBLE PRECISION, DIMENSION(3,3)                   :: box
  INTEGER, ALLOCATABLE                               :: ind_atom(:)
  INTEGER                                            :: natoms, nframes, nequil

END MODULE parameters

PROGRAM ZDF
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL PRINT_lammpstrj(frame)
  END DO

  CLOSE(1)
  CLOSE(2)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, box_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, box_file)
  CALL getarg(3, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(ind_atom(natoms))

  !Read the index file
  READ(1,*), ind_atom  
 
  CLOSE(1)

  OPEN(unit  =  1,file  =  pos_file)
  OPEN(unit  =  2,file  =  box_file)
  OPEN(unit  =  3,file  =  'out.lammpstrj')
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i,j
  DOUBLE PRECISION           :: invbox(3,3)  

  READ(1,*)
  READ(2,*)

  DO i = 1,natoms
    READ(1, fmt='(3(3X,E22.14))'), pos(:,i)
  END DO

  DO i = 1,3
    READ(2, *), box(i,1), box(i,2), box(i,3)
  END DO

  CALL MATINV3(box,invbox) 

  DO i = 1,natoms
    pos(:,i) = MATMUL(invbox,pos(:,i))
    pos(:,i) = pos(:,i) - nint(pos(:,i)) + 0.5
  END DO

END SUBROUTINE

SUBROUTINE REMOVE_EQUIL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+1)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE PRINT_Lammpstrj(timestep)
  !Works only for orthogonal cells
  USE parameters
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: timestep
  INTEGER                    :: i

  WRITE(3, *) 'ITEM: TIMESTEP'
  WRITE(3, *), timestep
  WRITE(3, *) 'ITEM: NUMBER OF ATOMS'
  WRITE(3, *), natoms
  WRITE(3, *) 'ITEM: BOX BOUNDS pp pp pp'
  WRITE(3, fmt='(2(F12.8,2X))') 0.00, box(1,1)*0.529177
  WRITE(3, fmt='(2(F12.8,2X))') 0.00, box(2,2)*0.529177
  WRITE(3, fmt='(2(F12.8,2X))') 0.00, box(3,3)*0.529177
  WRITE(3, *) 'ITEM: ATOMS id type xs ys zs'

  DO i = 1, natoms
    WRITE(3, fmt='(I5, I3, 3(F12.8, 3X))') i, ind_atom(i), pos(:,i)
  END DO

END SUBROUTINE

SUBROUTINE MATINV3(A,B)
  IMPLICIT NONE
  !! Performs a direct calculation of the inverse of a 3×3 matrix.
  DOUBLE PRECISION, intent(in) :: A(3,3)   !! Matrix
  DOUBLE PRECISION             :: B(3,3)   !! Inverse matrix
  DOUBLE PRECISION             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

END SUBROUTINE
