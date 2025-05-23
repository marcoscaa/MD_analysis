
MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: posm(:,:)
END MODULE histogram

PROGRAM VAC 
  USE parameters, ONLY : nframes, coarse
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_POS 
    CALL WRITE_VEL (frame)
  END DO

  CLOSE(1);CLOSE(2)

END PROGRAM VAC

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : posm
  USE parameters, ONLY : pos, vel, natoms, nframes
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, index_file)

  CALL READ_INDEX( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates at time t
  ALLOCATE(posm(natoms,3)) ! coordinates at time t-1
  ALLOCATE(vel(natoms,3)) ! velocities of water atoms 

  OPEN(unit = 1,file = file_pos)
  OPEN(unit = 2,file = 'out.vel')
  
END SUBROUTINE INITIALIZE

SUBROUTINE WRITE_VEL (step)
  USE histogram, ONLY : posm
  USE parameters, ONLY : pos, natoms, vel, dt, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: step
  INTEGER                    :: i,j

  IF (step==1) posm=pos
  WRITE(2,*) step

  DO i=1,natoms

    DO j=1,3
      vel(i,j) = pos(i,j) - posm(i,j)
      vel(i,j) = vel(i,j) - nint(vel(i,j)/box(j,j))*box(j,j)
      vel(i,j) = vel(i,j) / dt
    END DO

    WRITE(2, fmt='(3(E14.6,3X))') vel(i,1), vel(i,2), vel(i,3)

  END DO

  posm=pos

END SUBROUTINE WRITE_VEL
