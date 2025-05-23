
PROGRAM ZDF
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  
  CALL MAIN 

  CLOSE(1)
  CLOSE(2)
  CLOSE(3)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_ice, file_water

  CALL getarg(1, file_water)
  CALL getarg(2, file_ice)


  OPEN(unit  =  1,file  =  file_water)
  OPEN(unit  =  2,file  =  file_ice)
  OPEN(unit  =  3,file  =  'ice-water.gro')
 
END SUBROUTINE INITIALIZE

SUBROUTINE MAIN 
  !Subroutine reads gro files
  IMPLICIT NONE
  INTEGER                    :: i, j, junk1, iat
  INTEGER                    :: nat, natw, nati, imol
  CHARACTER(5),ALLOCATABLE   :: moltype(:),atomtype(:)
  REAL*8, ALLOCATABLE        :: pos(:,:), vel(:,:)
  REAL*8                     :: box(6)

  !Discarding first line
  READ(1,*) 
  READ(1,*) natw
  READ(2,*) 
  READ(2,*) nati
  nat = natw+nati

  ALLOCATE(pos(3,nat)); pos=0
  ALLOCATE(vel(3,nat)); vel=0
  ALLOCATE(moltype(nat))
  ALLOCATE(atomtype(nat))

  !Read water gro file
  DO i = 1, natw 
    READ(1, fmt="(i5,2a5,i5,6f8.3)"), junk1, moltype(i), atomtype(i), junk1, pos(:,i), vel(:,i)
  END DO

  READ(1,*), box(1), box(2), box(3)

  !Read ice gro file
  DO i = natw+1, nat 
    READ(2, fmt="(i5,2a5,i5,6f8.3)"), junk1, moltype(i), atomtype(i), junk1, pos(:,i), vel(:,i)
  END DO

  READ(2,*), box(4), box(5), box(6)

  !Start writing
  WRITE(3, fmt='(a5)'), "Interface"
  WRITE(3, fmt='(i5)'), 2*natw+nati

  iat = 1
  DO j = 1,1
    pos(3,:) = pos(3,:) + box(3)*(j-1)
    DO i = 1,natw
      imol = int( (iat-1)/4 ) + 1
      WRITE(3, fmt="(i5,2a5,i5,6f8.3)"), imol, moltype(i), atomtype(i), iat, pos(:,i), vel(:,i)
      iat = iat + 1
    END DO
  END DO

  pos(3,:) = pos(3,:) + box(3) + 0.1
  DO i = natw+1,nat
    imol = int( (iat-1)/4 ) + 1
    WRITE(3, fmt="(i5,2a5,i5,6f8.3)"), imol, moltype(i), atomtype(i), iat, pos(:,i), vel(:,i)
    iat = iat + 1
  END DO

  pos(3,:) = pos(3,:) + box(6) + 0.1
  DO j = 1,1
    pos(3,:) = pos(3,:) + box(3)*(j-1)
    DO i = 1,natw  
      imol = int( (iat-1)/4 ) + 1
      WRITE(3, fmt="(i5,2a5,i5,6f8.3)"), imol, moltype(i), atomtype(i), iat, pos(:,i), vel(:,i)
      iat = iat + 1
    END DO
  END DO

  WRITE(3, fmt='(3f10.5)'), box(4), box(5), 2*box(3) + box(6) + 0.2

END SUBROUTINE
