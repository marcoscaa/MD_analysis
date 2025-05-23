module myxyz
  CHARACTER(5), ALLOCATABLE :: type_map(:)
end module myxyz

PROGRAM XYZ
  IMPLICIT NONE
  INTEGER :: iostatus

  CALL INITIALIZE
  
  DO WHILE (iostatus==0) 
    CALL READ_ATOM_REDUCED_IO(iostatus)
    CALL PRINT_ATOMS 
  END DO

END PROGRAM XYZ

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, maxr, indexC
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file

  !User should provide filename as argument 
  CALL getarg(1, pos_file)
  READ(*,*) indexC, maxr

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = 'out.xyz')

END SUBROUTINE INITIALIZE

SUBROUTINE PRINT_ATOMS
  USE parameters, ONLY : natoms, pos, maxr, indexC, atype
  IMPLICIT NONE
  INTEGER                     :: iat
  INTEGER                     :: ind
  INTEGER                     :: n, new_type(natoms)
  REAL*8                      :: d(4), new_pos(3,natoms)

  new_pos=0
  new_type=0
  n = 1 !total number of atoms 
  DO iat = 1, natoms
    IF (iat==indexC) THEN
      new_type(1)=atype(iat)
    ELSE
      CALL VECTOR_DISTANCEi(iat,indexC,d)
      IF (d(4)<maxr) THEN
        n = n+1
        new_pos(:,n)=d(1:3)
        new_type(n)=atype(iat)
      END IF
    END IF
  END DO

  CALL WRITE_XYZ(new_pos,new_type,n)

END SUBROUTINE PRINT_ATOMS

SUBROUTINE WRITE_XYZ(pos,atype,n)
  use myxyz, only : type_map
  IMPLICIT NONE
  INTEGER, INTENT(IN)      :: n, atype(n)
  REAL*8, INTENT(IN)       :: pos(3,n)
  INTEGER                  :: iat

  IF (.not.allocated(type_map)) THEN
    OPEN(10,file='type_map')
    allocate(type_map(maxval(atype)))
    READ(10,*) type_map
    CLOSE(10)
  END IF

  WRITE(2,*) n
  WRITE(2,*)

  DO iat=1,n
    WRITE(2,fmt='(A5,3(3X,F12.8))') type_map(atype(iat)), pos(:,iat)
  END DO

END SUBROUTINE WRITE_XYZ
