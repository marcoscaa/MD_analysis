!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE parameters
  IMPLICIT NONE 
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE              :: pos, vel
  REAL*8, DIMENSION(3)                               :: box
  INTEGER                                            :: natoms, nmol 

END MODULE parameters

PROGRAM ZDF
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  
  CALL READ_FRAME
  CALL PRINT_GRO

  CLOSE(1)
  CLOSE(2)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  OPEN(unit  =  1,file  =  file_name)
  OPEN(unit  =  2,file  =  'out.gro')
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME 
  !Subroutine reads gro files
  USE parameters, ONLY : pos, vel, box, nmol, natoms
  IMPLICIT NONE
  INTEGER                    :: i, j, junk1
  REAL*8                     :: vtmp(3), ptmp(3) 
  CHARACTER(5)               :: moltype,atomtype
  LOGICAL                    :: in_slab

  nmol = 0
  !Discarding first line
  READ(1,*) 
  READ(1,*) natoms

  ALLOCATE(pos(3,4,natoms/4)); pos=0
  ALLOCATE(vel(3,4,natoms/4)); vel=0

  i = 1
  DO WHILE ( i <= natoms/4 )

    READ(1, fmt="(i5,2a5,i5,6f8.3)"), junk1, moltype, atomtype, junk1, ptmp(:), vtmp(:)

    if ( ( moltype.eq.'SOL  ' ) .and. ( in_slab(ptmp(3)) ) ) THEN
    
      DO j=1,3
        READ(1, *)
      END DO

    else

      nmol = nmol + 1
      !Oxygem coordinates
      pos(:,1,nmol) = ptmp(:)
      vel(:,1,nmol) = vtmp(:)
      !Hydrogen coordinates
      READ(1, fmt="(i5,2a5,i5,6f8.3)"), junk1,moltype,atomtype,junk1,pos(:,2,nmol),vel(:,2,nmol)
      READ(1, fmt="(i5,2a5,i5,6f8.3)"), junk1,moltype,atomtype,junk1,pos(:,3,nmol),vel(:,3,nmol)
      !Virtual site                                                               
      READ(1, fmt="(i5,2a5,i5,6f8.3)"), junk1,moltype,atomtype,junk1,pos(:,4,nmol),vel(:,4,nmol)


    end if

    i = i + 1

  END DO

  natoms = nmol*4

  READ(1,*), box(1), box(2), box(3)

END SUBROUTINE

SUBROUTINE PRINT_GRO
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: imol,iat

  WRITE(2, fmt='(a5)'), "Cleaned"
  WRITE(2, fmt='(i5)'), natoms

  iat = 1
  DO imol = 1, nmol
    WRITE(2, fmt='(i5,2a5,i5,6f8.3)'), imol, "SOL  ", "   OW", iat  , pos(:,1,imol), vel(:,1,imol)
    WRITE(2, fmt='(i5,2a5,i5,6f8.3)'), imol, "SOL  ", "  HW1", iat+1, pos(:,2,imol), vel(:,2,imol)
    WRITE(2, fmt='(i5,2a5,i5,6f8.3)'), imol, "SOL  ", "  HW2", iat+2, pos(:,3,imol), vel(:,3,imol)
    WRITE(2, fmt='(i5,2a5,i5,6f8.3)'), imol, "SOL  ", "   MW", iat+3, pos(:,4,imol), vel(:,4,imol)
    iat = iat + 4
  END DO

  WRITE(2, fmt='(3f10.5)'), box(1), box(2), box(3) 

END SUBROUTINE

LOGICAL FUNCTION in_slab ( z )
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: z
  REAL*8, PARAMETER          :: up_z=8.405, dwn_z=4.058

  IF ( ( z >= dwn_z ) .and. ( z <= up_z ) ) THEN
    in_slab = .true.
  ELSE
    in_slab = .false.
  END IF

END FUNCTION in_slab
