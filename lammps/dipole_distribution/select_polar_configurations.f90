MODULE wannier_mod 
  IMPLICIT NONE 
  INTEGER                      :: atypeO, atypeH, atypeC
  REAL*8, ALLOCATABLE          :: dipole(:,:) 
  REAL*8, ALLOCATABLE          :: ionic_dipole(:,:), electric_dipole(:,:)
  REAL*8                       :: T, dipole_boundary
  REAL*8, PARAMETER            :: convert_to_debye=1.d0/0.2081943d0
END MODULE wannier_mod

PROGRAM Dipole_Distribution
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL READ_WANNIER
    CALL COMPUTE_WATER_DIPOLE
    CALL PRINT_RESULTS
  END DO

END PROGRAM Dipole_Distribution

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist, nwann, wannier
  USE wannier_mod, ONLY : dipole,atypeO, atypeH, atypeC, T, &
                          ionic_dipole, electric_dipole, dipole_boundary
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, wfc_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, wfc_file)
  CALL getarg(3, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms,nwann, nframes, T!, nhist, rmax
  READ(1,*) atypeO, atypeH, atypeC
  READ(1,*) dipole_boundary
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(dipole(3,natoms)); dipole=0.d0
  ALLOCATE(ionic_dipole(3,natoms)); ionic_dipole=0.d0
  ALLOCATE(electric_dipole(3,natoms)); electric_dipole=0.d0
  ALLOCATE(wannier(3,natoms)); wannier=0.d0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
  OPEN(unit = 2, file = 'configurations_dipole_selected.xyz')
  OPEN(unit = 3, file =  wfc_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_WATER_DIPOLE
  USE parameters, ONLY : pos, atype, natoms, wannier
  USE wannier_mod, ONLY : dipole, convert_to_debye, atypeO, &
                          atypeH, electric_dipole, ionic_dipole
  IMPLICIT NONE
  INTEGER                    :: i, iH, ind_H(2)
  REAL*8                     :: d(4)

  dipole=0.d0
  ionic_dipole=0.d0
  electric_dipole=0.d0

  DO i = 1,natoms

    !Selecting only Ow
    IF ( atype(i) == atypeO )  THEN

      CALL get_H(i, ind_H, atypeH)
      DO iH=1,2
        CALL Distance_Vector_IND(ind_H(iH),i,d)
        ionic_dipole(:,i)=ionic_dipole(:,i)+d(1:3)
      END DO

      electric_dipole(:,i) = electric_dipole(:,i) -8*wannier(:,i)

    END IF 

  END DO !i 

  electric_dipole=electric_dipole*convert_to_debye
  ionic_dipole=ionic_dipole*convert_to_debye
  dipole=electric_dipole+ionic_dipole

END SUBROUTINE COMPUTE_WATER_DIPOLE 

SUBROUTINE get_H(ind_O, ind_H,atypeH)
  USE parameters, ONLY : atype, pos, cutoff_OH, natoms
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: atypeH,ind_O
  INTEGER                    :: i, k, ind_H(2)
  REAL*8                     :: Dist

  i=1;k=1

  ind_H = 0
  DO WHILE (( i <= natoms) .and. ( k <= 2))

    !Should be the index for H
    IF ( atype(i) == atypeH ) THEN

      IF ( Dist(ind_O,i) < cutoff_OH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

  IF (k.ne.3) THEN
    print *, 'Found water ion with charge ', k-3
  END IF
END SUBROUTINE get_H

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : atype, natoms, pos, box
  USE wannier_mod, ONLY : dipole, dipole_boundary
  IMPLICIT NONE
  INTEGER :: iat
  REAL*8  :: total_dipole
  CHARACTER(2) :: Symbol

  total_dipole = sum(dipole(3,:)) / box(3)
  IF (total_dipole>dipole_boundary) THEN
     WRITE(2,*) natoms
     WRITE(2,fmt="(A9,9(F9.5,1X),A1,2X,A10,F12.8,A1)") 'Lattice="', &
        box(1),0.0,0.0,0.0, &
        box(2),0.0,0.0,0.0, &
        box(3),'"', &
        'dipole_z="',total_dipole,'"'

     DO iat=1,natoms
       WRITE(2,fmt="(A1,3(3X,F12.8))") Symbol(atype(iat)), pos(:,iat)
     END DO
  END IF

END SUBROUTINE PRINT_RESULTS

CHARACTER(2) FUNCTION Symbol(atomic_index)
  USE wannier_mod, only : atypeO,atypeH,atypeC 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: atomic_index
  IF (atomic_index==atypeC) THEN
      Symbol = 'C'
  ELSEIF (atomic_index==atypeO) THEN
      Symbol = 'O'
  ELSEIF (atomic_index==atypeH) THEN
      Symbol = 'H'
  ELSE
      Symbol = 'X'
  ENDIF
END FUNCTION Symbol
