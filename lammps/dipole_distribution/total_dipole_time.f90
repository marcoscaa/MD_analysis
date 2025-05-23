MODULE wannier_mod 
  IMPLICIT NONE 
  INTEGER                      :: atypeO, atypeH, atypeC
  REAL*8                       :: T, center(3), dipole(3)
  REAL*8, PARAMETER            :: convert_to_debye=1.d0/0.2081943d0
END MODULE wannier_mod

PROGRAM Dipole_Distribution
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL READ_WANNIER
    CALL COMPUTE_WATER_DIPOLE
    CALL PRINT_RESULTS(frame)
  END DO

END PROGRAM Dipole_Distribution

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, dt, &
                         atype, stride, nhist, nwann, wannier
  USE wannier_mod, ONLY : dipole,atypeO, atypeH, atypeC, T
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, wfc_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, wfc_file)
  CALL getarg(3, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms,nwann, nframes, T, dt
  READ(1,*) atypeO, atypeH, atypeC
  !READ(1, *) center 
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(wannier(3,natoms)); wannier=0.d0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
  OPEN(unit = 3, file =  wfc_file)
  OPEN(unit = 2, file = "Total_dipole_vs_time.dat")
  WRITE (2,fmt="(A36)") '# Time(ps) Total_Water_Dipole(X,Y,Z)' 
 
END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_WATER_DIPOLE
  USE parameters, ONLY : pos, atype, natoms, wannier
  USE wannier_mod, ONLY : dipole, convert_to_debye, atypeO, &
                          atypeH
  IMPLICIT NONE
  INTEGER                    :: i, iH, ind_H(2)
  REAL*8                     :: d(4), ionic_dipole(3), electric_dipole(3)

  dipole=0.d0
  ionic_dipole=0.d0
  electric_dipole=0.d0

  DO i = 1,natoms

    !Selecting only Ow
    IF ( atype(i) == atypeO )  THEN

      CALL get_H(i, ind_H, atypeH)
      DO iH=1,2
        CALL Distance_Vector_IND(ind_H(iH),i,d)
        ionic_dipole=ionic_dipole+d(1:3)
      END DO

      electric_dipole = electric_dipole -8*wannier(:,i)

    END IF 

  END DO !i 

  dipole=electric_dipole!+ionic_dipole
  dipole=dipole*convert_to_debye

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

SUBROUTINE PRINT_RESULTS(step)
  USE parameters, ONLY : dt
  USE wannier_mod, ONLY : dipole
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: step

  WRITE (2,fmt="(F12.3,3(3X,E15.7))") dt*(float(step-1)), dipole 

END SUBROUTINE PRINT_RESULTS
