!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

PROGRAM ZDF
  USE parameters_pt, ONLY : READ_FRAME, REMOVE_EQUIL, &
                            nequil, nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL ANALYSIS(frame)
  END DO

  CLOSE(1)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters_pt, ONLY : natoms, nframes, nequil, pos, &
                            box, ind_atom, Find_Ti5c_and_O2c
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(ind_atom(natoms))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
  END DO
 
  CLOSE(1)

  CALL Find_Ti5c_and_O2c(pos_file)

  OPEN(unit  =  1,file  =  pos_file)
  print *, "# step path_O path_H"
  !print *, "# step dO2c_H1 dOw_H1 dOw_H2 path_O path_H"
  !print *, "# step dO2c_H1 dOw_H1 Angle(O2c_H2_Ow) dO2c_H1(3)"
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS(step)
  USE parameters_pt, ONLY : index_O2c, Count_hydroxyl, &
                            min_known_path_cvs, nTi5c, &
                            BiasPotential
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: step
  INTEGER                    :: iat, iend(1)
  INTEGER                    :: n_hydroxyl
  REAL*8                     :: path_O, path_H
  REAL*8                     :: TotalBias 

  PRINT *, "####", step

  n_hydroxyl = Count_Hydroxyl(index_O2c)

  DO iat=1,nTi5c

    CALL min_known_path_cvs(index_O2c(iat),path_O,path_H,iend(1))

    IF (iend(1).gt.0) THEN

      TotalBias = BiasPotential(index_O2c(iat),1)
      TotalBias = BiasPotential(iend,1) + TotalBias

      WRITE(*,fmt="(i5,3(3X,F14.8),3X, i3)"), & 
        iend(1), path_O, path_H, TotalBias, n_hydroxyl
    END IF

  END DO

END SUBROUTINE ANALYSIS

