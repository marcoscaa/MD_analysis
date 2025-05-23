PROGRAM ProtonTransfer
  USE parameters_pt, ONLY : nequil, nframes, &
                            REMOVE_EQUIL, READ_FRAME
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL ANALYSIS(frame)
  END DO

  CLOSE(1)

END PROGRAM ProtonTransfer

SUBROUTINE INITIALIZE
  USE parameters_pt, ONLY : natoms, nframes, nequil, &
                            pos, ind_atom, &
                            Find_Ti5c_and_S2c
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

  CLOSE(1)

  CALL Find_Ti5c_and_S2c(pos_file)

  OPEN(unit  =  1,file  =  pos_file)
  print *, "# index_S2c, d_SH-d_HO (cv1), d_SH+d_HO (cv2), n_hydroxyl"
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS(step)
  USE parameters_pt, ONLY : nTi5c, index_S2c, Count_Hydroxyl, &
                            compute_cvs
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: step
  INTEGER                    :: iat
  INTEGER                    :: n_hydroxyl
  REAL*8                     :: cv1,cv2

  PRINT *, "####", step

  n_hydroxyl = Count_Hydroxyl(index_S2c)

  DO iat=1,nTi5c 

    CALL compute_cvs(index_S2c(iat),cv1,cv2)

    WRITE(*,fmt="(i2,2(3X,E17.5),3X, i3)"), & 
      iat, cv1, cv2, n_hydroxyl

  END DO

END SUBROUTINE ANALYSIS

