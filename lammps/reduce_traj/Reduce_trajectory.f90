
MODULE parameters
  IMPLICIT NONE
  INTEGER :: natoms, nframes, nequil, nskip
END MODULE parameters

PROGRAM REDUCE_TRAJ
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  
  DO frame = 1,nframes
    IF ((frame.gt.nequil).and.(MOD(frame-nequil,nskip)==0)) THEN
      CALL READ_AND_PRINT_TRAJ
    ELSE
      CALL READ_ONLY      
    END IF
  END DO

  CLOSE(1);CLOSE(2)

END PROGRAM REDUCE_TRAJ

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, nskip
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nskip
  CLOSE(1)

  OPEN(unit  =  1,file  =  pos_file)
  OPEN(unit  =  2,file  =  "out.reduced")
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_AND_PRINT_TRAJ
  USE parameters, ONLY : natoms
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: line
   
  DO i =1,natoms+9

    READ(1,'(A)') line
    WRITE(2,'(100A)') line
    
  END DO  

END SUBROUTINE READ_AND_PRINT_TRAJ

SUBROUTINE READ_ONLY
  USE parameters, ONLY : natoms
  IMPLICIT NONE
  INTEGER                    :: i
   
  DO i =1,natoms+9
    READ(1,*) 
  END DO  

END SUBROUTINE READ_ONLY
