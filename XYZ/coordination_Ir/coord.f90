!Reads XYZ files and computes the number of O, OH and H2O 
!groups bound to a selected Ir atom.
!Run the code as:
!echo "Ir_index N_index" | ./coord.x traj.extxyz
!Ir_index is the index of the Ir atom to compute coordination number
!N_index is the index of the other atom used to compute the CV

PROGRAM CONVERT
 USE, intrinsic :: iso_fortran_env, Only : iostat_end
 INTEGER :: iostat, frame

 CALL INITIALIZE(frame)

 DO 
   CALL READ_EXTXYZ_IO (iostat)
   IF (iostat == iostat_end ) THEN
     EXIT
   END IF
   CALL COMPUTE_IR_COORDINATION
   CALL WRITE_OUTPUT(frame)
 END DO

 CALL FINALIZE

END PROGRAM

MODULE coordination
  IMPLICIT NONE 
  INTEGER           :: Ir_index, nOH(3), nOO, N_index
  REAL*8, PARAMETER :: cutoff_IrO=2.6d0, cutoff_OH=1.1d0
  REAL*8, PARAMETER :: cutoff_OO=1.6d0
END MODULE coordination

SUBROUTINE INITIALIZE(frame)
  USE coordination, only : Ir_index, N_index
  IMPLICIT NONE
  CHARACTER(100)             :: pos_file
  INTEGER :: i, frame

  CALL getarg(1, pos_file)
  READ *, Ir_index, N_index

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = "Ir_speciation.dat")
  WRITE(2, fmt="(A44)") "# Frame Ir-X_distance(A) n_O n_OH n_H2O n_OO"
  frame=0

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_IR_COORDINATION
  USE parameters, ONLY : natoms, atype
  USE coordination
  IMPLICIT NONE
  REAL*8                     :: d, Dist 
  INTEGER                    :: iOx, nH, CoordNumb

  nOH=0
  nOO=0

  DO iOx = 1,natoms

    IF ( trim(atype(iOx))=="O" ) THEN 

        d=Dist(iOx,Ir_index)

        IF (d < cutoff_IrO) THEN

          nH = CoordNumb(iOx,"H    ",cutoff_OH)
          nOO = nOO + CoordNumb(iOx,"O    ",cutoff_OO)
          nOH(nH+1) = nOH(nH+1) +1

        END IF

    END IF

  END DO

  !Avoid double counting O when OO is bound to Ir 
  nOH(1) = nOH(1) - nOO

END SUBROUTINE COMPUTE_IR_COORDINATION

SUBROUTINE WRITE_OUTPUT(frame)
  USE coordination, only: Ir_index, N_index, nOH, nOO
  IMPLICIT NONE
  REAL*8                    :: d, Dist
  INTEGER                   :: frame

  frame = frame+1
  d = Dist(Ir_index, N_index)

  WRITE(2, fmt="(I6,F14.10,4(1X,I4))") frame, d, nOH(1), nOH(2), nOH(3), nOO

END SUBROUTINE WRITE_OUTPUT

SUBROUTINE FINALIZE
  IMPLICIT NONE
  CLOSE(1);CLOSE(2)
END SUBROUTINE FINALIZE
