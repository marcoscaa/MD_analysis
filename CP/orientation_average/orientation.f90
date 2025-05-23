
MODULE histogram
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbin_h = 200
  DOUBLE PRECISION                         :: hist(nbin_h)
  INTEGER*8                                :: neq(nbin_h)
END MODULE histogram

PROGRAM RotCorr
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_POS
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RotCorr 

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : hist, neq
  USE parameters, ONLY : pos, natoms
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  CALL READ_INDEX( index_file ) 

  OPEN(unit  =  1,file  =  pos_file)
 
  ALLOCATE(pos(natoms,3))
  hist = 0.d0
  neq = 0

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE histogram,  ONLY : hist, nbin_h, neq
  USE parameters, ONLY : ind_atom, natoms, pos, box
  IMPLICIT NONE
  INTEGER                    :: i, j, ind
  INTEGER, DIMENSION(2)      :: ind_H
  DOUBLE PRECISION           :: theta, dipole, z
  LOGICAL                    :: is_water_oxygen

  DO i = 1,natoms

    !Selecting only Ow 
    IF ( is_water_oxygen(i) ) THEN
      
      !Applying PBC
      z = pos(i,3)  !- nint(pos(i,3)/box(3,3))*box(3,3) 
      !z = z + box(3,3)/2. ! for data analysis only
      z = z  - nint(z/box(3,3))*box(3,3) 
      z = z + box(3,3)/2. ! for data analysis only

      !Index for z layers histogram
      ind = int( z * DBLE(nbin_h) / box(3,3) ) + 1

      !Get the indexes of the H's bound to Ow
      ind_H=0
      CALL get_H(i, ind_H)

      IF ( (ind_H(1) .ne. 0) .and. (ind_H(2) .ne. 0) ) THEN 
  
        theta = dipole(i, ind_H(1),ind_H(2), 3)
        hist(ind) = hist(ind) + theta 
        neq(ind) = neq(ind) + 1 ! number of equivalent molecules to be averaged 
 
      ELSE
  
        PRINT *, "Some Oxygens have no H"
      END IF

    END IF

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram,  ONLY : hist, neq, nbin_h
  USE parameters, ONLY : box, nframes 
  IMPLICIT NONE
  INTEGER                    :: i, j
  DOUBLE PRECISION           :: aver, bin

  OPEN(unit  =  2,file  =  "orientationXz.dat")
  OPEN(unit  =  3,file  =  "debug_n.dat")

  bin = box(3,3) / dble(nbin_h)

  DO j = 1,nbin_h
    IF ( neq(j) > 0 ) THEN
      aver = hist(j) / dble( neq(j) ) 
    ELSE
      aver = 0
    END IF

    WRITE(2,fmt = "(F10.5, 3X, E12.5)"), bin/2.d0 + bin*dble(j-1), aver
    WRITE(3,fmt = "(F10.5, 3X, I7)"), bin/2.d0 + bin*dble(j-1), neq(j)

  END DO

  CLOSE(2)
  CLOSE(3)

END SUBROUTINE

