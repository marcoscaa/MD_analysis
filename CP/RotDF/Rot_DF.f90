
MODULE histogram
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbin_h = 30
  INTEGER*8, ALLOCATABLE                   :: hist(:,:)
  INTEGER*8, ALLOCATABLE                   :: neq(:)
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
  USE histogram,  ONLY : hist, neq, nbin_h
  USE parameters, ONLY : pos, natoms, nlayers
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  CALL READ_INDEX( index_file ) 

  OPEN(unit  =  1,file  =  pos_file)
 
  ALLOCATE(pos(natoms,3))
  ALLOCATE(hist(nlayers,nbin_h)); hist = 0
  ALLOCATE(neq(nlayers)); neq = 0

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE histogram,  ONLY : hist, neq, nbin_h
  USE parameters, ONLY : ind_atom, natoms, pos, box, nlayers
  IMPLICIT NONE
  INTEGER                    :: i, j, ind_hist, ind_layer
  INTEGER                    :: layer_index
  INTEGER, DIMENSION(2)      :: ind_H
  DOUBLE PRECISION           :: theta, dipole
  LOGICAL                    :: is_water_oxygen

  DO i = 1,natoms

    !Selecting only Ow 
    IF ( is_water_oxygen(i) ) THEN
      
      !Index for z layers histogram
      ind_layer = layer_index( pos(i,3) ) 

      !Not consider the particles outised the bins range
      IF ( ind_layer == 0 ) CYCLE

      !Get the indexes of the H's bound to Ow
      ind_H=0
      CALL get_H(i, ind_H)

      IF ( (ind_H(1) .ne. 0) .and. (ind_H(2) .ne. 0) ) THEN 
  
        theta = dipole(i, ind_H(1),ind_H(2), 3)
        ind_hist = int( (theta + 1.) * DBLE(nbin_h) / 2. ) + 1
        hist(ind_layer,ind_hist) = hist(ind_layer,ind_hist) + 1
        neq(ind_layer) = neq(ind_layer) + 1 ! number of equivalent molecules to be averaged 
 
      END IF

    END IF

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram,  ONLY : hist, neq, nbin_h
  USE parameters, ONLY : nlayers
  IMPLICIT NONE
  INTEGER                    :: i, j
  DOUBLE PRECISION           :: aver(nlayers)

  OPEN(unit  =  2,file  =  "orientation.dat")

  DO j = 1,nbin_h
    aver = float( hist(:,j) ) / float( neq(:) ) 

    WRITE(2,fmt = "(F10.5, *(3X, F12.7))"), ( float(j) - 0.5 ) *  & 
                                    2. / float(nbin_h) - 1., aver(:)

  END DO

  CLOSE(2)

END SUBROUTINE

