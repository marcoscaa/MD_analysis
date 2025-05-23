
MODULE histogram
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: oh_dist(:,:)
  INTEGER, PARAMETER                       :: nbins=30
END MODULE histogram

PROGRAM OHcorr
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_POS
    CALL IDENTIFY_OH_GROUPS
    CALL OH_DISTRIBUTION 
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM OHcorr

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : oh_dist, nbins
  USE parameters, ONLY : natoms, nframes, &
                         nOxygen, pos, coarse, &
                         nlayers, OH_index
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  CALL READ_INDEX ( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(oh_dist(nlayers,nbins)); oh_dist = 0.d0
  ALLOCATE(OH_index(natoms)); OH_index =0

  OPEN(unit = 1,file = pos_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE OH_DISTRIBUTION 
  !Build oh_dist matrix. oh_dist contains the normalized OH vectors 
  !for all OH groups in the system
  USE histogram,  ONLY : oh_dist, nbins 
  USE parameters, ONLY : natoms, nwater
  IMPLICIT NONE
  INTEGER                    :: i, ih, layer
  INTEGER                    :: ind_H(2)
  INTEGER                    :: ind, ind_layer
  INTEGER                    :: n_hydrogens, i_Ox, nOH
  INTEGER                    :: ihp 
  LOGICAL                    :: is_oxygen, switched_oh_index
  DOUBLE PRECISION           :: oh(2,3)

  !OH atom index
  nOH = 0
  i_Ox = 0 

  DO i = 1,natoms

    !Calculate the OH vector:
    IF ( is_oxygen(i) ) THEN

      i_Ox = i_Ox + 1
      CALL GET_IND_H(i, ind_H, n_hydrogens)

      IF ( n_hydrogens/=0 ) THEN 

        layer = ind_layer( i )

        CALL OH_VECT( i, ind_H, oh)

        DO ih = 1,n_hydrogens

          nOH = nOH + 1
          ind = int( nbins*( oh(ih,3) + 1. ) / 2. ) + 1 
          oh_dist(layer,ind) = oh_dist(layer,ind) + 1 
            
        END DO

      END IF !n_hydrogens

    END IF !is_oxygen

  ENDDO

  IF ( nOH /= 2*nwater ) THEN
    print *, 'WARNING: Wrong number of OH groups ', nOH
    !STOP
  END IF

END SUBROUTINE OH_DISTRIBUTION 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : oh_dist, nbins
  USE parameters, ONLY : nlayers 
  IMPLICIT NONE
  INTEGER                    :: i, j

  OPEN( unit = 2, file = "OH_dist.dat" )

  DO i = 1,nlayers

    oh_dist(i,:) = oh_dist(i,:) / sum(oh_dist(i,:))

  END DO

  DO i = 1,nbins

    write(2,fmt = '(E11.4,3X,*(E14.7,3X))'), -1.d0 + (dble(i-1)+0.5d0)*2.0/dble(nbins), &
         ( oh_dist(j,i) , j=1,nlayers) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
