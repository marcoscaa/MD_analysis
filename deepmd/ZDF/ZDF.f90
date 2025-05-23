!Computes the number density at different Z positions

MODULE histogram 
  IMPLICIT NONE 
  INTEGER                                  :: nhist
  INTEGER, ALLOCATABLE                     :: hist(:,:)
  INTEGER                                  :: ntype
  DOUBLE PRECISION                         :: zcenter
END MODULE histogram 

PROGRAM ZDF
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_RAW_POS_BOX
    IF (MOD(frame,nskip)==0) THEN
      !CALL SET_CENTER_OF_MASS_TO_ZERO
      CALL MAKE_HISTOGRAM
    END IF
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : hist, ntype, nhist, zcenter
  USE parameters, ONLY : natoms, nframes, nequil, box,&
                         atype, pos, nskip
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, ntype, nframes, nequil, nskip, nhist
  READ(1,*) zcenter
  CLOSE(1)

  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(hist(ntype,nhist)); hist=0

  !Read the index file
  OPEN(unit = 1,file = 'type.raw')
  READ(1, *) (atype(i), i=1,natoms)
  CLOSE(1)

  OPEN(unit  =  1,file  =  'coord.raw')
  OPEN(unit  =  2,file  =  'box.raw')
 
  hist = 0

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE histogram, ONLY  : hist, nhist, zcenter
  USE parameters, ONLY : pos, box, natoms, atype
  IMPLICIT NONE
  INTEGER                    :: i, ind
  DOUBLE PRECISION           :: z

  DO i = 1,natoms

      !Applying PBC
      z= pos(3,i) + zcenter
      z = z  - nint(z/box(3,3))*box(3,3) 
      z = z + box(3,3)/2. ! for data analysis only

      !Assigning the index to the histogram
      ind = int( z * DBLE(nhist) / box(3,3) ) + 1
      hist(atype(i)+1,ind) = hist(atype(i)+1,ind) + 1

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram, ONLY  : nhist, hist
  USE parameters, ONLY : box, nframes, nskip
  IMPLICIT NONE
  INTEGER                    :: i
  DOUBLE PRECISION                     :: bin,dV

  OPEN(unit = 10,file = "ZDF.dat")

  bin = box(3,3)/DBLE(nhist)
  dV = box(1,1)*box(2,2)*bin

    DO i = 1,nhist
      WRITE(10,fmt = "(F10.5, 5(3X, F12.7))"), bin/2.d0 + DBLE(i-1)*bin - box(3,3)/2.d0, &
                           & DBLE(hist(:,i))/(0.03316d0*dV*DBLE(nframes)/DBLE(nskip))
    END DO

  CLOSE(10)

END SUBROUTINE
