!Computes the number density at different Z positions
!Note that the xyz trajectory should be sorted by atoms
!Please check the order of atom index in the xyz file 
!in order to know the index of each element

MODULE histogram 
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbins = 200
  DOUBLE PRECISION, PARAMETER              :: D_conv=0.03316
  INTEGER, ALLOCATABLE                     :: hist(:,:)
  INTEGER                                  :: ntype
END MODULE histogram 

PROGRAM ZDF
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_POS
    !CALL SET_CENTER_OF_MASS_TO_ZERO
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : hist, ntype, nbins
  USE parameters, ONLY : natoms, nframes, nequil, box,&
                         ind_atom, pos
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)
  ntype = 4

  CALL READ_INDEX ( index_file )

  ALLOCATE(pos(natoms,3))
  ALLOCATE(hist(ntype,nbins)); hist=0

  OPEN(unit  =  1,file  =  pos_file)
 
  hist = 0

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE histogram, ONLY  : hist, nbins
  USE parameters, ONLY : pos, box, natoms, ind_atom
  IMPLICIT NONE
  INTEGER                    :: i, ind
  DOUBLE PRECISION                     :: z

  DO i = 1,natoms

      !Applying PBC
      z = pos(i,3) !+ box(3,3)/2. ! for data analysis only
      !z = pos(i,3)  - nint(pos(i,3)/box(3,3))*box(3,3) 
      z = z  - nint(z/box(3,3))*box(3,3) 
      z = z + box(3,3)/2. ! for data analysis only

      !print *, pos(i,1), pos(i,2), pos(i,3)
      !Assigning the index to the histogram
      ind = int( z * DBLE(nbins) / box(3,3) ) + 1
      if ((ind.le.0).or.(ind.gt.nbins)) THEN
        PRINT *, "There is something wrong with your input file"
        STOP
      END IF
      hist(ind_atom(i),ind) = hist(ind_atom(i),ind) + 1

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram, ONLY  : nbins, hist
  USE parameters, ONLY : box, nframes, inversion_symmetry
  IMPLICIT NONE
  INTEGER                    :: i
  DOUBLE PRECISION                     :: bin,dV

  OPEN(unit = 2,file = "ZDF.dat")

  bin = box(3,3)/DBLE(nbins)
  dV = box(1,1)*box(2,2)*bin

  IF ( inversion_symmetry ) THEN
    DO i = 1,nbins/2
      WRITE(2,fmt = "(F10.5, 5(3X, F12.7))"), bin/2.d0 + DBLE(i-1)*bin, &
                           & DBLE( hist(:,i)+hist(:,nbins-i+1) ) / 2.d0 / &
                           & ( 0.03316d0*dV*DBLE(nframes) )
    END DO
  ELSE
    DO i = 1,nbins
!      WRITE(2,fmt = "(F10.5, 5(3X, F12.7))"), bin/2.d0 + DBLE(i-1)*bin - box(3,3)/2.d0, &
!                           & DBLE(hist(:,i))/(0.03316d0*dV*DBLE(nframes))
      WRITE(2,fmt = "(F10.5, 5(3X, F12.7))"), bin/2.d0 + DBLE(i-1)*bin, &
                           & DBLE(hist(:,i))/(0.03316d0*dV*DBLE(nframes))
    END DO
  END IF

  CLOSE(2)

END SUBROUTINE



