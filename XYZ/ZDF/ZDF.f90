!Computes the number density at different Z positions

MODULE histogram 
  IMPLICIT NONE 
  INTEGER                                  :: nhist, zdir
  INTEGER, ALLOCATABLE                     :: hist(:)
  CHARACTER(5)                             :: ztype
END MODULE histogram 

PROGRAM ZDF
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
      CALL READ_EXTXYZ
    IF (MOD(frame,nskip)==0) THEN
      !CALL SET_CENTER_OF_MASS_TO_ZERO
      CALL MAKE_HISTOGRAM
    END IF
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : hist, nhist, zdir, ztype
  USE parameters, ONLY : natoms, nframes, nequil, box,&
                         atype, pos, nskip, zoffset
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nskip, nhist, zdir
  READ(1,*) zoffset
  READ(*,*) ztype
  CLOSE(1)

  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(hist(nhist)); hist=0

  OPEN(unit  =  1,file  =  pos_file)
  hist = 0

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE histogram, ONLY  : hist, nhist, zdir, ztype
  USE parameters, ONLY : pos, box, natoms, atype, zoffset
  IMPLICIT NONE
  INTEGER                    :: i, ind
  REAL*8           :: z

  DO i = 1,natoms

    IF (atype(i) == ztype) THEN
      !Applying PBC
      z = pos(zdir,i) + zoffset
      z = z  - nint(z/box(zdir,zdir))*box(zdir,zdir) 
      z = z + box(zdir,zdir)/2. ! for data analysis only

      !Assigning the index to the histogram
      ind = int( z * float(nhist) / box(zdir,zdir) ) + 1
      hist(ind) = hist(ind) + 1
    END IF

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram, ONLY  : nhist, hist, zdir
  USE parameters, ONLY : box, nframes, nskip
  IMPLICIT NONE
  INTEGER                    :: i
  REAL*8                     :: bin,dV, M33DET

  OPEN(unit = 2,file = "ZDF.dat")

  bin = box(zdir,zdir)/float(nhist)
  dV = M33DET(box)
  dV = dV/box(zdir,zdir)*bin

    DO i = 1,nhist
      WRITE(2,fmt = "(F10.5, 3X, F12.7)") bin/2.d0 + float(i-1)*bin - box(zdir,zdir)/2.d0, &
                           & float(hist(i))/(0.03316d0*dV*float(nframes)/float(nskip))
    END DO

  CLOSE(2)

END SUBROUTINE
