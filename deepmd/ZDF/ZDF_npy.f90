!Computes the number density at different Z positions

MODULE histogram 
  IMPLICIT NONE 
  INTEGER                                  :: nhist, counter
  INTEGER                                  :: ntype
  REAL, ALLOCATABLE            :: hist(:,:)
  REAL                         :: zcenter, zmax
END MODULE histogram 

PROGRAM ZDF
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: iset

  CALL INITIALIZE
  
  DO iset = 1,nsets
    CALL OPEN_SET(iset-1)
    CALL MAKE_HISTOGRAM
    CALL CLOSE_SET
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : hist, ntype, nhist, zcenter, counter, zmax
  USE parameters, ONLY : natoms, nframes, nequil, boxt,&
                         atype, nskip, nsets
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)
  ntype = 4

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nsets, nequil, nskip, nhist
  READ(1,*) zcenter, zmax
  CLOSE(1)

  ALLOCATE(atype(natoms))
  ALLOCATE(hist(ntype,nhist)); hist=0

  !Read the index file
  OPEN(unit = 1,file = 'type.raw')
  READ(1, *) (atype(i), i=1,natoms)
  CLOSE(1)

  hist = 0.d0
  counter = 0

END SUBROUTINE INITIALIZE

SUBROUTINE OPEN_SET(myset) 
  USE parameters, ONLY : natoms, boxt, post, nframes
  use iso_fortran_env
  IMPLICIT NONE
  INTEGER                    :: i,myset
  CHARACTER(100)              :: filename
  !REAL*4, ALLOCATABLE   :: p(:,:)

  !Read the index file
  write(filename, "(A4,I0.3,A8)") "set.", myset, "/nframes"
  print *, filename
  OPEN(unit=1, file = filename)
  READ(1,*) nframes
  CLOSE(1)

  ALLOCATE(post(3,natoms,nframes))
  ALLOCATE(boxt(9,nframes))

  write(filename, "(A4,I0.3,A10)") "set.", myset, "/coord.bin"
  OPEN(1, file=filename, action="read", form='unformatted', access='stream')
  write(filename, "(A4,I0.3,A8)") "set.", myset, "/box.bin"
  OPEN(2, file=filename, action="read", form='unformatted', access='stream')
  READ(1) post
  READ(2) boxt

  CLOSE(1);CLOSE(2)

END SUBROUTINE OPEN_SET

SUBROUTINE CLOSE_SET
  USE parameters, only : post,boxt
  IMPLICIT NONE

  DEALLOCATE(post,boxt)

END SUBROUTINE CLOSE_SET

SUBROUTINE MAKE_HISTOGRAM
  USE histogram, ONLY  : hist, nhist, zcenter, counter, zmax, ntype
  USE parameters, ONLY : post, boxt, natoms, atype, nskip, nframes
  IMPLICIT NONE
  INTEGER                    :: i, ind, n
  INTEGER                    :: hist_tmp(ntype,nhist)
  REAL           :: z, dV

  DO n = 1,nframes,nskip

    hist_tmp = 0

    DO i = 1,natoms

        !Applying PBC
        z= post(3,i,n) - zcenter
        z = z  - nint(z/boxt(9,n))*boxt(9,n) 
        z = z + boxt(9,n)/2. ! for data analysis only

        !Assigning the index to the histogram
        ind = int( z * float(nhist) / zmax ) + 1
        hist_tmp(atype(i)+1,ind) = hist_tmp(atype(i)+1,ind) + 1

    END DO

    dV = boxt(1,n)*boxt(5,n)*zmax/float(nhist)
    hist = hist + float(hist_tmp)/dV
    counter = counter + 1

  END DO

  print *, dV, counter, nframes, nskip
  print *, boxt(1,1), boxt(5,1), boxt(9,1)
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram, ONLY  : nhist, hist, counter, zmax
  USE parameters, ONLY : nframes, nskip
  IMPLICIT NONE
  INTEGER                    :: i
  REAL                     :: bin,dV

  OPEN(unit = 10,file = "ZDF.dat")

  bin = zmax/DBLE(nhist)

    DO i = 1,nhist
      WRITE(10,fmt = "(F10.5, 5(3X, F12.7))"), bin/2.d0 + DBLE(i-1)*bin - zmax/2.d0, &
                           & DBLE(hist(:,i))/(0.03316d0*DBLE(counter))
    END DO

  CLOSE(10)

END SUBROUTINE
