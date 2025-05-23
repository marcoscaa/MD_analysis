!Computes the average and instantaneous Hbond
!Returns these values with respect to layers in the z axis

MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbins = 25
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: pos
  REAL*8, DIMENSION(9)                     :: box
  REAL*8, PARAMETER                        :: cutoff_OH = 1.3D0 
  INTEGER*8, DIMENSION(:), ALLOCATABLE     :: ind_atom
  INTEGER*8, DIMENSION(nbins)              :: n_inlayer
  REAL*8, DIMENSION(nbins)                 :: n_H, n_H2
  INTEGER                                  :: natoms, nframes
  INTEGER                                  :: frame, nequil

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE

  CALL INITIALIZE

  !Remove equilibration part
  DO frame = 1,nequil*(natoms+2)
    READ(1, *)
  END DO
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  READ(1, *), box
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(ind_atom(natoms))

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
  END DO
 
  !Initialize the counters
  n_H=0
  n_H2=0
  n_inlayer=0

  CLOSE(1)

  OPEN(unit = 1,file = file_name)

END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(5)               :: junk

  READ(1,*)
  READ(1,*)

  DO i = 1,natoms
    READ(1, *), junk, pos(i,1), pos(i,2), pos(i,3)
  END DO

  !Shifting everything up in the z direction
  pos(:,3) = pos(:,3) + box(9)/2.

END SUBROUTINE

SUBROUTINE MAKE_HISTOGRAM
  USE parameters
  IMPLICIT NONE
  INTEGER                      :: i, j, ind_i, ind_j
  INTEGER, DIMENSION(2)        :: ind_H 
  INTEGER*8, DIMENSION(nbins)  :: hist, n_in
  REAL*8                       :: d_OOw, Dist, Angle
  REAL*8, DIMENSION(2)         :: ang_OHO

  hist = 0
  n_in = 0

  DO i = 1,natoms

    !Selecting only O from H2O 
    IF ( ind_atom(i) == 4 ) THEN

      !Counting total number of atoms inside a layer
      ind_i = int( pos(i,3) * DBLE(nbins) / box(9) ) + 1
      n_in(ind_i) = n_in(ind_i) + 1
      !Get the indexes of H bonded to Ow
      ind_H = 0
      CALL get_H(i, ind_H)

      DO j = 1,natoms

        !Selecting only O2c from TiO2 and OW (water)
        IF ( (ind_atom(j) == 3) .or. (ind_atom(j) == 4) ) THEN

          !Getting the layer index for atom j 
          ind_j = int( pos(j,3) * DBLE(nbins) / box(9) ) + 1

          !Calculate distace between Ow - Ow 
          d_OOw = Dist(pos,i,j)

          !Definition of H-bond
          IF ( (d_OOw < 3.5) .and.  (ind_H(1) /= 0) .and. (ind_H(2) /= 0) ) THEN

            !Calculates the H-bond angle
            ang_OHO(1) = Angle(pos,i,j,ind_H(1))
            ang_OHO(2) = Angle(pos,i,j,ind_H(2))

            !Definition of H-bond (includes also d < 3.5) - Luzar and Chandler
            IF (maxval(ang_OHO) > 0.8660) THEN

              hist(ind_i) = hist(ind_i) + 1
              IF (ind_atom(j) == 3) hist(ind_j) = hist(ind_j) + 1

            END IF !d_OOw
 
          END IF !ind_H

        END IF !ind_atom

      END DO !j

    END IF !ind_atom

  END DO !i 
  do i=1,nbins
    IF (n_in(i) .ne. 0) THEN
      n_H(i) = n_H(i) + DBLE(hist(i))/DBLE(n_in(i))
      n_H2(i) = n_H2(i) + DBLE(hist(i)*hist(i))/DBLE(n_in(i)**2)
      n_inlayer(i) = n_inlayer(i) + 1
    ENDIF
  enddo

END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE parameters
  IMPLICIT NONE
  REAL*8                     :: mean, stdev, bin, dV
  INTEGER                    :: i

  OPEN(unit = 2,file = "Hbond_layer.dat")

  bin = box(9)/DBLE(nbins)
  dV = box(1)*box(5)*bin

  DO i = 1,nbins
    IF ( n_inlayer(i) .ne. 0 ) THEN
      mean = n_H(i) / DBLE(n_inlayer(i))
      stdev = SQRT( n_H2(i) / DBLE(n_inlayer(i)) - mean**2 )
    ELSE
      mean = 0.
      stdev = 0.
    ENDIF

    WRITE(2,fmt = "(F10.5, 5(3X, F12.7))"), bin/2 + DBLE(i-1)*bin, mean, stdev
  END DO

  CLOSE(2)

END SUBROUTINE
 
SUBROUTINE get_H(ind_O, ind_H)
  USE parameters
  IMPLICIT NONE
  INTEGER, DIMENSION(2)      :: ind_H
  INTEGER                    :: i, k, ind_O
  REAL*8                     :: Dist

  i=1;k=1

  DO WHILE (( i <= natoms) .and. ( k <= 2))

    !Should be the index for H
    IF ( ind_atom(i) == 2 ) THEN

      IF ( Dist(pos,ind_O,i) < cutoff_OH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

END SUBROUTINE get_H

REAL*8 FUNCTION Dist(pos, ind1, ind2)
  ! Distance between particles with indexes 1 and 2
  USE parameters, ONLY : box, natoms
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: pos(natoms,3)
  INTEGER, INTENT(IN)        :: ind1, ind2
  REAL*8                     :: xyz(3)
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = pos(ind1,i) - pos(ind2,i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(4*i-3) ) * box(4*i-3)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

REAL*8 FUNCTION Dist2(pos, ind1, ind2)
  ! Squared distance between two particles including pbc
  USE parameters, ONLY : box, natoms
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: pos(natoms,3)
  INTEGER, INTENT(IN)        :: ind1, ind2
  REAL*8                     :: xyz(3)
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = pos(ind1,i) - pos(ind2,i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(4*i-3) ) * box(4*i-3)

  END DO

  Dist2 = SUM(xyz*xyz)

END FUNCTION Dist2

REAL*8 FUNCTION Angle(pos,ind1, ind2, ind3)
  !Angle between particles with index ind1-ind2 and ind1-ind3
  !Return value in cos(Angle)
  USE parameters, ONLY : natoms
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: pos(natoms,3)
  INTEGER, INTENT(IN)        :: ind1, ind2, ind3
  REAL*8                     :: v1(3), v2(3)

  v1(:) = pos(ind2,:) - pos(ind1,:)
  v2(:) = pos(ind3,:) - pos(ind1,:)

  Angle = dot_product(v1,v2) / ( norm2(v1) * norm2(v2) )
  !Angle = Dist2(C,a1) + Dist2(C,a2) - Dist2(a1,a2)
  !Angle = Angle / ( 2 * Dist(C,a1) * Dist(C,a2) )

END FUNCTION Angle
