
MODULE histogram
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbins = 250
  DOUBLE PRECISION, ALLOCATABLE            :: hist_angle(:), hist_bond(:)
  DOUBLE PRECISION, ALLOCATABLE            :: hist_angle2(:), hist_bond2(:)
  INTEGER, ALLOCATABLE                     :: n_water(:)

END MODULE histogram 

PROGRAM angbond
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

END PROGRAM angbond

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, pos
  USE histogram,  ONLY : hist_angle, hist_bond, hist_angle2, &
                         hist_bond2, n_water, nbins
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  CALL READ_INDEX( index_file )
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(hist_angle(nbins)); hist_angle = 0.d0
  ALLOCATE(hist_bond(nbins)); hist_bond = 0.d0
  ALLOCATE(hist_angle2(nbins)); hist_angle2 = 0.d0
  ALLOCATE(hist_bond2(nbins)); hist_bond2 = 0.d0
  ALLOCATE(n_water(nbins)); n_water=0

  OPEN(unit = 1,file = file_name)
 
END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE parameters, ONLY : pos, box, natoms, ind_atom
  USE histogram,  ONLY : hist_angle, hist_angle2, hist_bond, &
                         hist_bond2, n_water, nbins
  IMPLICIT NONE
  INTEGER                    :: iat, ind, i_bond, ind_H(2)
  DOUBLE PRECISION           :: z, angle_h2o, OH_len
  DOUBLE PRECISION           :: Angle, Dist

  DO iat = 1,natoms

    !Selecting only oxygen atoms from water
    IF (ind_atom(iat) == 4) THEN

      !Z position of water oxygen
      z = pos(iat,3)  - nint(pos(iat,3)/box(3,3))*box(3,3) 
      z = z + box(3,3)/2. ! for data analysis only

      !Assigning the water index to the histogram
      ind = int( z * DBLE(nbins) / box(3,3) ) + 1
      n_water(ind) = n_water(ind) + 1

      !Get the indexes of the two H bond to OW
      CALL get_H(iat, ind_H)

      !Bond angle of water molecule
      angle_h2o = Angle(iat, ind_H(1), ind_H(2))
      hist_angle(ind) = hist_angle(ind) + angle_h2o 
      hist_angle2(ind) = hist_angle2(ind) + angle_h2o*angle_h2o 

      !OH bond length in water
      DO i_bond = 1,2
        OH_len = Dist(iat,ind_H(i_bond))
        hist_bond(ind) = hist_bond(ind) + OH_len 
        hist_bond2(ind) = hist_bond2(ind) + OH_len*OH_len
      END DO

    END IF

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : box, inversion_symmetry
  USE histogram,  ONLY : n_water, nbins, hist_angle, hist_angle2, &
                         hist_bond, hist_bond2
  IMPLICIT NONE
  INTEGER                    :: i
  DOUBLE PRECISION           :: aver, std, bin 

  OPEN(unit = 2,file = "angle.dat")
  OPEN(unit = 3,file = "bond.dat")

  bin = box(3,3)/DBLE(nbins)

  IF ( inversion_symmetry ) THEN
    DO i = 1,nbins/2

      IF ( ( n_water(i) == 0 ) .or. ( n_water(nbins-i+1) == 0) ) CYCLE

      aver = ( hist_angle(i) + hist_angle(nbins-i+1) ) / &
              DBLE( n_water(i) + n_water(nbins-i+1) ) 
      std = ( hist_angle2(i) + hist_angle2(nbins-i+1) ) / &
             DBLE( n_water(i) + n_water(nbins-i+1) ) 
      std = SQRT( std - aver*aver )

      WRITE(2,fmt = "(E12.5, 2(3X, E20.12))"), bin/2.d0 + DBLE(i-1)*bin, &
                                               aver, std

      aver = ( hist_bond(i) + hist_bond(nbins-i+1) ) / &
              DBLE( n_water(i) + n_water(nbins-i+1) ) / 2.d0
      std = ( hist_bond2(i) + hist_bond2(nbins-i+1) ) / &
             DBLE( n_water(i) + n_water(nbins-i+1) ) / 2.d0
      std = SQRT( std - aver*aver )

      WRITE(3,fmt = "(E12.5, 2(3X, E20.12))"), bin/2.d0 + DBLE(i-1)*bin, &
                                            aver, std 
    END DO
  ELSE
    DO i = 1,nbins

      IF ( n_water(i) == 0 ) CYCLE

      aver = hist_angle(i) / DBLE( n_water(i) ) 
      std =  hist_angle2(i) / DBLE( n_water(i) ) 
      std = SQRT( std - aver*aver )

      WRITE(2,fmt = "(E12.5, 2(3X, E20.12))"), bin/2.d0 + DBLE(i-1)*bin - box(3,3)/2.d0, &
                                            aver, std 

      aver = hist_bond(i) / DBLE( n_water(i) ) / 2.D0 
      std =  hist_bond2(i) / DBLE( n_water(i) ) / 2.d0 
      std = SQRT( std - aver*aver )

      WRITE(3,fmt = "(E12.5, 2(3X, E20.12))"), bin/2.d0 + DBLE(i-1)*bin - box(3,3)/2.d0, &
                                            aver, std
    END DO
  END IF

  CLOSE(2);CLOSE(3)

END SUBROUTINE
