!Angle and bond length distribution. Length unis in nm

MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER             :: nbins = 300
  REAL*8, ALLOCATABLE            :: pos(:,:)
  REAL*8                         :: box(3)
  REAL*8, ALLOCATABLE            :: hist_angle(:), hist_bond(:)
  REAL*8, ALLOCATABLE            :: hist_angle2(:), hist_bond2(:)
  INTEGER, ALLOCATABLE           :: ind_atom(:), ind_mol(:), n_water(:) 
  INTEGER                        :: natoms, nframes, nequil
  LOGICAL, PARAMETER             :: inversion_symmetry=.false. 

END MODULE parameters

PROGRAM angbond
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM angbond

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(ind_atom(natoms))
  ALLOCATE(ind_mol(natoms))
  ALLOCATE(hist_angle(nbins)); hist_angle = 0.d0
  ALLOCATE(hist_bond(nbins)); hist_bond = 0.d0
  ALLOCATE(hist_angle2(nbins)); hist_angle2 = 0.d0
  ALLOCATE(hist_bond2(nbins)); hist_bond2 = 0.d0
  ALLOCATE(n_water(nbins)); n_water=0

  CLOSE(1)

  OPEN(unit = 1,file = file_name)
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters, ONLY : pos, ind_mol, ind_atom, natoms, box
  IMPLICIT NONE
  INTEGER                    :: i, junk
  CHARACTER(100)             :: atom_type, junk2

  READ(1,*)
  READ(1,*) natoms

  DO i = 1,natoms
    READ(1, fmt='(i5,2a5,i5,3f8.3)'), ind_mol(i), junk2, atom_type, junk, &
                                      pos(i,1), pos(i,2), pos(i,3)
    CALL convert_atom_type(atom_type,ind_atom(i))
  END DO

  READ(1,*) box(1), box(2), box(3)

END SUBROUTINE

SUBROUTINE MAKE_HISTOGRAM
  USE parameters, ONLY : pos, box, natoms, ind_atom, hist_bond, &
                         hist_angle, n_water, nbins, &
                         hist_angle2, hist_bond2
  IMPLICIT NONE
  INTEGER          :: iat, ind, i_bond, ind_H(2)
  REAL*8           :: z, angle_h2o, OH_len
  REAL*8           :: Angle, Dist

  DO iat = 1,natoms

    !Selecting only oxygen atoms from water
    IF (ind_atom(iat) == 1) THEN

      !Z position of water oxygen
      z = pos(iat,3)  - nint(pos(iat,3)/box(3))*box(3) 
      z = z + box(3)/2. ! for data analysis only

      !Assigning the water index to the histogram
      ind = int( z * float(nbins) / box(3) ) + 1
      n_water(ind) = n_water(ind) + 1

      !For this particular case, H atoms of a water molecule appear just after
      !the water oxygen atom
      ind_H(1) = iat+1
      ind_H(2) = iat+2

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

REAL*8 FUNCTION Angle(ind1, ind2, ind3)
  !Angle between particles with index ind1-ind2 and ind1-ind3
  !Return value in cos(Angle)
  USE parameters, ONLY : natoms, pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind1, ind2, ind3
  INTEGER                    :: ipol
  REAL*8                     :: v1(3), v2(3)

  DO ipol = 1,3

    v1(ipol) = pos(ind2,ipol) - pos(ind1,ipol)
    v1(ipol) = v1(ipol) - nint( v1(ipol) / box(ipol) ) * box(ipol)
    v2(ipol) = pos(ind3,ipol) - pos(ind1,ipol)
    v2(ipol) = v2(ipol) - nint( v2(ipol) / box(ipol) ) * box(ipol)

  END DO

  Angle = dot_product(v1,v2) / ( norm2(v1) * norm2(v2) )

END FUNCTION Angle

REAL*8 FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  REAL*8                     :: xyz(3)
  INTEGER                    :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos(ind1,i) - pos(ind2,i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

SUBROUTINE REMOVE_EQUIL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+3)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE PRINT_RESULTS
  USE parameters
  IMPLICIT NONE
  INTEGER          :: i
  REAL*8           :: aver, std, bin 

  OPEN(unit = 2,file = "angle.dat")
  OPEN(unit = 3,file = "bond.dat")

  bin = box(3)/float(nbins)

  IF ( inversion_symmetry ) THEN
    DO i = 1,nbins/2

      !Only selecting bins with enough sample
      IF ( ( n_water(i) == 0  ) .or. ( n_water(nbins-i+1) == 0) ) CYCLE

      aver = ( hist_angle(i) + hist_angle(nbins-i+1) ) / &
              float( n_water(i) + n_water(nbins-i+1) ) 
      std = ( hist_angle2(i) + hist_angle2(nbins-i+1) ) / &
             float( n_water(i) + n_water(nbins-i+1) ) 
      std = SQRT( std - aver*aver )

      WRITE(2,fmt = "(E12.5, 2(3X, E20.12))"), bin/2.d0 + float(i-1)*bin, &
                                               aver, std

      aver = ( hist_bond(i) + hist_bond(nbins-i+1) ) / &
              float( n_water(i) + n_water(nbins-i+1) ) / 2.d0
      std = ( hist_bond2(i) + hist_bond2(nbins-i+1) ) / &
             float( n_water(i) + n_water(nbins-i+1) ) / 2.d0
      std = SQRT( std - aver*aver )

      WRITE(3,fmt = "(E12.5, 2(3X, E20.12))"), bin/2.d0 + float(i-1)*bin, &
                                            aver, std 
    END DO
  ELSE
    DO i = 1,nbins

      !Only selecting bins with proper statistics
      IF ( n_water(i) == 0 ) CYCLE

      aver = hist_angle(i) / float( n_water(i) ) 
      std =  hist_angle2(i) / float( n_water(i) ) 
      std = SQRT( std - aver*aver )

      WRITE(2,fmt = "(f9.5, 2(3X, f12.8),2X,f12.2)"), bin/2.d0 + float(i-1)*bin,& !- box(3)/2.d0, &
                                            aver, std, float(n_water(i))

      aver = hist_bond(i) / float( n_water(i) ) / 2.D0 
      std =  hist_bond2(i) / float( n_water(i) ) / 2.d0 
      std = SQRT( std - aver*aver )

      WRITE(3,fmt = "(f9.5, 2(3X, f12.8))"), bin/2.d0 + float(i-1)*bin,& !- box(3)/2.d0, &
                                            aver, std
    END DO
  END IF

  CLOSE(2);CLOSE(3)

END SUBROUTINE

SUBROUTINE convert_atom_type(atom_type, atom_ind)
  !Convert char type to int type
  IMPLICIT NONE
  CHARACTER(100), INTENT(in)      :: atom_type
  INTEGER, INTENT(INOUT)          :: atom_ind

  SELECT CASE ( trim(atom_type) )
    CASE( "   OW" )
      atom_ind = 1
    CASE( "  HW1" )
      atom_ind = 2
    CASE( "  HW2" )
      atom_ind = 2
  END SELECT

END SUBROUTINE convert_atom_type
