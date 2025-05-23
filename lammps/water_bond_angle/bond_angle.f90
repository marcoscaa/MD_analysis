MODULE histogram 
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE         :: ntot
  REAL*8, ALLOCATABLE          :: OH(:), HOH(:) 
  REAL*8, PARAMETER            :: rmax=1.3d0
  INTEGER                      :: typeO, typeH
END MODULE histogram

PROGRAM Orientation
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    if (MOD(frame,stride)==0) CALL ANALYSIS
  END DO

  CALL PRINT_RESULTS

END PROGRAM Orientation

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist
  USE histogram, ONLY : OH,HOH, ntot, rmax, typeO, typeH
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, stride, nhist
  READ(1, *) typeO, typeH 
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(OH(nhist)); OH=0.0
  ALLOCATE(HOH(nhist)); HOH=0.0

  ntot=0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS
  USE parameters, ONLY : pos, atype, natoms, nhist 
  USE histogram, ONLY : OH,HOH, ntot,typeO,rmax
  IMPLICIT NONE
  INTEGER                    :: i, ind, iH
  INTEGER, DIMENSION(2)      :: ind_H 
  REAL*8                     :: Dist, Angle
  REAL*8                     :: d_, a_

  DO i = 1,natoms

    !Selecting only Ow
    IF ( atype(i) == typeO )  THEN

      ntot = ntot+1
      CALL get_H(i, ind_H)
      a_ = Angle(pos(:,i),pos(:,ind_H(1)),pos(:,ind_H(2)))
      ind = (a_+1.0)/2.*float(nhist)
      HOH(ind)=HOH(ind)+1

      DO iH=1,2
        d_ = Dist(i,ind_H(iH))
        ind = d_/rmax * float(nhist)
        OH(ind) = OH(ind)+1
      END DO

    END IF !

  END DO !i 

END SUBROUTINE ANALYSIS

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : atype, pos, cutoff_OH, natoms
  USE histogram, ONLY: typeH
  IMPLICIT NONE
  INTEGER, DIMENSION(2)      :: ind_H
  INTEGER                    :: i, k, ind_O
  REAL*8                     :: Dist

  i=1;k=1

  ind_H = 0
  DO WHILE (( i <= natoms) .and. ( k <= 2))

    !Should be the index for H
    IF ( atype(i) == typeH ) THEN

      IF ( Dist(ind_O,i) < cutoff_OH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

  IF (k.ne.3) THEN
    print *, 'Found water ion with charge ', k-3
  END IF
END SUBROUTINE get_H

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist
  USE histogram, ONLY : rmax,OH,HOH,ntot,rmax
  IMPLICIT NONE
  INTEGER :: ihist

  OPEN(unit = 2, file = "OH_bond_distribution.dat")
  OPEN(unit = 3, file = "HOH_angle_distribution.dat")

  WRITE(2, *) "# OH_bond(Angstrom) Probability"
  WRITE(3, *) "# HOH_angle(cos) Probability"

  DO ihist=1,nhist
    WRITE(2,fmt='(F12.8,3X,F12.8)') (float(ihist))*rmax/float(nhist), &
        OH(ihist)/float(ntot)/2. !2 OH bonds per water
    WRITE(3,fmt='(F12.8,3X,F12.8)') 2*float(ihist)/float(nhist)-1, &
        HOH(ihist)/float(ntot)
  END DO

  CLOSE(2)
  CLOSE(3)

END SUBROUTINE PRINT_RESULTS
