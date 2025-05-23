!This code assumes that the radial direction is along XY

MODULE wannier_mod 
  IMPLICIT NONE 
  INTEGER, PARAMETER           :: nhist_dipole=2000, nlayers=3
  INTEGER                      :: ntot(nlayers)
  INTEGER                      :: atypeO, atypeH, atypeC
  REAL*8, ALLOCATABLE          :: dipole(:,:) 
  REAL*8                       :: T, center(3), hist_dipole(7,nhist_dipole,nlayers)
  REAL*8                       :: layers(nlayers+1)
  REAL*8, PARAMETER            :: maxR=50.0d0
  REAL*8, PARAMETER            :: min_dipole=-20.d0, max_dipole=20.d0
  REAL*8, PARAMETER            :: convert_to_debye=1.d0/0.2081943d0
END MODULE wannier_mod

PROGRAM Dipole_Distribution
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  !CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL READ_WANNIER
    if (frame>100) then
      CALL FIND_CENTER_CNT
      CALL COMPUTE_WATER_DIPOLE
    endif
    CALL DISTRIBUTION_MOLECULAR_DIPOLE
  END DO

  CALL NORMALIZE_HISTOGRAM
  CALL PRINT_RESULTS

END PROGRAM Dipole_Distribution

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist, nwann, wannier
  USE wannier_mod, ONLY : dipole,atypeO, atypeH, atypeC, T, &
                          hist_dipole, maxr, ntot, layers
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, wfc_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, wfc_file)
  CALL getarg(3, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms,nwann, nframes, T
  READ(1,*) atypeO, atypeH, atypeC
  READ(1,*) layers
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(dipole(3,natoms)); dipole=0.d0
  ALLOCATE(wannier(3,natoms)); wannier=0.d0

  hist_dipole=0
  ntot=0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
  OPEN(unit = 3, file =  wfc_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_WATER_DIPOLE
  USE parameters, ONLY : pos, atype, natoms, wannier
  USE wannier_mod, ONLY : dipole, convert_to_debye, atypeO, &
                          atypeH
  IMPLICIT NONE
  INTEGER                    :: i, iH, ind_H(2)
  REAL*8                     :: d(4)

  dipole=0.d0

  DO i = 1,natoms

    !Selecting only Ow
    IF ( atype(i) == atypeO )  THEN

      CALL get_H(i, ind_H, atypeH)
      DO iH=1,2
        CALL Distance_Vector_IND(ind_H(iH),i,d)
        dipole(:,i)=dipole(:,i)+d(1:3)
      END DO

      dipole(:,i) = dipole(:,i) -8*wannier(:,i)

    END IF 

  END DO !i 

  dipole=dipole*convert_to_debye

END SUBROUTINE COMPUTE_WATER_DIPOLE 

SUBROUTINE FIND_CENTER_CNT
  USE parameters, only : pos, box, natoms, atype
  USE wannier_mod, only : atypeC, center
  IMPLICIT NONE
  REAL*8 :: sum_coord(3),d(4)
  INTEGER :: i, nc

  center=0.d0
  sum_coord = 0.0
  nc=0
  DO i = 1,natoms
    IF (atype(i)==atypeC) THEN
      CALL DISTANCE_VECTOR(pos(:,i),center,d)
      sum_coord = sum_coord + d(1:3)
      nc = nc+1
    END IF
  END DO

  center = sum_coord/float(nc)

END SUBROUTINE FIND_CENTER_CNT

INTEGER FUNCTION Index_Histogram_Radius(pos)
  USE wannier_mod, ONLY : maxR,nhist_dipole,center,layers,nlayers
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: pos(3)
  REAL*8                     :: d(4),r
  INTEGER                    :: i,ind

  ind=0
  CALL DISTANCE_VECTOR(pos,center,d)
  r = sqrt(d(1)**2+d(2)**2)
  DO i=1,nlayers
    IF ((r>=layers(i)).AND.(r<layers(i+1))) THEN
        ind=i
    END IF
  END DO

  IF (ind==0) THEN
      print *, r, layers
      STOP
  END IF

  Index_Histogram_Radius = ind

END FUNCTION Index_Histogram_Radius

INTEGER FUNCTION Index_Histogram_Dipole(val)
  USE wannier_mod, ONLY : min_dipole,max_dipole,nhist_dipole
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: val
  REAL*8                     :: hrange,frac

  hrange=max_dipole-min_dipole
  frac=(val-min_dipole)/hrange*float(nhist_dipole)
  IF ((int(frac)+1>nhist_dipole).or.(int(frac)<0)) THEN
    PRINT *, 'Value ', val, ' outside histogram range'
    STOP
  END IF
  Index_Histogram_Dipole = int(frac) + 1

END FUNCTION Index_Histogram_Dipole

SUBROUTINE DISTRIBUTION_MOLECULAR_DIPOLE
  USE parameters, ONLY : natoms, atype, box, pos
  USE wannier_mod, ONLY : dipole, hist_dipole, atypeO, &
                       nhist_dipole, ntot, nlayers
  IMPLICIT NONE
  INTEGER                    :: iat, iR, idip, nwater, ipol
  INTEGER                    :: Index_Histogram_Radius
  INTEGER                    :: Index_Histogram_Dipole
  REAL*8                     :: dipole_moment, total_dipole(3,nlayers)

  total_dipole = 0.0
  DO iat=1,natoms

    IF (atype(iat)==atypeO) THEN

      iR = Index_Histogram_Radius(pos(:,iat)) 

      DO ipol=1,3
        idip  = Index_Histogram_Dipole(dipole(ipol,iat)) 
        hist_dipole(ipol,idip,iR)= hist_dipole(ipol,idip,iR)+1.0
      END DO
      
      dipole_moment = norm2(dipole(:,iat))
      idip  = Index_Histogram_Dipole(dipole_moment) 
      hist_dipole(4,idip,iR)=hist_dipole(4,idip,iR)+dipole_moment
      total_dipole(:,iR) = total_dipole(:,iR) + dipole(:,iat)
      ntot(iR)=ntot(iR)+1

    END IF

  END DO

  DO iR = 1,nlayers
    DO ipol = 1,3
      idip  = Index_Histogram_Dipole(total_dipole(iR,ipol)/box(ipol)) 
      hist_dipole(4+ipol,idip,iR)=hist_dipole(4+ipol,idip,iR)+1.0
    END DO
  END DO


END SUBROUTINE DISTRIBUTION_MOLECULAR_DIPOLE

SUBROUTINE get_H(ind_O, ind_H,atypeH)
  USE parameters, ONLY : atype, pos, cutoff_OH, natoms
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: atypeH,ind_O
  INTEGER                    :: i, k, ind_H(2)
  REAL*8                     :: Dist

  i=1;k=1

  ind_H = 0
  DO WHILE (( i <= natoms) .and. ( k <= 2))

    !Should be the index for H
    IF ( atype(i) == atypeH ) THEN

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

SUBROUTINE NORMALIZE_HISTOGRAM
  USE wannier_mod, ONLY : hist_dipole, nhist_dipole, ntot, nlayers
  IMPLICIT NONE
  INTEGER                    :: i,j,k

  DO i=1,7
    DO j=1,nlayers
      if (ntot(j)>100) then
        DO k=1,nhist_dipole
          hist_dipole(i,k,j)          = hist_dipole(i,k,j)/ntot(j)
        END DO
      else
        DO k=1,nhist_dipole
          hist_dipole(i,k,j)          = 0.d0 
        END DO
      end if
    END DO
  END DO

END SUBROUTINE NORMALIZE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE wannier_mod, ONLY : hist_dipole, maxR, &
                          nhist_dipole, nlayers, layers, &
                          max_dipole, min_dipole
  IMPLICIT NONE
  INTEGER :: ihist, ilayer
  REAL*8 :: hist_stride

  hist_stride = (max_dipole-min_dipole)/float(nhist_dipole)

  OPEN(unit = 2, file = "Dipole_distribution_Radial_distribution.dat")

  DO ilayer=1,nlayers

    WRITE (2,fmt="(A131,2(2X,F12.8))") '# Radius(A) Molecular_Dipole(X,Y,Z) &
                            & Molecular_Dipole(norm) Total_Dipole(X,Y,Z)/box(X,Y,Z) &
                            & Layer boundaries: ', layers(ilayer), layers(ilayer+1)

    DO ihist=1,nhist_dipole

      WRITE (2,fmt="(8(F12.8,3X))") hist_stride*(float(ihist)-0.5d0)+min_dipole, &
                                      hist_dipole(:,ihist,ilayer)
    END DO

  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
