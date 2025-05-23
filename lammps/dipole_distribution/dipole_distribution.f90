MODULE wannier_mod 
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE         :: ntot(:)
  INTEGER, PARAMETER           :: nhist_dipole=2000
  INTEGER                      :: atypeO, atypeH, atypeC
  REAL*8, ALLOCATABLE          :: dipole(:,:) 
  REAL*8, ALLOCATABLE          :: ionic_dipole(:,:), electric_dipole(:,:)
  REAL*8                       :: T, center(3), hist_dipole(7,nhist_dipole)
  REAL*8                       :: hist_electric_dipole(4,nhist_dipole)
  REAL*8                       :: hist_ionic_dipole(4,nhist_dipole)
  REAL*8, PARAMETER            :: min_dipole=-20.d0, max_dipole=20.d0
  REAL*8, PARAMETER            :: convert_to_debye=1.d0/0.2081943d0
END MODULE wannier_mod

PROGRAM Dipole_Distribution
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL READ_WANNIER
    CALL COMPUTE_WATER_DIPOLE
    CALL DISTRIBUTION_MOLECULAR_DIPOLE
    CALL DISTRIBUTION_IONIC_ELECTRIC_DIPOLE
    !CALL RADIAL_DISTRIBUTION_DIPOLE
  END DO

  CALL NORMALIZE_HISTOGRAM
  CALL PRINT_RESULTS

END PROGRAM Dipole_Distribution

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist, nwann, wannier
  USE wannier_mod, ONLY : dipole,atypeO, atypeH, atypeC, T, &
                          ionic_dipole, electric_dipole, &
                          hist_dipole, hist_electric_dipole, &
                          hist_ionic_dipole
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, wfc_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, wfc_file)
  CALL getarg(3, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms,nwann, nframes, T!, nhist, rmax
  READ(1,*) atypeO, atypeH, atypeC
  !READ(1, *) center 
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(dipole(3,natoms)); dipole=0.d0
  ALLOCATE(ionic_dipole(3,natoms)); ionic_dipole=0.d0
  ALLOCATE(electric_dipole(3,natoms)); electric_dipole=0.d0
  ALLOCATE(wannier(3,natoms)); wannier=0.d0

  hist_dipole=0
  hist_electric_dipole=0
  hist_ionic_dipole=0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
  OPEN(unit = 3, file =  wfc_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_WATER_DIPOLE
  USE parameters, ONLY : pos, atype, natoms, wannier
  USE wannier_mod, ONLY : dipole, convert_to_debye, atypeO, &
                          atypeH, electric_dipole, ionic_dipole
  IMPLICIT NONE
  INTEGER                    :: i, iH, ind_H(2)
  REAL*8                     :: d(4)

  dipole=0.d0
  ionic_dipole=0.d0
  electric_dipole=0.d0

  DO i = 1,natoms

    !Selecting only Ow
    IF ( atype(i) == atypeO )  THEN

      CALL get_H(i, ind_H, atypeH)
      DO iH=1,2
        CALL Distance_Vector_IND(ind_H(iH),i,d)
        ionic_dipole(:,i)=ionic_dipole(:,i)+d(1:3)
      END DO

      electric_dipole(:,i) = electric_dipole(:,i) -8*wannier(:,i)

    END IF 

  END DO !i 

  electric_dipole=electric_dipole*convert_to_debye
  ionic_dipole=ionic_dipole*convert_to_debye
  dipole=electric_dipole+ionic_dipole

END SUBROUTINE COMPUTE_WATER_DIPOLE 

INTEGER FUNCTION Index_Histogram(val)
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
  Index_Histogram = int(frac) + 1

END FUNCTION Index_Histogram

SUBROUTINE DISTRIBUTION_MOLECULAR_DIPOLE
  USE parameters, ONLY : natoms, atype, box
  USE wannier_mod, ONLY : dipole, hist_dipole, atypeO, &
                      max_dipole,min_dipole, nhist_dipole
  IMPLICIT NONE
  INTEGER                    :: iat, ihist, nwater, ipol
  INTEGER                    :: Index_Histogram
  REAL*8                     :: total_dipole(3), dipole_range
  REAL*8                     :: dipole_moment

  total_dipole=0.d0
  dipole_range=max_dipole-min_dipole

  DO iat=1,natoms

    IF (atype(iat)==atypeO) THEN

      DO ipol=1,3
        ihist = Index_Histogram(dipole(ipol,iat)) 
        hist_dipole(ipol,ihist)= hist_dipole(ipol,ihist)+1.d0
      END DO
      
      dipole_moment = norm2(dipole(:,iat))
      ihist = Index_Histogram(dipole_moment) 
      hist_dipole(4,ihist)=hist_dipole(4,ihist)+1.d0

      total_dipole(:)=total_dipole(:)+dipole(:,iat)

    END IF

  END DO

  DO ipol=1,3
    ihist = Index_Histogram(total_dipole(ipol)/box(ipol)) 
    hist_dipole(ipol+4,ihist)= hist_dipole(ipol+4,ihist)+1.d0
  END DO

END SUBROUTINE DISTRIBUTION_MOLECULAR_DIPOLE

SUBROUTINE DISTRIBUTION_IONIC_ELECTRIC_DIPOLE
  USE parameters, ONLY : natoms, atype, box
  USE wannier_mod, ONLY : electric_dipole, hist_electric_dipole, &
                          ionic_dipole, hist_ionic_dipole, atypeO, &
                          max_dipole, min_dipole, nhist_dipole, dipole
  IMPLICIT NONE
  INTEGER                    :: iat, ihist, nwater, ipol
  INTEGER                    :: Index_Histogram
  REAL*8                     :: dipole_range
  REAL*8                     :: dipole_moment

  dipole_range=max_dipole-min_dipole

  DO iat=1,natoms
    !
    IF (atype(iat)==atypeO) THEN
      !
      DO ipol=1,3
        ihist = Index_Histogram(electric_dipole(ipol,iat)) 
        hist_electric_dipole(ipol,ihist)= hist_electric_dipole(ipol,ihist)+1.d0
        !
        ihist = Index_Histogram(ionic_dipole(ipol,iat)) 
        hist_ionic_dipole(ipol,ihist)= hist_ionic_dipole(ipol,ihist)+1.d0
      END DO
      !
      dipole_moment =sum(electric_dipole(:,iat)*dipole(:,iat)) / &
                        norm2(dipole(:,iat))
      ihist = Index_Histogram(dipole_moment) 
      hist_electric_dipole(4,ihist)=hist_electric_dipole(4,ihist)+1.d0
      !
      dipole_moment =sum(ionic_dipole(:,iat)*dipole(:,iat)) / &
                        norm2(dipole(:,iat))
      ihist = Index_Histogram(dipole_moment) 
      hist_ionic_dipole(4,ihist)=hist_ionic_dipole(4,ihist)+1.d0
      !
    END IF
  !
  END DO

END SUBROUTINE DISTRIBUTION_IONIC_ELECTRIC_DIPOLE

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
  USE wannier_mod, ONLY : hist_dipole, nhist_dipole, &
                          hist_electric_dipole, hist_ionic_dipole
  IMPLICIT NONE
  REAL*8                     :: sum_hist, sum_el_hist, sum_ion_hist
  INTEGER                    :: i,j

  DO i = 1,7
    sum_hist=0.d0
    DO j=1,nhist_dipole
      sum_hist=sum_hist+hist_dipole(i,j)
    END DO
    DO j=1,nhist_dipole
      hist_dipole(i,j)=hist_dipole(i,j)/sum_hist
    END DO
  END DO

  DO i = 1,4
    sum_ion_hist=0.d0
    sum_el_hist=0.d0
    DO j=1,nhist_dipole
      sum_ion_hist=sum_ion_hist+hist_ionic_dipole(i,j)
      sum_el_hist=sum_el_hist+hist_electric_dipole(i,j)
    END DO
    DO j=1,nhist_dipole
      hist_ionic_dipole(i,j)=hist_ionic_dipole(i,j)/sum_ion_hist
      hist_electric_dipole(i,j)=hist_electric_dipole(i,j)/sum_el_hist
    END DO
  END DO

END SUBROUTINE NORMALIZE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist
  USE wannier_mod, ONLY : hist_dipole, min_dipole, max_dipole, &
                          nhist_dipole, hist_electric_dipole, &
                          hist_ionic_dipole
  IMPLICIT NONE
  INTEGER :: ihist
  REAL*8 :: hist_stride

  hist_stride = (max_dipole-min_dipole)/nhist_dipole

  OPEN(unit = 2, file = "Dipole_distribution.dat")
  OPEN(unit = 5, file = "Electric_dipole_distribution.dat")
  OPEN(unit = 6, file = "Ionic_dipole_distribution.dat")

  WRITE (2,fmt="(A111)") '# Dipole(Debye) Hist_Water_Molecular_Dipole(X,Y,X) &
                          & Hist_Water_Molecular_Dipole(norm) Total_Water_Dipole(X,Y,Z)' 
  WRITE (5,fmt="(A111)") '# Dipole(Debye) Hist_Water_Electric_Dipole(X,Y,X) &
                          & Hist_Water_Electric_Dipole(norm)' 
  WRITE (6,fmt="(A111)") '# Dipole(Debye) Hist_Water_Ionic_Dipole(X,Y,X) &
                          & Hist_Water_Ionic_Dipole(norm) ' 

  DO ihist=1,nhist_dipole

    WRITE (2,fmt="(8(F12.8,3X))") hist_stride*(float(ihist)-0.5d0)+min_dipole, &
                                    hist_dipole(:,ihist)
    WRITE (5,fmt="(5(F12.8,3X))") hist_stride*(float(ihist)-0.5d0)+min_dipole, &
                                    hist_electric_dipole(:,ihist)
    WRITE (6,fmt="(5(F12.8,3X))") hist_stride*(float(ihist)-0.5d0)+min_dipole, &
                                    hist_ionic_dipole(:,ihist)

  END DO

  CLOSE(2)
  CLOSE(5)
  CLOSE(6)

END SUBROUTINE PRINT_RESULTS
