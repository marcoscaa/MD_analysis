MODULE local 
  IMPLICIT NONE 
  LOGICAL, ALLOCATABLE               :: is_water_ion(:)
  INTEGER                            :: nhb_tot, nhb_acc, nhb_don, ntot 
  INTEGER, ALLOCATABLE               :: ind_Ow(:)
  INTEGER, PARAMETER                 :: nhbhist=100
  REAL*8                             :: hb_length_acc(nhbhist)
  REAL*8                             :: hb_length_don(nhbhist)
  REAL*8, PARAMETER                  :: rmin=0.9, rmax=3.0
END MODULE local 

PROGRAM Hbond 
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    if (MOD(frame,stride)==0) THEN
      CALL FIND_WATER_ION
      CALL ANALYSIS 
    end if
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, atype, atypeO, &
                         nframes, nequil, atypeH, &
                         atype, stride
  USE local,      ONLY : is_water_ion, nhb_tot, nhb_acc, ntot, &
                         nhb_don, hb_length_acc, hb_length_don, &
                         ind_Ow
  IMPLICIT NONE
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil,stride
  READ(1,*) atypeO, atypeH 

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(ind_Ow(natoms))
  ALLOCATE(is_water_ion(natoms))

  nhb_tot=0
  nhb_acc=0
  nhb_don=0
  hb_length_acc=0.d0
  hb_length_don=0.d0
  ntot=0

  CLOSE(1)
  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS
  USE parameters, ONLY : pos, atype, natoms, atypeO, atypeH 
  USE local,      ONLY : nhb_acc, nhb_don, nhb_tot, ntot, &
                         hb_length_acc, hb_length_don, is_water_ion, &
                         ind_Ow
  IMPLICIT NONE
  INTEGER                    :: i, j, iH
  REAL*8                     :: d_OOw, Dist,  Angle
  REAL*8                     :: ang_OHO, d1, d2 

  DO i = 1,natoms
    !Selecting only Ow from water ion
    IF ( is_water_ion(i) )  THEN

      ntot = ntot + 1
      DO j = 1,natoms

        !Selecting only Ow
        IF ((atype(j) == atypeO).and.(i.ne.j)) THEN

          d_OOw = Dist(i, j)

          !Definition of H-bond - Luzar and Chandler
          IF (d_OOw<3.5) THEN

            DO iH=1,natoms

              IF (atype(iH) == atypeH) THEN

                IF ((ind_Ow(iH)==i).or.(ind_Ow(iH)==j)) THEN

                  d1 = Dist(iH,i)
                  d2 = Dist(iH,j)
                  ang_OHO = Angle(pos(:,iH), pos(:,i), pos(:,j))

                  !Definition of H-bond - Luzar and Chandler
                  IF (ang_OHO < -0.8660) THEN

                    IF (d1<d2) THEN
                      nhb_don = nhb_don + 1 
                      CALL ADD_ONE(hb_length_don,d2)
                    ELSE
                      nhb_acc = nhb_acc + 1 
                      CALL ADD_ONE(hb_length_acc,d1)
                    END IF

                    nhb_tot = nhb_tot + 1

                  END IF !and_OHO

                END IF !d1,d2

              END IF !atype(iH)

            END DO !iH

          END IF !dOOw

        END IF !atype

      END DO !j

    END IF !is_water_ion

  END DO !i 

END SUBROUTINE ANALYSIS

SUBROUTINE ADD_ONE(hist,val)
  USE local, only : nhbhist, rmax, rmin
  IMPLICIT NONE
  INTEGER                    :: ihb
  REAL*8                     :: hist(nhbhist)
  REAL*8, INTENT(IN)         :: val

  ihb = int( ( val - rmin )*nhbhist / ( rmax - rmin ) ) + 1
  hist(ihb) = hist(ihb) + 1.0

END SUBROUTINE ADD_ONE

SUBROUTINE FIND_WATER_ION
  USE parameters, ONLY : atype, natoms, atypeO, atypeH
  USE local, ONLY : is_water_ion, ind_Ow
  IMPLICIT NONE
  INTEGER :: iat, cn_ow(natoms)
  INTEGER :: closest_atom

  cn_ow=0
  ind_Ow=0
  is_water_ion=.false.
  DO iat=1,natoms
    if (atype(iat)==atypeH) then !2: H index
      ind_Ow(iat) = closest_atom(iat,atypeO)
      cn_ow(ind_Ow(iat)) = cn_ow(ind_Ow(iat)) + 1
    END IF
  END DO
      
  DO iat=1,natoms
    if (atype(iat)==atypeO) then !1: O index
      if(cn_Ow(iat)/=2) then
        is_water_ion(iat)=.True.
        !DO jat=1,natoms
        !if ((atype(jat)==atypeH).and.(ind_Ow(jat)==iat)) then
        !    is_water_ion(jat)=.True.
        !  end if
        !END DO
      end if
    end if
  END DO

END SUBROUTINE FIND_WATER_ION

SUBROUTINE PRINT_RESULTS
  USE local, ONLY : ntot, nhb_acc, nhb_don, hb_length_acc, &
                    hb_length_don, nhb_tot, rmin, rmax, nhbhist 
  IMPLICIT NONE
  REAL*8                     :: mean_hb_len_don, mean_hb_len_acc
  REAL*8                     :: compute_mean 
  INTEGER                    :: ihb
  
  OPEN(unit = 2,file = "Hbonds.dat")
  OPEN(unit = 3,file = "Hbond_length_distribution.dat")

  write(2,*) '#', 'n_hb_donated', 'n_hb_accepted', 'n_hb_tot', &
             'mean_hb_length_donated', 'mean_hb_length_accepted'

  mean_hb_len_don=compute_mean(hb_length_don)
  mean_hb_len_acc=compute_mean(hb_length_acc)

  WRITE(2,fmt="(5(F12.8,3X))") &
             float(nhb_don)/float(ntot), float(nhb_acc)/float(ntot), &
             float(nhb_tot)/float(ntot), mean_hb_len_don, &
             mean_hb_len_acc 

  CLOSE(2)

  write(3,*) "#", "Hb_Length(A)", "P(Hb_donated)", "P(Hb_accepted)"

  DO ihb = 1, nhbhist
    write(3,fmt="(3(F12.8,3X))") &
             float(ihb-1)/float(nhbhist)*(rmax-rmin)+rmin, &
             hb_length_don(ihb)/float(ntot), &
             hb_length_acc(ihb)/float(ntot)
  END DO

  CLOSE(3)

END SUBROUTINE

INTEGER FUNCTION closest_atom(icenter,atyp)
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: icenter, atyp
  INTEGER :: iat, iat_min
  REAL*8 :: d, min_d, Dist

  min_d = 100.
  iat_min = 0

  DO iat=1,natoms

    IF (atype(iat)==atyp) THEN
      d = Dist(iat,icenter)
      IF (d < min_d) THEN
        min_d=d
        iat_min=iat
      END IF
    END IF

  END DO

  IF (iat_min==0) THEN
    PRINT *, "Could not find closest atom to atom", icenter
    PRINT *, "STOP!!!!"
    STOP
  END IF

  closest_atom = iat_min

END FUNCTION closest_atom

REAL*8 FUNCTION compute_mean(array)
  USE local, only : nhbhist, rmin, rmax
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: array(nhbhist)
  REAL*8                     :: mean, total
  INTEGER                    :: i

  mean = 0.d0
  total = 0.d0
  DO i=1,nhbhist
    total = total + array(i)
    mean = mean + array(i)*(float(i-1)*(rmax-rmin)/float(nhbhist)+rmin)
  END DO

  compute_mean = mean / total 

END FUNCTION compute_mean
