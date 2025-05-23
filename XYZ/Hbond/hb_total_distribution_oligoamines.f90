MODULE histogram 
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE         :: NC_coordnumb(:), H_id(:)
  INTEGER                      :: hist_hb_peratom(2,10)
  REAL*8                       :: hb_length(100)
  CHARACTER(5)                 :: ind_hb(2)
END MODULE histogram

PROGRAM ZDF
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_EXTXYZ
    if (MOD(frame,stride)==0) THEN
      CALL ANALYSIS
    end if
  END DO

  CALL PRINT_RESULTS

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist 
  USE histogram, ONLY : ind_hb
  IMPLICIT NONE
  INTEGER                    :: i,j
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, stride
  READ(*,*) ind_hb
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS
  USE parameters, ONLY : pos, atype, natoms, nhist 
  USE histogram, ONLY : hist_hb_peratom, ind_hb 
  IMPLICIT NONE
  INTEGER                    :: i, j, k
  INTEGER                    :: nhb_peratom(2,natoms)
  INTEGER                    :: intra_or_inter, hbtype
  INTEGER, DIMENSION(4)      :: ind_H 
  REAL*8                     :: d_OOw, Dist, Angle, n(nhist)
  REAL*8                     :: ang_OHO
  LOGICAL                    :: is_hb_type

  nhb_peratom=0

  DO i = 1,natoms

    IF ( is_hb_type(i) ) THEN

      ind_H = 0
      CALL get_H(i, ind_H)

      DO j = 1,natoms

        IF ( ( is_hb_type(j) ) .and. (i.ne.j) ) THEN

          d_OOw = Dist(i, j)

          DO k = 1,4

            IF (ind_H(k) /= 0) THEN

              ang_OHO = Angle(pos(:,i), pos(:,j), pos(:,ind_H(k)))

              !Definition of H-bond - Luzar and Chandler
              IF ( (d_OOw < 3.5) .and. (ang_OHO > 0.8660) ) THEN

                IF (trim(atype(i))==trim(ind_hb(1))) THEN
                  nhb_peratom(1,i) = nhb_peratom(1,i) + 1
                  CALL UPDATE_HBLEN_HISTOGRAM(j,ind_H(k))
                END IF
                IF (trim(atype(j))==trim(ind_hb(1))) THEN
                  nhb_peratom(2,j) = nhb_peratom(2,j) + 1
                END IF

              END IF !d_OOw

            END IF ! ind_H

          END DO ! k

        END IF ! atype
 
      END DO !j

    END IF !ind_atom

  END DO !i 

  DO i=1,natoms
    IF ( trim(atype(i)) == trim(ind_hb(1)) ) THEN
      CALL ADD_ONE(hist_hb_peratom(1,nhb_peratom(1,i)+1))
      CALL ADD_ONE(hist_hb_peratom(2,nhb_peratom(2,i)+1))
    END IF
  END DO

END SUBROUTINE ANALYSIS

SUBROUTINE ADD_ONE(hist)
  IMPLICIT NONE
  INTEGER :: hist
  hist=hist+1
END SUBROUTINE ADD_ONE

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : atype, pos, cutoff_NH, natoms
  IMPLICIT NONE
  INTEGER, DIMENSION(4)      :: ind_H
  INTEGER                    :: i, k, ind_O
  REAL*8                     :: Dist

  i=1;k=1

  DO WHILE (( i <= natoms) .and. ( k <= 4))

    !Should be the index for H
    IF ( trim(atype(i)) == "H" ) THEN

      IF ( Dist(ind_O,i) < cutoff_NH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

END SUBROUTINE get_H

SUBROUTINE UPDATE_HBLEN_HISTOGRAM(Acceptor, H)
  USE histogram, only : hb_length
  IMPLICIT NONE
  INTEGER, INTENT(IN)         :: Acceptor, H
  REAL*8                      :: Dist, d1
  INTEGER                     :: ind

  d1 = Dist(Acceptor, H)
  ind = int((d1-1.d0)/2.d0*100)+1
  hb_length(ind) = hb_length(ind) + 1

END SUBROUTINE UPDATE_HBLEN_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram, ONLY : hist_hb_peratom, hb_length, ind_hb
  IMPLICIT NONE
  INTEGER :: ihist

  OPEN(unit = 2, file = "HB_Number_"//trim(ind_hb(1))//".dat")
  OPEN(unit = 3, file = "HB_Length_"//trim(ind_hb(1))//".dat")

  WRITE(2, *) "# N_Hbonds Count(Donated) Count(Accepted)"
  WRITE(3, *) "# HBond_Length_Donated(A) Count"

  DO ihist=1,10
    WRITE(2,fmt='(I3,3X,2(I10,3X))') ihist-1, &
            hist_hb_peratom(1,ihist), &
            hist_hb_peratom(2,ihist)
  END DO

  DO ihist=1,100
    WRITE(3,fmt="(F12.8,3X,F12.3)") float(ihist-1)/50.+1.d0, &
           hb_length(ihist) 
  END DO

  CLOSE(2);CLOSE(3)

END SUBROUTINE PRINT_RESULTS

LOGICAL FUNCTION is_hb_type(ind)
  USE parameters, only : atype
  USE histogram, only : ind_hb
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind
  is_hb_type=.false.
  IF (trim(atype(ind))==trim(ind_hb(1))) THEN
      is_hb_type=.true.
  ELSE IF (trim(atype(ind))==trim(ind_hb(2))) THEN
      is_hb_type=.true.
  END IF
END FUNCTION is_hb_type
