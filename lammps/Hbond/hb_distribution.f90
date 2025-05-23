MODULE histogram 
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE         :: nhb_acc(:), nhb_don(:), ntot(:)
  INTEGER                      :: hist_hb_peratom(3,10)
  REAL*8                       :: rmax, center(3)
  REAL*8, ALLOCATABLE          :: hb_length(:)
END MODULE histogram

PROGRAM ZDF
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

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist, typeO, typeH
  USE histogram, ONLY : nhb_acc, nhb_don, ntot, rmax, center, &
                        hb_length
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, stride, nhist, rmax
  READ(1, *) typeO, typeH 
  READ(1, *) center 
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(nhb_acc(nhist)); nhb_acc=0
  ALLOCATE(nhb_don(nhist)); nhb_don=0
  ALLOCATE(hb_length(nhist)); hb_length=0.d0
  ALLOCATE(ntot(nhist)); ntot=0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS
  USE parameters, ONLY : pos, atype, natoms, nhist, typeO, typeH 
  USE histogram, ONLY : nhb_acc, nhb_don, ntot, rmax, center, hist_hb_peratom,&
                        hb_length
  IMPLICIT NONE
  INTEGER                    :: i, j, indri, indrj, iH
  INTEGER                    :: nhb_peratom(2,natoms)
  INTEGER, DIMENSION(3)      :: ind_H 
  REAL*8                     :: d_OOw, Dist, Dist_Center, Angle, n(nhist)
  REAL*8                     :: ri, rj, ang_OHO, d 

  nhb_peratom=0

  DO i = 1,natoms

    !Selecting only Ow
    IF ( atype(i) == typeO )  THEN

      ri = Dist_Center(center,pos(:,i))
      indri = int(nhist*ri/rmax)+1
      ntot(indri) = ntot(indri)+1
      ind_H = 0
      CALL get_H2(i, ind_H)

      DO iH=1,3
      
        IF (ind_H(iH) /= 0) THEN

          DO j = 1,natoms

            !Selecting only Ow
            IF ((atype(j) == typeO).and.(i.ne.j)) THEN

              d_OOw = Dist(i, j)
              rj = Dist_Center(center,pos(:,j))
              indrj = int(nhist*rj/rmax)+1

              ang_OHO = Angle(pos(:,i), pos(:,j), pos(:,ind_H(iH)))

              !Definition of H-bond - Luzar and Chandler
              IF ( (d_OOw < 3.5) .and. (ang_OHO > 0.8660) ) THEN
              !IF ( (d_OOw < 3.5) .and. (maxval(ang_OHO) > 0.6427) ) THEN

                nhb_acc(indrj) = nhb_acc(indrj) + 1 
                nhb_don(indri) = nhb_don(indri) + 1 
                d = Dist(ind_H(iH),j)
                hb_length(indri) = hb_length(indri) + d
                hb_length(indrj) = hb_length(indrj) + d
                nhb_peratom(1,i) = nhb_peratom(1,i) + 1
                nhb_peratom(2,j) = nhb_peratom(2,j) + 1

              END IF !d_OOw

            END IF !ind_atom

          END DO !j

        END IF !ind_H

      END DO

    END IF !ind_atom

  END DO !i 

  DO i=1,natoms
  IF ((atype(i)==typeO)) THEN
      CALL ADD_ONE(hist_hb_peratom(1,nhb_peratom(1,i)+1))
      CALL ADD_ONE(hist_hb_peratom(2,nhb_peratom(2,i)+1))
      CALL ADD_ONE(hist_hb_peratom(3,sum(nhb_peratom(:,i))+1))
    ENDIF
  END DO

END SUBROUTINE ANALYSIS

SUBROUTINE ADD_ONE(hist)
  IMPLICIT NONE
  INTEGER :: hist
  hist=hist+1
END SUBROUTINE ADD_ONE

SUBROUTINE get_H2(ind_O, ind_H)
  USE parameters, ONLY : atype, typeH, pos, cutoff_OH, natoms
  IMPLICIT NONE
  INTEGER, DIMENSION(2)      :: ind_H
  INTEGER                    :: i, k, ind_O
  REAL*8                     :: Dist

  i=1;k=1

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

  !IF (k.ne.3) THEN
  !  print *, 'Found water ion with charge ', k-3
  !END IF
END SUBROUTINE get_H2

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist
  USE histogram, ONLY : rmax,nhb_acc,nhb_don,ntot,hist_hb_peratom, &
                        hb_length
  IMPLICIT NONE
  INTEGER :: ihist, nhb

  OPEN(unit = 2, file = "HB_distribution.dat")

  WRITE(2, *) "# R HB_acc HB_don HB_length N_particles"

  DO ihist=1,nhist
    nhb=nhb_acc(ihist)+nhb_don(ihist)
    if (nhb==0) nhb=1
    if (ntot(ihist).ne.0) then
      WRITE(2,fmt='(F12.8,3(3X,F12.8),3X,I10)') (float(ihist-1))*rmax/float(nhist), &
              float(nhb_acc(ihist))/float(ntot(ihist)), &
              float(nhb_don(ihist))/float(ntot(ihist)), &
              hb_length(ihist)/float(nhb), &
              ntot(ihist)
    else
      WRITE(2,fmt='(F12.8,3(3X,F12.8),3X,I10)') (float(ihist-1))*rmax/float(nhist),0.0,0.0,0.0,0
    end if
  END DO

  CLOSE(2)

  OPEN(unit = 2, file = "HB_Number.dat")

  WRITE(2, *) "# N_Hbonds Count(Donated) Count(Accepted) Count(Total)"

  DO ihist=1,10
    WRITE(2,fmt='(I3,3X,3(I10,3X))') ihist-1, &
            hist_hb_peratom(1,ihist), &
            hist_hb_peratom(2,ihist), &
            hist_hb_peratom(3,ihist)
  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
