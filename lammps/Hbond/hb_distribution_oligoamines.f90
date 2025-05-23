MODULE histogram 
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE         :: NC_coordnumb(:), H_id(:)
  INTEGER                      :: hist_hb_peratom(6,10), ind_hb(2)
END MODULE histogram

PROGRAM ZDF
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    if (frame==1) CALL DEFINE_MOLECULE 
    if (MOD(frame,stride)==0) THEN
      CALL CLASSIFY_N_AND_H_TYPES
      CALL ANALYSIS
    end if
  END DO

  CALL PRINT_RESULTS

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist, rcut, molid, &
                         index_bonded, nattype
  USE histogram, ONLY :  NC_coordnumb, H_id, ind_hb
  IMPLICIT NONE
  INTEGER                    :: i,j
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, stride, nattype
  ALLOCATE(rcut(nattype,nattype))
  READ(1, *) ( ( rcut(i,j), j=i,nattype ), i=1,nattype )
  READ(*,*) ind_hb
  do i=1,nattype
    do j=i+1,nattype
      rcut(j,i)=rcut(i,j)
    end do
  end do
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(H_id(natoms)); H_id=0
  ALLOCATE(NC_coordnumb(natoms)); NC_coordnumb=0
  ALLOCATE(molid(natoms))
  ALLOCATE(index_bonded(natoms,10))

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE CLASSIFY_N_AND_H_TYPES
  USE parameters, only : natoms, atype
  USE histogram, only : NC_coordnumb, H_id
  IMPLICIT NONE
  INTEGER :: iat, jat
  REAL*8 :: Dist, d
  REAL*8, parameter :: rcut_H = 1.5d0, rcut_CN=1.8d0
  
  H_id=0
  NC_coordnumb=0
  do iat=1,natoms
    IF (atype(iat)==2) then
      do jat=1,natoms
        if (atype(jat)==1) then 
          d=Dist(iat,jat)
          if (d<rcut_CN) NC_coordnumb(iat) = NC_coordnumb(iat) + 1
        end if
      end do
    ELSE IF (atype(iat)==3) then
      do jat=1,natoms
        if ((atype(jat)==2).or.(atype(jat)==4)) then 
          d=Dist(iat,jat)
          !if (d<rcut_H) H_id(iat) = atype(jat)
          if (d<rcut_H) H_id(iat) = jat
        end if
      end do
    END IF
  end do

END SUBROUTINE CLASSIFY_N_AND_H_TYPES

SUBROUTINE ANALYSIS
  USE parameters, ONLY : pos, atype, natoms, nhist 
  USE histogram, ONLY : hist_hb_peratom, &
                        NC_coordnumb, H_id, ind_hb
  IMPLICIT NONE
  INTEGER                    :: i, j, k
  INTEGER                    :: nhb_peratom(2,2,natoms)
  INTEGER                    :: intra_or_inter, hbtype
  INTEGER, DIMENSION(4)      :: ind_H 
  REAL*8                     :: d_OOw, Dist, Angle, n(nhist)
  REAL*8                     :: ang_OHO

  nhb_peratom=0

  DO i = 1,natoms

    !Selecting only N with CN coordination number ind_hb(1)
    IF ( (( atype(i) == 2 ).or.(atype(i) == 4)) &
      .and. ( NC_coordnumb(i) == ind_hb(1)) )  THEN

      ind_H = 0
      CALL get_H(i, ind_H)

      DO j = 1,natoms

        !Selecting only N and O atoms 
        IF (((atype(j) == 2).or.(atype(j)==4)).and.(i.ne.j)) THEN

          d_OOw = Dist(i, j)
          hbtype=intra_or_inter(i,j)

          DO k = 1,4

            IF (ind_H(k) /= 0) THEN

              ang_OHO = Angle(pos(:,i), pos(:,j), pos(:,ind_H(k)))

              !Definition of H-bond - Luzar and Chandler
              !IF ( (d_OOw < 3.5) .and. (ang_OHO > 0.8660) ) THEN
              IF ( (d_OOw < 4.0) .and. (ang_OHO > 0.8660) ) THEN

                nhb_peratom(1,hbtype,i) = nhb_peratom(1,hbtype,i) + 1
                nhb_peratom(2,hbtype,j) = nhb_peratom(2,hbtype,j) + 1

              END IF !d_OOw

            END IF ! ind_H

          END DO ! k

        END IF ! atype
 
      END DO !j

    END IF !ind_atom

  END DO !i 

  DO i=1,natoms
    IF ((atype(i)==2).or.(atype(i)==4)) THEN
      IF (NC_coordnumb(i) == ind_hb(1)) THEN
        CALL ADD_ONE(hist_hb_peratom(1,nhb_peratom(1,1,i)+1))
        CALL ADD_ONE(hist_hb_peratom(2,nhb_peratom(2,1,i)+1))
        CALL ADD_ONE(hist_hb_peratom(3,nhb_peratom(1,2,i)+1))
        CALL ADD_ONE(hist_hb_peratom(4,nhb_peratom(2,2,i)+1))
        CALL ADD_ONE(hist_hb_peratom(5,sum(nhb_peratom(1,:,i))+1))
        CALL ADD_ONE(hist_hb_peratom(6,sum(nhb_peratom(2,:,i))+1))
      END IF
    ENDIF
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
    IF ( atype(i) == 3 ) THEN

      IF ( Dist(ind_O,i) < cutoff_NH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

END SUBROUTINE get_H

SUBROUTINE PRINT_RESULTS
  USE histogram, ONLY : hist_hb_peratom
  IMPLICIT NONE
  INTEGER :: ihist

  OPEN(unit = 2, file = "HB_Number.dat")

  WRITE(2, *) "# N_Hbonds Count(Don, Intra) Count(Acc, Intra) Count(Don, Inter) Count(Acc, Inter) Count(Total Don) Count(Total Acc)"

  DO ihist=1,10
    WRITE(2,fmt='(I3,3X,6(I10,3X))') ihist-1, &
            hist_hb_peratom(1,ihist), &
            hist_hb_peratom(2,ihist), &
            hist_hb_peratom(3,ihist), &
            hist_hb_peratom(4,ihist), &
            hist_hb_peratom(5,ihist), &
            hist_hb_peratom(6,ihist)
  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

INTEGER FUNCTION intra_or_inter(i,j)
  USE parameters, only : molid
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  if (molid(i)==molid(j)) then
    intra_or_inter=1
  else
    intra_or_inter=2
  end if
END FUNCTION intra_or_inter
