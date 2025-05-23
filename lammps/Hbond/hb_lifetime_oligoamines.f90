MODULE histogram 
  IMPLICIT NONE 
  INTEGER*2, ALLOCATABLE       :: hb_matrix(:,:,:)
  INTEGER, ALLOCATABLE         :: NC_coordnumb(:), H_id(:)
  INTEGER                      :: ind_hb(2), nat1,nat2
  REAL*8, ALLOCATABLE          :: lifetime(:)
END MODULE histogram

PROGRAM ZDF
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL CLASSIFY_N_AND_H_TYPES
    CALL MAKE_HB_MATRIX(frame)
  END DO

  CALL COMPUTE_LIFETIME
  CALL PRINT_RESULTS

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist
  USE histogram, ONLY : hb_matrix, H_id, NC_coordnumb, &
                         lifetime, ind_hb, nat1, nat2
  IMPLICIT NONE
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, stride, nhist
  READ(*,*) ind_hb 
  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))

  CALL GET_NAT1_NAT2
  
  ALLOCATE(hb_matrix(nframes,nat1,nat2)); hb_matrix=0
  ALLOCATE(lifetime(nframes)); lifetime=0.d0
  ALLOCATE(H_id(natoms)); H_id=0
  ALLOCATE(NC_coordnumb(natoms)); NC_coordnumb=0

END SUBROUTINE INITIALIZE

SUBROUTINE GET_NAT1_NAT2
  USE parameters, only : natoms,atype
  USE histogram, only : nat1,nat2
  IMPLICIT NONE
  INTEGER :: i

  nat1=0
  CALL READ_ATOM_REDUCED
  do i=1,natoms
    if (atype(i)==2) nat1=nat1+1
  end do
  nat2=nat1
  REWIND(1)

END SUBROUTINE GET_NAT1_NAT2

SUBROUTINE MAKE_HB_MATRIX(frame)
  USE parameters, ONLY : pos, atype, natoms, nhist 
  USE histogram, ONLY : hb_matrix, NC_coordnumb, ind_hb 
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: frame
  INTEGER                    :: i, j, k, iat,jat
  INTEGER, DIMENSION(4)      :: ind_H 
  REAL*8                     :: d_OOw, Dist, Angle
  REAL*8                     :: ang_OHO 

  iat=1
  DO i = 1,natoms

    !Selecting only Ow
    IF ( ( ( atype(i) == 2 ).or.(atype(i) == 4)) .and. ( NC_coordnumb(i) == ind_hb(1)) )  THEN

      CALL get_H2(i, ind_H)

      jat=1
      DO j = 1,natoms

        !Selecting only Ow
        IF (((atype(j) == 2).or.(atype(i) == 4)).and.(NC_coordnumb(j) == ind_hb(2)).and.(i.ne.j)) THEN

          d_OOw = Dist(i, j)

          DO k = 1,4

            IF (ind_H(k) /= 0) THEN

              ang_OHO = Angle(pos(:,i), pos(:,j), pos(:,ind_H(k)))

              !Definition of H-bond - Luzar and Chandler
              !IF ( (d_OOw < 3.5) .and. (ang_OHO > 0.8660) ) THEN
              IF ( (d_OOw < 4.0) .and. (ang_OHO > 0.8660) ) THEN

                hb_matrix(frame,iat,jat) = 1
                hb_matrix(frame,jat,iat) = 1

              END IF !d_OOw

            END IF ! ind_H

          END DO ! k

          jat=jat+1

        END IF !ind_atom

      END DO !j

      iat=iat+1

    END IF !ind_atom

  END DO !i 

END SUBROUTINE MAKE_HB_MATRIX

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

SUBROUTINE get_H2(ind_O, ind_H)
  USE parameters, ONLY : atype, cutoff_OH, natoms
  IMPLICIT NONE
  INTEGER, DIMENSION(4)      :: ind_H
  INTEGER                    :: i, k, ind_O
  REAL*8                     :: Dist

  i=1;k=1
  ind_H=0

  DO WHILE (( i <= natoms) .and. ( k <= 4))

    !Should be the index for H
    IF ( atype(i) == 3 ) THEN

      IF ( Dist(ind_O,i) < cutoff_OH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

END SUBROUTINE get_H2

SUBROUTINE COMPUTE_LIFETIME
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : hb_matrix,lifetime,nat1,nat2
  USE parameters, ONLY : natoms, nframes
  IMPLICIT NONE
  INTEGER                    :: frame1, frame2, i, j
  INTEGER                    :: lt, c
  INTEGER, PARAMETER         :: nsep=100

  !frame1 is the duration of the interval 
  DO frame1 = 0,nframes/2

    lt=0
    c=0
!$omp parallel do reduction(+:lt)
    DO i = 1,nat1

      DO j = 1,nat2

        !Change initial condition, look for combinations with same interval
        DO frame2 = 1,nframes-frame1,nsep

          IF ( ALL( hb_matrix(frame2:frame2+frame1,i,j) == 1 ) ) THEN
         
            lt = lt + 1
          
          ENDIF

          c = c + hb_matrix(frame2,i,j)

        END DO

      END DO

    END DO
!$omp end parallel do

    lifetime(frame1+1)   = dble(lt) / dble(c)

  ENDDO

END SUBROUTINE COMPUTE_LIFETIME 

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nframes
  USE histogram, ONLY : lifetime 
  IMPLICIT NONE
  INTEGER :: n

  OPEN(unit = 2, file = "HB_lifetime.dat")

  WRITE(2, *) "# MD_Step  Lifetime_probability"

  DO n=1,nframes/2
      WRITE(2,fmt='(I10,1X,F12.8)') n, lifetime(n)
  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
