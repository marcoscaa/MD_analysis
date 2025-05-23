!Set of subroutines used to analyse the lipid trajectory

MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: rotcorr(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: vect(:,:,:), rotfreq(:,:)
  INTEGER                                  :: stride, tcorr, tlinear
  INTEGER                                  :: typerot1,typerot2
  INTEGER                                  :: nattype1,nattype2
  INTEGER, ALLOCATABLE                     :: is_bound(:,:,:), rindex(:,:,:)
  INTEGER                                  :: nhist1,nhist2, central_atom(3)

END MODULE histogram

PROGRAM Hbond
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL ASSIGN_VECTOR (frame)
  END DO

  CALL ROT_CORR 
  CALL COMPUTE_ROTATION_FREQUENCY
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : vect, rotcorr, stride, tcorr, is_bound, &
                         nattype1, nattype2, typerot1, typerot2, &
                         nhist1,nhist2, central_atom, rindex, rotfreq, &
                         tlinear
  USE parameters, ONLY : natoms, nframes, nattype,&
                         pos, nequil, atype, dt
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nattype1, nattype2, nframes, nequil, dt, stride
  READ(1,*) typerot1, typerot2, tcorr, tlinear, nhist1,nhist2 
  READ(*,*) central_atom
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(vect(3,nframes,nattype2)); vect=0 ! normal bond vector
  ALLOCATE(is_bound(nframes,nattype1,nattype2)); is_bound=0
  ALLOCATE(rindex(2,nattype1,nframes)); rindex=0
  ALLOCATE(rotcorr(tcorr,nhist1,nhist2)); rotcorr = 0
  ALLOCATE(rotfreq(nhist1,nhist2)); rotfreq=0.d0
  ALLOCATE(atype(natoms))

  OPEN(unit = 1,file = file_pos)
  
END SUBROUTINE INITIALIZE

SUBROUTINE ROT_CORR
  ! Calculate the velocity msdelation function of water atoms 
  ! Generates a z resolved MSD
  USE histogram,  ONLY : rotcorr, vect, stride, tcorr, nattype1, &
                         nattype2, is_bound, rindex, nhist1, nhist2
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  DOUBLE PRECISION           :: dotprod(nhist1,nhist2)
  INTEGER                    :: frame1, frame2, i_at, j_at, i
  INTEGER                    :: c(nhist1,nhist2), ind(2)
  INTEGER                    :: ihist1,ihist2

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    c=0
    dotprod = 0.D0
    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,stride

!$omp parallel do private(ind) reduction(+:dotprod,c)
      DO i_at = 1,nattype1

        ind=rindex(:,i_at,frame2)

        if (((ind(1)>0).and.(ind(1) <= nhist1)).and.&
          ((ind(2)>0).and.(ind(2)<=nhist2))) then

          DO j_at = 1,nattype2

            !IF (ALL(is_bound(frame2:frame2+frame1,i_at,j_at)==1)) THEN
            IF ((is_bound(frame2,i_at,j_at)==1).and.(is_bound(frame2+frame1,i_at,j_at)==1)) THEN

              dotprod(ind(1),ind(2)) = dotprod(ind(1),ind(2)) + &
                             sum(vect(:,frame2,j_at)*vect(:,frame2+frame1,j_at))
              c(ind(1),ind(2)) = c(ind(1),ind(2)) + 1

            END IF

          END DO

        end if

      ENDDO
!$omp end parallel do

    ENDDO

    DO ihist1=1,nhist1
      DO ihist2=1,nhist2
        if (c(ihist1,ihist2).gt.0) then
          rotcorr(frame1+1,ihist1,ihist2) = dotprod(ihist1,ihist2) &
                                       / dble(c(ihist1,ihist2))
        end if
      END DO
    END DO

  ENDDO

END SUBROUTINE ROT_CORR

SUBROUTINE ASSIGN_VECTOR ( frame )
  !Select on OW atoms to compute MSD
  USE histogram,  ONLY : vect, typerot1, typerot2, is_bound, rindex, &
                         central_atom, nhist1,nhist2
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat, jat, iw, jw
  REAL*8                     :: d(4)
  REAL*8                     :: molframe(3,3)
  REAL*8, PARAMETER          :: rcut=1.4d0

  CALL Molecular_Frame(molframe)

  iw = 1
  DO iat=1,natoms

    IF ( atype(iat)==typerot1 ) THEN 

      jw=1
      CALL Vector_Distancei(iat,central_atom(2),d)
      CALL Assign_Hist_Index(d(1:3),nhist1,nhist2,molframe,rindex(:,iw,frame))

      DO jat=1,natoms

        IF ( atype(jat)==typerot2 ) THEN

          CALL DIST_UNIT_VECTOR(iat,jat,d(1:3),d(4))

          IF (d(4)<rcut) THEN

            vect(:,frame,jw) = d(1:3)
            is_bound(frame,iw,jw)=1

          END IF
          jw=jw+1

        END IF

      END DO
      iw=iw+1

    END IF

  END DO

END SUBROUTINE ASSIGN_VECTOR

SUBROUTINE Molecular_Frame(molframe)
  USE histogram, ONLY : central_atom
  IMPLICIT NONE
  REAL*8                     :: molframe(3,3), d(4)
  REAL*8                     :: r(3,3)
  LOGICAL, PARAMETER         :: is_linear=.false.

  CALL Vector_Distancei(central_atom(1),central_atom(2),d)
  r(:,1)=d(1:3)/d(4)
  CALL Vector_Distancei(central_atom(3),central_atom(2),d)
  r(:,2)=d(1:3)/d(4)
  IF (.not. is_linear ) THEN
    r(:,1) = r(:,1) + r(:,2) 
    r(:,1) = r(:,1) / norm2(r(:,1))
  END IF
  CALL CROSSPROD( r(:,1), r(:,2), r(:,3))
  CALL CROSSPROD( r(:,3), r(:,1), r(:,2))
  r(:,2)=r(:,2)/norm2(r(:,2))
  r(:,3)=r(:,3)/norm2(r(:,3))
  molframe=r

END SUBROUTINE Molecular_Frame

SUBROUTINE Assign_Hist_Index(d,nhist1,nhist2,molframe,hist_index)
  IMPLICIT NONE
  REAL*8,INTENT(IN) :: d(3), molframe(3,3)
  INTEGER, INTENT(IN) :: nhist1,nhist2
  REAL*8, PARAMETER :: dmax=6.0d0
  INTEGER :: hist_index(2), i
  REAL*8 :: r(3)

  DO i=1,2
    r(i) = sum(d*molframe(:,i+1))
  END DO
  r(3) = sqrt(r(1)**2+r(2)**2)
  hist_index(1) = int(float(nhist1)*(dmax+sum(d*molframe(:,1)))/(2*dmax)) + 1
  hist_index(2) = int(float(nhist2)*(dmax+r(3))/(2*dmax)) + 1

END SUBROUTINE assign_hist_index

SUBROUTINE COMPUTE_ROTATION_FREQUENCY
  USE histogram, only : nhist1,nhist2,rotcorr,tcorr,rotfreq
  USE parameters, only : dt
  INTEGER :: ih1,ih2
  REAL*8 :: rotation_frequency

  dif_coef=0.d0

  DO ih1=1,nhist1
    DO ih2=1,nhist2
      if (rotcorr(tcorr,ih1,ih2).ne.0) then
        rotfreq(ih1,ih2) = rotation_frequency(rotcorr(:,ih1,ih2),tcorr,dt)
      end if
    END DO
  END DO

END SUBROUTINE COMPUTE_ROTATION_FREQUENCY

REAL*8 FUNCTION rotation_frequency(rotcorr,tcorr,dt)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: tcorr
  INTEGER :: t, tlinear
  REAL*8, INTENT(IN) :: rotcorr(tcorr), dt
  REAL*8 :: mean_rot, nominator, denominator, b
  REAL*8 :: mean_time, logrotcorr(tcorr)

  logrotcorr=0.d0

  do t=1,tcorr
    if (rotcorr(t)<0.2) then
      tlinear=t
      EXIT
    end if
  end do

  do t=1,tlinear
    logrotcorr(t)=log(rotcorr(t))
  end do

  mean_rot=0.d0
  mean_time=0.d0
  nominator=0.d0
  denominator=0.d0

  do t=1,tlinear
    nominator=nominator+logrotcorr(t)*dble(t-1)*dt
    denominator=denominator+dble(t-1)**2*dt**2
    mean_time=mean_time+dble(t-1)*dt
    mean_rot=mean_rot+logrotcorr(t)
  end do

  mean_time=mean_time/dble(tlinear)
  mean_rot=mean_rot/dble(tlinear)

  b = nominator - dble(tlinear)*mean_rot*mean_time
  b = b / (denominator -dble(tlinear)*mean_time**2)

  rotation_frequency = -b

END FUNCTION rotation_frequency 

SUBROUTINE PRINT_RESULTS 
  USE histogram,  ONLY : nhist1,nhist2, rotfreq
  IMPLICIT NONE
  INTEGER                    :: ih1,ih2
  DOUBLE PRECISION           :: dh1,dh2 

  OPEN(unit = 2,file = "Rotational_Frequency.dat")

  dh1 = 12.0d0/dble(nhist1)
  dh2 = 12.0d0/dble(nhist2)

  do ih2=1,int(nhist2/2)
    rotfreq(:,ih2) = rotfreq(:,nhist2-ih2+1)
  end do

  DO ih1 = 1,nhist1
    DO ih2 = 1,nhist2
      WRITE(2,fmt = '(3(E12.5,3X))') dble(2*ih1-1)/2.d0*dh1-6.d0,&
                                     dble(2*ih2-1)/2.d0*dh2-6.d0,&
                                     rotfreq(ih1,ih2) 
    END DO
  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
