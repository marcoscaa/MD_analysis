
MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: msd(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: post(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: dif_coef(:,:)
  INTEGER                                  :: stride, typmsd, tcorr
  INTEGER                                  :: nhist1, nhist2, tlinear
  INTEGER                                  :: central_atom(3)
  INTEGER, ALLOCATABLE                     :: rindex(:,:,:)

END MODULE histogram

PROGRAM Hbond
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL ASSIGN_POS (frame)
  END DO

  CALL MSD_CORR 
  CALL COMPUTE_DIFFUSION_COEFFICIENT
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : post, msd, stride, typmsd, central_atom, &
                         rindex, tcorr, nhist1, nhist2, tlinear, &
                         dif_coef
  USE parameters, ONLY : natoms, nframes, nattype,&
                         pos, nequil, atype, dt
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nattype, nframes, nequil, dt
  READ(1,*) stride, typmsd, nhist1, nhist2, tcorr, tlinear
  READ(*,*) central_atom
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(post(3,nattype,nframes)); post=0 ! coordinates
  ALLOCATE(msd(tcorr,nhist1,nhist2)); msd = 0
  ALLOCATE(dif_coef(nhist1,nhist2)); dif_coef=0.d0
  ALLOCATE(atype(natoms))
  ALLOCATE(rindex(2,nattype,nframes)); rindex=0

  OPEN(unit = 1,file = file_pos)
  
END SUBROUTINE INITIALIZE

SUBROUTINE MSD_CORR
  ! Calculate the velocity msdelation function of water atoms 
  ! Generates a z resolved MSD
  USE histogram,  ONLY : msd, post, stride, tcorr, nhist1,nhist2, rindex
  USE parameters, ONLY : nframes, box, nattype
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_disp(nhist1,nhist2), dist2
  INTEGER                    :: frame1, frame2, i_at, ix
  INTEGER                    :: ind(2), ninlayer(nhist1,nhist2), ihist1, ihist2

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    sum_disp = 0.D0
    ninlayer=0
    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,stride

      DO i_at = 1,nattype

          ind=rindex(:,i_at,frame2)

          if (((ind(1)>0).and.(ind(1) <= nhist1)).and.&
            ((ind(2)>0).and.(ind(2)<=nhist2))) then

            DO ix = 1,3
              dist2 = post(ix,i_at,frame2) - post(ix, i_at, frame2+frame1)
              dist2 = dist2 - nint( dist2 / box(ix) ) * box(ix)
              dist2 = dist2*dist2

              sum_disp(ind(1),ind(2)) = sum_disp(ind(1),ind(2)) + dist2 
            END DO

            ninlayer(ind(1),ind(2))=ninlayer(ind(1),ind(2))+1

          end if

      ENDDO

    ENDDO

    DO ihist1=1,nhist1
      DO ihist2=1,nhist2
        if (ninlayer(ihist1,ihist2).gt.0) then
          msd(frame1+1,ihist1,ihist2) = sum_disp(ihist1,ihist2) &
                                       / dble(ninlayer(ihist1,ihist2))
        end if
      END DO
    END DO

  ENDDO

END SUBROUTINE MSD_CORR

SUBROUTINE COMPUTE_DIFFUSION_COEFFICIENT
  USE histogram, only : nhist1,nhist2,msd,tcorr,tlinear,dif_coef
  USE parameters, only : dt
  INTEGER :: ih1,ih2
  REAL*8 :: diffusion_coefficient

  dif_coef=0.d0

  DO ih1=1,nhist1
    DO ih2=1,nhist2
      if (msd(tcorr,ih1,ih2).ne.0) then
        dif_coef(ih1,ih2) = diffusion_coefficient(msd(:,ih1,ih2),tcorr,tlinear,dt)
      end if
    END DO
  END DO

END SUBROUTINE COMPUTE_DIFFUSION_COEFFICIENT

SUBROUTINE ASSIGN_POS ( frame )
  !Select on OW atoms to compute MSD
  USE histogram,  ONLY : post, typmsd, central_atom, rindex, nhist1, nhist2
  USE parameters, ONLY : pos, natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat, iw
  REAL*8                     :: d(4), molframe(3,3)

  CALL Molecular_Frame(molframe)

  iw = 1
 
  DO iat=1,natoms

    IF ( atype(iat)==typmsd ) THEN 
      post(:,iw,frame) = pos(:,iat)
      CALL Vector_Distancei(iat,central_atom(2),d)
      CALL Assign_Hist_Index(d(1:3),nhist1,nhist2,molframe,rindex(:,iw,frame))
      iw = iw +1
    END IF

  END DO

END SUBROUTINE ASSIGN_POS

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

REAL*8 FUNCTION diffusion_coefficient(msd,tcorr,tlinear,dt)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: tcorr, tlinear
  REAL*8, INTENT(IN) :: msd(tcorr), dt
  REAL*8 :: mean_msd, nominator, denominator, b
  REAL*8 :: mean_time
  INTEGER :: t

  mean_msd=0.d0
  mean_time=0.d0
  nominator=0.d0
  denominator=0.d0

  do t=1,tlinear
    nominator=nominator+msd(t)*dble(t-1)*dt
    denominator=denominator+dble(t-1)**2*dt**2
    mean_msd=mean_msd+msd(t)
    mean_time=mean_time+dble(t-1)*dt
  end do

  mean_msd = mean_msd/dble(tlinear)
  mean_time = mean_time/dble(tlinear)

  b = nominator - dble(tlinear)*mean_msd*mean_time
  b = b / (denominator -dble(tlinear)*mean_time**2)

  diffusion_coefficient = b/6.d0

END FUNCTION diffusion_coefficient

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : dif_coef,nhist1,nhist2
  IMPLICIT NONE
  INTEGER                    :: ih1,ih2
  REAL*8                     :: dh1,dh2

  OPEN(unit = 2,file = "diffusion_coefficient.dat")

  dh1 = 12.0d0/dble(nhist1)
  dh2 = 12.0d0/dble(nhist2)
  do ih2=1,int(nhist2/2)
    dif_coef(:,ih2) = dif_coef(:,nhist2-ih2+1)
  end do

  DO ih1 = 1,nhist1
    DO ih2 = 1,nhist2
      WRITE(2,fmt = '(3(E12.5,3X))') dble(2*ih1-1)/2.d0*dh1-6.d0,&
                                     dble(2*ih2-1)/2.d0*dh2-6.d0,&
                                     dif_coef(ih1,ih2) 
    END DO
  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
