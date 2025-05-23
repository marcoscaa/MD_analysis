MODULE histogram 
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: velt(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: posm(:,:)
  CHARACTER(5)                             :: atype_vac
END MODULE histogram

PROGRAM VAC 
  USE parameters, ONLY : nframes, corr_time
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_EXTXYZ 
    CALL COMPUTE_VEL(frame)
    CALL ASSIGN_VEL (frame) 
  END DO

  CALL VEL_AUTOCORR 
  CALL PRINT_RESULTS

  CLOSE(2)

END PROGRAM VAC

SUBROUTINE INITIALIZE
  USE histogram, ONLY  : velt, autocorr, atype_vac, posm
  USE parameters, ONLY : pos, vel, natoms, atype, &
                         nframes, dt, corr_time, nequil
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_pos, index_file

  CALL getarg(1, file_pos)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, corr_time, dt
  READ(*,*) atype_vac
  CLOSE(1)

  ALLOCATE(vel(3,natoms)) 
  ALLOCATE(pos(3,natoms)) 
  ALLOCATE(atype(natoms)) 
  ALLOCATE(posm(3,natoms)) 
  ALLOCATE(velt(3,natoms,nframes)) ! velocities of all atoms
  ALLOCATE(autocorr(3,corr_time)); autocorr = 0

  OPEN(unit = 1,file = file_pos)
  
END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_VEL(frame)
  USE parameters, ONLY : vel, pos, natoms, dt
  USE histogram, ONLY : posm
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER :: iat
  
  IF (frame==1) THEN
    vel=0.0
  ELSE
    DO iat=1,natoms
      CALL DISTANCE_VECTOR(pos(:,iat),posm(:,iat),vel(:,iat))
      vel(:,iat)=vel(:,iat)/dt
    END DO
  END IF
 
  posm=pos

END SUBROUTINE COMPUTE_VEL

SUBROUTINE ASSIGN_VEL(frame)
  USE histogram, ONLY : velt
  USE parameters, ONLY : vel
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  
  velt(:,:,frame) = vel

END SUBROUTINE ASSIGN_VEL

SUBROUTINE VEL_AUTOCORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE histogram,  ONLY : velt, autocorr, atype_vac
  USE parameters, ONLY : nframes, natoms, corr_time, atype 
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_vel(3), dotprod(3)
  INTEGER                    :: frame1, frame2, i_at, i_coord, c_ini, c_tmp

  !frame1 is the duration of the interval 
  DO frame1 = 0,corr_time-1

    c_ini=0
    !e initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,1

      sum_vel = 0.D0
      c_tmp = 0

!$omp parallel do private(dotprod) reduction(+:sum_vel,c_tmp)
      DO i_at = 1,natoms

        if ((trim(atype(i_at)).eq.trim(atype_vac)).or.(trim(atype_vac).eq."All")) then
          dotprod = 0.d0

          DO i_coord = 1,3
            dotprod(i_coord) = dotprod(i_coord) + &
                    & velt(i_coord,i_at,frame2) * &
                    & velt(i_coord,i_at,frame2+frame1)
          END DO

          sum_vel = sum_vel + dotprod
          c_tmp = c_tmp + 1
        END IF

      ENDDO
!$omp end parallel do

      autocorr(:,frame1+1) = autocorr(:,frame1+1) + sum_vel
      c_ini = c_ini + c_tmp

    ENDDO

      autocorr(:,frame1+1) = autocorr(:,frame1+1) / dble( c_ini ) 

  ENDDO

END SUBROUTINE VEL_AUTOCORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,   ONLY : autocorr
  USE parameters, ONLY : box, dt, corr_time
  IMPLICIT NONE
  INTEGER                    :: i, j

  OPEN(unit = 3,file = 'vac.dat')

  DO i = 1,corr_time

    write(3,fmt = '(E11.4,3(3X,E12.5))') dble(i-1)*dt, autocorr(:,i) 

  ENDDO

  CLOSE(3)

END SUBROUTINE PRINT_RESULTS

