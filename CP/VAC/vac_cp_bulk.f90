!Set of subroutines used to analyse the lipid trajectory

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: vel(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:)
  DOUBLE PRECISION                         :: box(9)
  INTEGER*8, ALLOCATABLE                   :: ind_atom(:)
  DOUBLE PRECISION                         :: dt ! in ps
  INTEGER                                  :: natoms, nframes 
  INTEGER                                  :: nequil, index_equil 

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  !Remove equilibration part
  CALL REMOVE_EQUIL 

  DO frame = 1,nframes
    CALL READ_FRAME (frame)
  END DO

  CALL VAC 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, notused
  CHARACTER(100)             :: file_vel, index_file

  CALL getarg(1, file_vel)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil, notused, dt
  READ(1, *), box

  ALLOCATE( ind_atom(natoms) )

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
  END DO

  CLOSE(1)

  nframes = nframes - nequil
  ALLOCATE(vel(3,natoms,nframes)) ! velocities of water atoms
  ALLOCATE(autocorr(nframes)); autocorr = 0

  OPEN(unit = 1,file = file_vel)
  
END SUBROUTINE INITIALIZE

SUBROUTINE REMOVE_EQUIL
  !Subroutine reads gro files
  USE parameters, ONLY : nequil, natoms 
  IMPLICIT NONE
  INTEGER                    :: i,j 

  DO i = 1, nequil
    DO j = 1,natoms+1
      READ(1, *)
    END DO
  END DO

END SUBROUTINE

SUBROUTINE READ_FRAME (index_frame)
  !Subroutine reads gro files
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, index_frame

  READ(1,*)

  DO i = 1,natoms

    !Store velocities and positions of water atoms only
      READ(1, *), & 
          vel(1,i,index_frame), vel(2,i,index_frame),vel(3,i,index_frame)

  END DO

END SUBROUTINE

SUBROUTINE VAC
  ! Calculate the velocity autocorrelation function of water atoms 
  ! Generates a z resolved VAC
  USE parameters
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_vel, sum_vel_temp
  INTEGER                    :: frame1, frame2, i_at, i_coord
  INTEGER                    :: counter, counter_tmp 

    !frame1 is the duration of the interval 
    DO frame1 = 0,nframes-1

      counter = 0
      sum_vel = 0.D0

      !Change initial condition, look for combinations with same interval
      DO frame2 = 1,nframes-frame1

        counter_tmp = 0
        sum_vel_temp = 0.D0

!$omp parallel do reduction(+:counter_tmp, sum_vel_temp)
        DO i_at = 1,natoms

            counter_tmp = counter_tmp + 1

            DO i_coord = 1,3
              sum_vel_temp = sum_vel_temp + &
                           & vel(i_coord, i_at, frame2) * &
                           & vel(i_coord, i_at, frame2+frame1)
            END DO

        ENDDO
!$omp end parallel do

        counter = counter + counter_tmp
        sum_vel = sum_vel + sum_vel_temp

      ENDDO

      IF ( counter /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(frame1+1) = sum_vel / ( dble( counter )  )!dble( nframes-frame1 ) )
      END IF

    ENDDO

END SUBROUTINE VAC 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j

  OPEN (unit=10, file='vac.dat')

  DO i = 0,nframes-1

    write(10,fmt = '(E11.4,3X,E12.5)'), dble(i)*dt, autocorr(i+1) / autocorr(1) 

  ENDDO

  CLOSE(10)

END SUBROUTINE PRINT_RESULTS
