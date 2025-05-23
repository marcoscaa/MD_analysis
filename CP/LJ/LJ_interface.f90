MODULE NoseHoover
  ! Variables to perform Nose-Hoover dynamics
  IMPLICIT NONE
  REAL*8, ALLOCATABLE :: x_NH(:) ! Position of NH particle
  REAL*8, ALLOCATABLE :: v_NH(:) ! Velocity of NH particle
  REAL*8, ALLOCATABLE :: Q(:) ! Mass of NH particle
  INTEGER, PARAMETER  :: n_ys=3,n_c=1 ! Constants from MTTK paper
  INTEGER,PARAMETER   :: n_chains=3 ! Number of NH particles in a chain
  REAL*8, PARAMETER   :: T=1. ! Target temperature
  REAL*8, PARAMETER   :: tau=11. ! Frequency of NH oscillator
  REAL*8,PARAMETER,DIMENSION(3) :: w = (/ 1.351207192 , -1.702414384 , 1.351207192 /) !Martyna
  INTEGER             :: ndeg ! Number of degrees of freedom
  LOGICAL  :: tnose
  
END MODULE NoseHoover

PROGRAM MD
  !Using LJ units: sigma = 1, epsilon = 1, T = 1 (119.8 K)
  !Time: 2.18E-12s, sigma = 3.405E-10m
  USE NoseHoover, ONLY : tnose
  IMPLICIT NONE
  REAL*8,PARAMETER                    :: dt=0.005 !Both in LJ units
  INTEGER,PARAMETER                   :: nsteps=10000
  REAL*8,ALLOCATABLE                  :: pos(:,:),vel(:,:),forces(:,:)
  REAL*8,ALLOCATABLE                  :: applyforce(:,:) 
  REAL*8                              :: time,box(3),sumv2,energy,ecut,dist2
  INTEGER                             :: step,nat
  !
  ! Initialize variables for LJ MD
  !
  CALL Init(pos,vel,forces,sumv2,n_dens,nat,box,ecut,nsteps,time) 
  !
  ! Main loop for MD
  !
  do step = 1,nsteps
    !
    ! Thermostat
    !
    IF ( tnose ) CALL NH_CHAIN(sumv2,vel,nat,dt)
    !
    ! Velocity Verlet part 1 
    !
    CALL Velocity_Verlet(pos,vel,forces,nat,dt,sumv2,1)
    !
    ! Compute forces
    !
    CALL Force(forces,pos,nat,box,energy,ecut,applyforce)
    !
    ! Velocity Verlet part 2
    !
    CALL Velocity_Verlet(pos,vel,forces,nat,dt,sumv2,2)
    !
    ! Thermostat
    !
    IF ( tnose ) CALL NH_CHAIN(sumv2,vel,nat,dt)
    !
    ! Printing section
    !
    if (mod(step,10) == 0) then
      CALL Print_Results(sumv2,energy,step,nat,dt)
      CALL PRINT_POS(pos,vel,nat,step)
    end if 
    !
    ! Increment simulation time
    !
    time = time + dt
    !
  end do ! nsteps
  !
  CALL Close_Run()
  !

END PROGRAM MD

SUBROUTINE Init(pos,vel,forces,applyforce,sumv2,nat,box,ecut,nsteps,time)
  ! Generate initial pos and velocities for the system
  ! Set up initial configuration for thermostat, if needed
  USE NoseHoover
  IMPLICIT NONE
  INTEGER,INTENT(IN)                  :: nsteps
  INTEGER,INTENT(INOUT)               :: nat
  REAL*8,ALLOCATABLE                  :: pos(:,:),vel(:,:),forces(:,:)
  REAL*8,ALLOCATABLE                  :: applyforce(:,:) 
  REAL*8                              :: boxmuller, time, sumv2
  REAL*8,DIMENSION(3)                 :: sumvel
  REAL*8,INTENT(inout)                :: box(3),ecut
  INTEGER*8                           :: i,k,j
  
  tnose=.true.
  time = 0.
  forces=0.
  !
  CALL READ_INIGEO(pos,nat,box,applyforce) 
  ! 
  ALLOCATE(vel(nat,3))
  ALLOCATE(forces(nat,3))
  ndeg = sum(applyforce)
  !
  ! Gaussian distribution of initial velocities
  do i = 1,nat
    do j = 1,3
      vel(i,j) = boxmuller(DBLE(0.9),DBLE(0.)) 
      vel(i,j) = vel(i,j) * applyforce(i,j)
    end do !j
  end do !i
  !
  ! Remove the center of mass motion
  CALL Remove_COM_vel(vel,nat)
  vel = vel * applyforce
  !
  sumv2 = sum(vel*vel)
  !
  ! Define the cutoff for LJ interactions
  ecut = 4. * (1./(min(box)/2.)**6)*( (1./(min(box)/2.)**6) - 1. )
  !
  ! Starting Nose-Hoover variables
  !
  IF (tnose) THEN
    !
    ALLOCATE( x_NH(n_chains) ); x_NH=0.
    ALLOCATE( v_NH(n_chains) ); v_NH=boxmuller(DBLE(T),0.)
    ALLOCATE( Q(n_chains) ); Q=float(ndeg/3)*T/(tau*tau) 
    !
  ENDIF
  ! 
  ! Open proper output files
  OPEN(unit=1,file='energy.dat')
  OPEN(unit=5,file='pos.dat')
  OPEN(unit=6,file='vel.dat')

END SUBROUTINE Init

SUBROUTINE Close_Run()
  ! Finish everything
  USE NoseHoover
  IMPLICIT NONE
  
  CLOSE(1) ! OPENED in Init()
  CLOSE(5)
  CLOSE(6)
  
  IF (tnose) DEALLOCATE( Q, x_NH, v_NH ) 

END SUBROUTINE Close_Run

SUBROUTINE READ_INIGEO(pos,n,box,applyforce)
  !This section generates ordered initial positions, in a simple cubic
  !Cube root of n should be integer
  IMPLICIT NONE
  REAL*8, ALLOCATABLE                 :: pos(:,:)
  REAL*8, ALLOCATABLE                 :: applyforce(:,:)
  REAL*8                              :: r,box(3)
  INTEGER*8                           :: i,j,k,l
  INTEGER                             :: n,n_row

  OPEN(10, file = 'inipos.in')
  !
  READ(10,*) n
  !
  ALLOCATE(pos(n,3))
  ALLOCATE(applyforce(n,3))
  !
  DO i = 1,n
    !
    READ(10,*) pos(i,:), applyforce(i,:)
    !
  END DO
  !
  READ(10,*) box(:)
  !
  CLOSE(10)
  !
END SUBROUTINE INIGEO

REAL*8 FUNCTION boxmuller(var,x_aver)
  !Generates a random number within a gaussian distribution
  IMPLICIT NONE
  REAL*8,INTENT(IN)                   :: var,x_aver
  REAL*8                              :: rnd1,rnd2
  REAL*8,PARAMETER                    :: Pi=3.1415
  
  CALL RANDOM_NUMBER(rnd1)
  CALL RANDOM_NUMBER(rnd2)

  boxmuller = var*(x_aver + sqrt(-2.*log(rnd1))*cos(2.*Pi*rnd2))

END FUNCTION boxmuller

SUBROUTINE Force(forces,pos,nat,box,energy,ecut,applyforce)
  !Calculates Lennard-Jones forces for a single MD step
  IMPLICIT NONE
  REAL*8,DIMENSION(nat,3)             :: forces,pos,applyforce
  REAL*8,INTENT(IN)                   :: box(3)
  INTEGER,INTENT(IN)                  :: nat
  INTEGER*8                           :: i,j,k
  REAL*8                              :: dx,R2,R6,ecut,Rcut2,dist2,ff
  REAL*8,INTENT(inout)                :: energy

  Rcut2=(min(box)/2.)**2
  energy = 0.
  forces = 0.

  do i = 1,nat-1
    do j = i+1,nat

      R2 = dist2(pos(i,:),pos(j,:),box)

        if (R2 < Rcut2) then
          R6 = (1./R2)**3

          do k = 1,3

            dx = pos(i,k) - pos(j,k)
            dx = dx - box(k) * nint(dx/box(k)) !PBC
            ff = 48.*dx*R6*(R6-0.5)/R2 
            forces(i,k) = forces(i,k) + ff 
            forces(j,k) = forces(j,k) - ff 

          end do

          energy = energy + (4.*R6*(R6-1.) - ecut) !/ float(nat)

        end if ! r < rcut

    end do !j
  end do !i

  force = force * applyforce

END SUBROUTINE Force

SUBROUTINE Remove_COM_Vel(vel,nat)
  ! Remove velocity of COM
  !
  USE NoseHoover, ONLY : ndeg
  IMPLICIT NONE
  ! Parse-in variables
  INTEGER,INTENT(IN)       :: nat,ndeg
  REAL*8                   :: vel(nat,3)
  ! Internal variables
  INTEGER                  :: k
  REAL*8                   :: sumvel(3)
  !
  do k=1,3
    sumvel(k) = sum(vel(:,k)) / float(ndeg/3)
    vel(:,k) = vel(:,k) - sumvel(k)
  end do !k
  !
END SUBROUTINE Remove_COM_Vel 

REAL*8 FUNCTION dist2(c1,c2,box)
  ! Computes the square distance between two particles
  ! PBC considered
  IMPLICIT NONE
  REAL*8,DIMENSION(3)                 :: c1,c2
  INTEGER*8                           :: i
  REAL*8                              :: d,d2,box(3)

  d2 = 0.
  do i=1,3
    d = c1(i) - c2(i)
    d = d - box(i) * nint(d/box(i))
    d2 = d2 + d**2
  end do

  dist2 = d2

END FUNCTION dist2

SUBROUTINE Velocity_Verlet(pos,vel,forces,nat,dt,sumv2,flag)
  !Uses the Verlet Algorithm to propagate the equations of motion
  IMPLICIT NONE
  REAL*8,DIMENSION(nat,3)                :: pos,vel,forces
  REAL*8,INTENT(IN)                     :: dt
  REAL*8,INTENT(INOUT)                  :: sumv2
  REAL*8                                :: pos_temp
  INTEGER,INTENT(IN)                    :: nat, flag
  INTEGER*8                             :: i,j

  IF (flag .eq. 2) sumv2 = 0.

  do i = 1,nat
    !
    do j = 1,3
      !
      IF ( flag .eq. 1) THEN
        ! Propagate with old forces
        vel(i,j) = vel(i,j) + dt/2. * forces(i,j)
        pos(i,j) = pos(i,j) + vel(i,j)*dt
      ELSE
        ! Propagate with new forces
        vel(i,j) = vel(i,j) + dt/2. * forces(i,j) 
        !
        CALL Remove_COM_Vel(vel,nat)
        !
        sumv2 = sumv2 + vel(i,j)**2 
      ENDIF
      !
    end do ! j
    !
  end do ! i

END SUBROUTINE Velocity_Verlet 

SUBROUTINE NH_CHAIN(sumv2,vel,nat,dt)
  !Following the implementation of Martyna
  USE NoseHoover
  IMPLICIT NONE
  ! Parsed-in variables
  REAL*8                                  :: sumv2, vel(nat,3)
  REAL*8, INTENT(IN)                      :: dt
  INTEGER, INTENT(IN)                     :: nat
  ! Internal variables
  REAL*8,PARAMETER                        :: kb=1.
  REAL*8                                  :: dts,SCAL,G,AKIN
  INTEGER                                 :: i,j,k

  AKIN = sumv2/2.
  SCAL=1.
  ndeg = ndeg -3

  DO i=1,n_c
    DO j=1,n_ys
      dts = w(j)*dt/float(n_ys)

      IF (n_chains >= 2) THEN
        
        !Start from the end
        G = (Q(n_chains-1)*(v_NH(n_chains-1)**2) - kb*T)/Q(n_chains)
        v_NH(n_chains) = v_NH(n_chains) + G*(dts/4.)

        !until the beginning of the chain
        DO k=n_chains-1,1,-1
          v_NH(k) = v_NH(k) * exp(-dts*v_NH(k+1)/8.)
          if (k == 1 ) then
            G = (2.*AKIN - ndeg*kb*T) / Q(1)
          else 
            G = (Q(k-1)*((v_NH(k-1)**2)) - kb*T) / Q(k)
          end if
          v_NH(k) = v_NH(k) + (dts*G/4.)
          v_NH(k) = v_NH(k) * exp(-dts*v_NH(k+1)/8.)
        END DO

        !Scale factor for ekin and velocities
        SCAL = SCAL * exp(-dts*v_NH(1)/2.)
        AKIN = AKIN * exp(-dts*v_NH(1))

        !Update NH positions
        DO k=1,n_chains
          x_NH(k) = x_NH(k) + (dts*v_NH(k)/2.)
        END DO

        !Do procedure above in reverse
        DO k=1,n_chains-1
          v_NH(k) = v_NH(k) * exp(-dts*v_NH(k+1)/8.)
          if (k == 1 ) then
            G = (2*AKIN - ndeg*kb*T) / Q(1)
          else 
            G = (Q(k-1)*(v_NH(k-1)**2) - kb*T) / Q(k)
          end if
          v_NH(k) = v_NH(k) + (dts*G/4.)
          v_NH(k) = v_NH(k) * exp(-dts*v_NH(k+1)/8.)
        END DO

        G = (Q(n_chains-1)*(v_NH(n_chains-1)**2) - kb*T)/Q(n_chains)
        v_NH(n_chains) = v_NH(n_chains) + G*(dts/4.)

      ELSE
        
        G = (2.*AKIN - ndeg*kb*T) / Q(1)
        v_NH(1) = v_NH(1) + (dts*G/4.)

        SCAL = SCAL * exp(-dts*v_NH(1)/2.)
        AKIN = AKIN * exp(-dts*v_NH(1))

        x_NH(1) = x_NH(1) + (dts*v_NH(1)/2.)

        G = (2.*AKIN - ndeg*kb*T) / Q(1)
        v_NH(1) = v_NH(1) + (dts*G/4.)

      END IF
    END DO !j
  END DO !i

  !Update particle velocities
  vel = vel*SCAL
  sumv2 = AKIN*2.
  ndeg = ndeg + 3

END SUBROUTINE NH_CHAIN

SUBROUTINE Print_Results(sumv2,energy,step,nat,dt)
  USE NoseHoover
  IMPLICIT NONE
  REAL*8,INTENT(IN)                   :: sumv2,dt
  INTEGER,INTENT(IN)                  :: nat,step
  REAL*8                              :: temp,ekin,energy,econs

  ! Taking into account removal of COM velocity
  temp = sumv2 / 3. / float(nat-1)
  ekin = sumv2 / 2. 
  econs = energy + ekin 

  IF (tnose) THEN 
    econs = econs + sum(v_NH*v_NH*Q/2.) + float(ndeg-3)*T*x_NH(1) 
    IF (n_chains > 1) econs = econs + T*sum(x_NH(2:)) 
  ENDIF

  WRITE(1,fmt='(5(E14.6,3X))') float(step)*dt, temp, energy, ekin, econs

END SUBROUTINE Print_Results

SUBROUTINE PRINT_POS(pos,vel,nat,step)
  IMPLICIT NONE
  REAL*8,DIMENSION(nat,3)                       :: pos,vel
  INTEGER,INTENT(IN)                            :: step,nat
  INTEGER                                       :: j

  WRITE(5,*) step
  WRITE(6,*) step

  do j=1,nat
    WRITE(5,*) pos(j,:)
    WRITE(6,*) vel(j,:)
  end do

END SUBROUTINE PRINT_POS

REAL*8 FUNCTION norm(vector)
  IMPLICIT NONE
  REAL*8,DIMENSION(3)                           :: vector
  INTEGER                                       :: i
  REAL*8                                        :: s

  s=0.
  do i=1,3
    s = s + vector(i)**2
  end do

  norm = SQRT(s)

END FUNCTION norm
