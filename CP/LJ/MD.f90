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
  LOGICAL  :: tnose
  
END MODULE NoseHoover

PROGRAM MD
  !Using LJ units: sigma = 1, epsilon = 1, T = 1 (119.8 K)
  !Time: 2.18E-12s, sigma = 3.405E-10m
  USE NoseHoover, ONLY : tnose
  IMPLICIT NONE
  REAL*8,PARAMETER                    :: dt=0.005,n_dens=0.8178 !Both in LJ units
  REAL*8                              :: pbox 
  INTEGER,PARAMETER                   :: nsteps=10000,nat=125
  REAL*8,DIMENSION(nat,3)             :: pos,vel,forces
  REAL*8                              :: time,box,sumv2,energy,ecut,dist2
  INTEGER                             :: step
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
    CALL Force(forces,pos,nat,box,energy,ecut)
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

SUBROUTINE Init(pos,vel,forces,sumv2,n_dens,nat,box,ecut,nsteps,time)
  ! Generate initial coordinates and velocities for the system
  ! Set up initial configuration for thermostat, if needed
  USE NoseHoover
  IMPLICIT NONE
  INTEGER,INTENT(IN)                  :: nat,nsteps
  REAL*8,DIMENSION(nat,3)             :: pos,vel,forces
  REAL*8,INTENT(IN)                   :: n_dens
  REAL*8                              :: boxmuller, time, sumv2
  REAL*8,DIMENSION(3)                 :: sumvel
  REAL*8,INTENT(inout)                :: box,ecut
  INTEGER*8                           :: i,k,j
  
  tnose=.true.
  time = 0.
  forces=0.
  box = (float(nat) / n_dens) ** (1./3.)
  !
  ! Simple cubic initial configuration
  CALL INIGEO(pos,nat,box) 
  ! 
  ! Gaussian distribution of initial velocities
  do i = 1,nat
    do j = 1,3
      vel(i,j) = boxmuller(DBLE(0.9),DBLE(0.)) 
    end do !j
  end do !i
  !
  ! Remove the center of mass motion
  CALL Remove_COM_vel(vel,nat)
  !
  sumv2 = sum(vel*vel)
  !
  ! Define the cutoff for LJ interactions
  ecut = 4. * (1./(box/2.)**6)*( (1./(box/2.)**6) - 1. )
  !
  ! Starting Nose-Hoover variables
  !
  IF (tnose) THEN
    !
    ALLOCATE( x_NH(n_chains) ); x_NH=0.
    ALLOCATE( v_NH(n_chains) ); v_NH=boxmuller(DBLE(T),0.)
    ALLOCATE( Q(n_chains) ); Q=float(nat)*T/(tau*tau) 
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

SUBROUTINE INIGEO(coordinates,n,box)
  !This section generates ordered initial positions, in a simple cubic
  !Cube root of n should be integer
  IMPLICIT NONE
  REAL*8, DIMENSION(n,3)              :: coordinates
  REAL*8                              :: r,box
  INTEGER*8                           :: i,j,k,l
  INTEGER                             :: n,n_row

  n_row = n**(1./3.)
  l=1
  do i=0,n_row-1
    do j=0,n_row-1
      do k=0,n_row-1
        coordinates(l,1) = box/(2.*float(n_row)) + (box/float(n_row))*i
        coordinates(l,3) = box/(2.*float(n_row)) + (box/float(n_row))*j
        coordinates(l,2) = box/(2.*float(n_row)) + (box/float(n_row))*k
        l=l+1
      end do
    end do
  end do 

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

SUBROUTINE Force(forces,pos,nat,box,energy,ecut)
  !Calculates Lennard-Jones forces for a single MD step
  IMPLICIT NONE
  REAL*8,DIMENSION(nat,3)              :: forces,pos
  REAL*8,INTENT(IN)                   :: box
  INTEGER,INTENT(IN)                  :: nat
  INTEGER*8                           :: i,j,k
  REAL*8                              :: dx,R2,R6,ecut,Rcut2,dist2,ff
  REAL*8,INTENT(inout)                :: energy

  Rcut2=(box/2.)**2
  energy = 0.
  forces = 0.

  do i = 1,nat-1
    do j = i+1,nat

      R2 = dist2(pos(i,:),pos(j,:),box)

        if (R2 < Rcut2) then
          R6 = (1./R2)**3

          do k = 1,3

            dx = pos(i,k) - pos(j,k)
            dx = dx - box * nint(dx/box) !PBC
            ff = 48.*dx*R6*(R6-0.5)/R2 
            forces(i,k) = forces(i,k) + ff 
            forces(j,k) = forces(j,k) - ff 

          end do

          energy = energy + (4.*R6*(R6-1.) - ecut) !/ float(nat)

        end if ! r < rcut

    end do !j
  end do !i

END SUBROUTINE Force

SUBROUTINE Remove_COM_Vel(vel,nat)
  ! Remove velocity of COM
  !
  IMPLICIT NONE
  ! Parse-in variables
  INTEGER,INTENT(IN)       :: nat
  REAL*8                   :: vel(nat,3)
  ! Internal variables
  INTEGER                  :: k
  REAL*8                   :: sumvel(3)
  !
  do k=1,3
    sumvel(k) = sum(vel(:,k)) / float(nat)
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
  REAL*8                              :: d,d2,box

  d2 = 0.
  do i=1,3
    d = c1(i) - c2(i)
    d = d - box * nint(d/box)
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
  REAL*8                                  :: dts,SCAL,G,AKIN,ndeg
  INTEGER                                 :: i,j,k

  AKIN = sumv2/2.
  ! Number of degrees of freedom. Removing COM displacement
  ndeg = float(nat-1)*3.

  SCAL=1.

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
    econs = econs + sum(v_NH*v_NH*Q/2.) + 3.*float(nat-1)*T*x_NH(1) 
    IF (n_chains > 1) econs = econs + T*sum(x_NH(2:)) 
  ENDIF

  WRITE(1,fmt='(5(E14.6,3X))') float(step)*dt, temp, energy, ekin, econs

END SUBROUTINE Print_Results

SUBROUTINE PRINT_POS(pos,vel,nat,step)
  IMPLICIT NONE
  REAL*8,DIMENSION(nat,3)                        :: pos,vel
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ANALYSIS SECTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RDF(n,xyz,gof,n_bins,box)
  !Calculate the RDF of a single step and average over all hard spheres
  IMPLICIT NONE
  REAL*8                               :: r,dist2,box,bin_size
  INTEGER*8                            :: i,j,w,n_bins
  REAL*8, DIMENSION(n, 3)              :: xyz
  INTEGER*8, DIMENSION(n_bins)         :: gof 
  INTEGER, INTENT(IN)                :: n

  bin_size = box / (2.*n_bins) 
  do i = 1,n-1            ! considering all possible distances
    do j = i+1,n          ! between particles with no double couting

      r = SQRT(dist2(xyz(i,:),xyz(j,:),box)) 
      w = nint(r/bin_size) !index for the rdf matrix
      gof(w) = gof(w) + 2 ! counting i,j and j,i

    end do ! j
  end do ! i

END SUBROUTINE RDF

SUBROUTINE AVRG_RDF(gof,Nsteps,N_particles,n_bins,box,n_dens)
  !Average the rdf over all steps
  IMPLICIT NONE
  INTEGER*8, DIMENSION(n_bins)             :: gof
  REAL*8                                   :: avr_rdf,dV,bin_size
  REAL*8, INTENT(IN)                       :: box,n_dens
  INTEGER, INTENT(IN)                    :: Nsteps,N_particles,n_bins
  INTEGER*8                                :: i
  REAL*8, PARAMETER                        :: Pi=3.1415926

  OPEN(unit=2,file='rdf.dat')

  bin_size = box / (2.*n_bins)

  do i = 1, n_bins 
    dV = (4./3.)* Pi * float(i**3 - (i-1)**3) * bin_size**3
    avr_rdf = gof(i) / (dV * N_particles * n_dens * Nsteps)
    WRITE(2,*) (float(i-1)+0.5) * bin_size , avr_rdf
  end do ! i

  CLOSE(2)

END SUBROUTINE AVRG_RDF

SUBROUTINE MSD(pos,nat,box,nsteps,dt)
  !Calculate the Mean Squere Displacement of all particles 
  IMPLICIT NONE
  REAL*8,DIMENSION(nsteps,nat,3),INTENT(in)   :: pos
  REAL*8,INTENT(IN)                          :: box,dt
  INTEGER,INTENT(IN)                       :: nat,nsteps
  REAL*8                                     :: dist2,msdisp
  INTEGER*8                                  :: i,j,k

  OPEN(unit=3,file='msd.dat')

  do i=1,nsteps
    msdisp = 0.
    do j=1,nsteps-i
      do k=1,nat
        msdisp = msdisp + dist2(pos(j,k,:),pos(i+j,k,:),box)
      end do ! k
    end do !i

    WRITE(3,*) i*dt, msdisp/float(nat*j)
  end do !j

  CLOSE(3) ! OPEN in Init()


END SUBROUTINE MSD

SUBROUTINE Kubo(vel_corr,nsteps,nat,dt)
  !Calculate the velocity auto-corr function and determine D Using
  !Kubo's formula
  IMPLICIT NONE
  REAL*8,DIMENSION(nsteps,nat,3)           :: vel_corr
  INTEGER,INTENT(IN)                    :: nsteps,nat
  INTEGER*8                               :: i,j,k,l
  REAL*8,DIMENSION(nsteps)                :: v_c
  REAL*8                                  :: trapz,vc_temp,av_vel,dt

  OPEN(unit=4,file='vac.dat')
  v_c = 0.
  trapz = 0.

  do i=1,nsteps
    do j=1,nsteps-i
      do k=1,nat 
    
        do l=1,3
          v_c(i) = v_c(i) + vel_corr(j,k,l)*vel_corr(i+j,k,l)
        end do

      end do !k
    end do !j
    v_c(i) = v_c(i) / float(3*nat*(j))
  end do !i
  
  do i=1,nsteps-1000
    WRITE(4,*) i*dt, v_c(i)
    trapz = trapz + dt*(v_c(i) + v_c(i+1)) / 2.
  end do !i

  PRINT *, trapz*dt

  CLOSE(4)

END SUBROUTINE Kubo

SUBROUTINE Struc_fac(pos,nat,nsteps,box)
  IMPLICIT NONE
  REAL*8,DIMENSION(nsteps,nat,3),INTENT(IN)         :: pos
  INTEGER,INTENT(IN)                             :: nat,nsteps
  INTEGER*8                                        :: i,j,k,g,h,c
  REAL*8                                           :: sf,sf1,sf2,dist2,box,norm
  REAL*8,DIMENSION(3)                              :: qv
  REAL*8,DIMENSION(1771,2)                         :: total
  REAL,PARAMETER                                   :: Pi=3.1415926

  c=0
  !The results will be put altogether in a matrix
  !c is the index counter for this matrix
  do g=0,20
    do h=g,20
      do i=h+1,21
        c = c+1
        !Only allowed q vectors are computed 
        qv = (/ float(g), float(h), float(i) /) *2*Pi/box
        sf=0.
        do j=1,nsteps
          sf1 = 0.
          sf2 = 0.
          !No need to include imaginary numbers in the structure factor formula
          do k=1,nat
            sf1 = sf1 + cos(dot_product(qv,pos(j,k,:)))
            sf2 = sf2 + sin(dot_product(qv,pos(j,k,:)))
          end do !k
          sf = sf + sf1**2 + sf2**2
        end do !j
         total(c,1) = norm(qv)
         total(c,2) = sf / float(nat*nsteps)
         !PRINT *, total(c,1),total(c,2)
      end do !i
    end do !h
  end do !g

  CALL Make_Histogram(total,1771,200)

END SUBROUTINE Struc_fac

SUBROUTINE Make_Histogram(vector,len_vec,nw)
  !This histogram will compute the average of the vector values inside
  !the determined window
  IMPLICIT NONE
  REAL*8,DIMENSION(len_vec,2)                     :: vector
  REAL*8,DIMENSION(nw,2)                          :: hist
  INTEGER,INTENT(IN)                              :: len_vec
  INTEGER                                         :: i,j,nw
  REAL*8                                          :: mini,maxi,dw,lower,upper

  OPEN(unit=5,file='str_fac.dat')
 
  mini = minval(vector(:,1))
  maxi = maxval(vector(:,1))
  dw = (maxi - mini) / float(nw)

  hist = 0.

  !Check each window for the points inside it
  do i=1,nw
    lower = mini + float(i-1)*dw
    upper = mini + float(i)*dw
    !For each point in the matrix
    do j=1,len_vec

        if (vector(j,1) >= lower .AND. vector(j,1) < upper) then
          hist(i,1) = hist(i,1) + vector(j,2)
          hist(i,2) = hist(i,2) + 1.
        end if

    end do

    if (hist(i,2) > 0.) then
      WRITE(5,*) (lower+upper)/2., hist(i,1) / hist(i,2)
    end if
  end do

  CLOSE(5)
  
END SUBROUTINE Make_Histogram 


