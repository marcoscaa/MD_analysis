PROGRAM MD
  !Using LJ units: sigma = 1, epsilon = 1, T = 1 (119.8 K)
  !Time: 2.18E-12s, sigma = 3.405E-10m
  IMPLICIT NONE
  REAL*8,PARAMETER                    :: dt=0.005,n_dens=0.8178 !Both in LJ units
  REAL*8                              :: pbox 
  INTEGER,PARAMETER                   :: nsteps=100000,np=125,n_bins=100,n_eq=1000
  REAL*8,DIMENSION(np,3)              :: old_pos,pos,vel,forces
  REAL*8,DIMENSION(nsteps,np,3)       :: vel_corr,pos_corr
  INTEGER*8,DIMENSION(n_bins)         :: gof
  REAL*8                              :: time,box,sumv2,energy,ecut,dist2,msdisp
  INTEGER                             :: step
   
  time = 0.
  gof = 0.
  msdisp = 0.
  CALL Init(pos,old_pos,vel,n_dens,np,box,ecut,dt,nsteps) 
  do step = 1,nsteps
    CALL Force(forces,pos,np,box,energy,ecut)
    CALL New_Coordinates(pos,old_pos,vel,forces,np,dt,sumv2)
    CALL Print_Results(sumv2,energy,step,np,vel,dt)
    if ((step >= n_eq) .and. (mod(step,10) == 0)) then
      CALL PRINT_POS(pos,vel,np,step)
    end if 
    time = time + dt
  end do ! nsteps
  CLOSE(1) ! OPEN in Init()
  CLOSE(5)
  CLOSE(6)

END PROGRAM MD

SUBROUTINE Init(pos,old_pos,vel,n_dens,np,box,ecut,dt,nsteps)
  !Generate initial coordinates and velocities for the system
  IMPLICIT NONE
  REAL*8,DIMENSION(np,3)              :: pos,old_pos,vel
  REAL*8,INTENT(IN)                   :: n_dens,dt
  INTEGER,INTENT(IN)                  :: np,nsteps
  REAL*8                              :: boxmuller
  REAL*8,DIMENSION(3)                 :: sumvel
  REAL*8,INTENT(inout)                :: box,ecut
  INTEGER*8                           :: i,j,k,l
  
  pos=0.
  vel=0.
  box = (DBLE(np) / n_dens) ** (1./3.)
  !Dimple cubic initial configuration
  CALL INIGEO(pos,np,box) 
  
  do i = 1,np
    do j = 1,3
      !Here I define the initial velocity
      vel(i,j) = boxmuller(DBLE(0.9),DBLE(0.)) !125:0.9,343:0.83,512:0.85
    end do !j
  end do !i

  !Remove the center of mass motion
  do k=1,3
    sumvel(k) = sum(vel(:,k)) / DBLE(np)
    vel(:,k) = vel(:,k) - sumvel(k)
  end do !k

  do i = 1,np
    do j = 1,3
      old_pos(i,j) = pos(i,j) - vel(i,j)*dt
    end do !j
  end do !i

  ecut = 4. * (1./(box/2.)**6)*( (1./(box/2.)**6) - 1. )

  OPEN(unit=1,file='energy.dat')
  OPEN(unit=5,file='pos.dat')
  OPEN(unit=6,file='vel.dat')

  WRITE(5,*) nsteps/10,np,dt*10,box
  WRITE(6,*) nsteps/10,np,dt*10,box

END SUBROUTINE Init

SUBROUTINE INIGEO(coordinates,n,box)
  !This section generates ordered initial positions, in a simple cubic
  !Cube root of n should be integer
  IMPLICIT NONE
  REAL*8, DIMENSION(n,3)              :: coordinates
  REAL*8                              :: r,box
  INTEGER*8                           :: i,j,k,l
  INTEGER                           :: n,n_row

  n_row = n**(1./3.)
  l=1
  do i=0,n_row-1
    do j=0,n_row-1
      do k=0,n_row-1
        coordinates(l,1) = box/(2.*DBLE(n_row)) + (box/DBLE(n_row))*i
        coordinates(l,3) = box/(2.*DBLE(n_row)) + (box/DBLE(n_row))*j
        coordinates(l,2) = box/(2.*DBLE(n_row)) + (box/DBLE(n_row))*k
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

SUBROUTINE Force(forces,pos,np,box,energy,ecut)
  !Calculates Lennard-Jones forces for a single MD step
  IMPLICIT NONE
  REAL*8,DIMENSION(np,3)              :: forces,pos
  REAL*8,INTENT(IN)                   :: box
  INTEGER,INTENT(IN)                :: np
  INTEGER*8                           :: i,j,k
  REAL*8                              :: dx,R2,R6,ecut,Rcut2,dist2,ff
  REAL*8,INTENT(inout)                :: energy

  Rcut2=(box/2.)**2
  energy = 0.
  forces = 0.

  do i = 1,np-1
    do j = i+1,np

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

          energy = energy + (4*R6*(R6-1.) - ecut) / DBLE(np)

        end if ! r < rcut

    end do !j
  end do !i

END SUBROUTINE Force

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

SUBROUTINE New_Coordinates(pos,old_pos,vel,forces,np,dt,sumv2)
  !Uses the Verlet Algorithm to propagate the equations of motion
  IMPLICIT NONE
  REAL*8,DIMENSION(np,3)                :: pos,old_pos,vel,forces
  REAL*8,INTENT(IN)                     :: dt
  INTEGER,INTENT(IN)                    :: np                            
  INTEGER*8                             :: i,j
  REAL*8,INTENT(INOUT)                  :: sumv2
  REAL*8                                :: pos_temp

  sumv2 = 0.

  do i = 1,np
    do j = 1,3
      pos_temp = pos(i,j)
      pos(i,j) = 2*pos(i,j) - old_pos(i,j) + forces(i,j)*dt*dt
      vel(i,j) = (pos(i,j) - old_pos(i,j)) / (2.*dt)
      old_pos(i,j) = pos_temp
      sumv2 = sumv2 + vel(i,j)**2 
    end do ! j
  end do ! i

END SUBROUTINE New_Coordinates

SUBROUTINE Print_Results(sumv2,energy,step,np,vel,dt)
  IMPLICIT NONE
  REAL*8,INTENT(IN)                   :: sumv2,dt
  INTEGER,INTENT(IN)                  :: np,step
  REAL*8                              :: temp,ekin,energy,mod_vel
  REAL*8,DIMENSION(np,3),INTENT(IN)   :: vel
  INTEGER*8                           :: i,j

  temp = sumv2 / 3. / DBLE(np)
  ekin = sumv2 / 2. 
  energy = energy 

  WRITE(1,fmt='(5(F12.8,3X))') step*dt, temp, energy, ekin, energy+ekin

END SUBROUTINE Print_Results

SUBROUTINE PRINT_POS(pos,vel,np,step)
  IMPLICIT NONE
  REAL*8,DIMENSION(np,3)                        :: pos,vel
  INTEGER,INTENT(IN)                            :: step,np
  INTEGER                                       :: j

  WRITE(5,*) step
  WRITE(6,*) step

  do j=1,np
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
    dV = (4./3.)* Pi * DBLE(i**3 - (i-1)**3) * bin_size**3
    avr_rdf = gof(i) / (dV * N_particles * n_dens * Nsteps)
    WRITE(2,*) (DBLE(i-1)+0.5) * bin_size , avr_rdf
  end do ! i

  CLOSE(2)

END SUBROUTINE AVRG_RDF

SUBROUTINE MSD(pos,np,box,nsteps,dt)
  !Calculate the Mean Squere Displacement of all particles 
  IMPLICIT NONE
  REAL*8,DIMENSION(nsteps,np,3),INTENT(in)   :: pos
  REAL*8,INTENT(IN)                          :: box,dt
  INTEGER,INTENT(IN)                       :: np,nsteps
  REAL*8                                     :: dist2,msdisp
  INTEGER*8                                  :: i,j,k

  OPEN(unit=3,file='msd.dat')

  do i=1,nsteps
    msdisp = 0.
    do j=1,nsteps-i
      do k=1,np
        msdisp = msdisp + dist2(pos(j,k,:),pos(i+j,k,:),box)
      end do ! k
    end do !i

    WRITE(3,*) i*dt, msdisp/DBLE(np*j)
  end do !j

  CLOSE(3) ! OPEN in Init()


END SUBROUTINE MSD

SUBROUTINE Kubo(vel_corr,nsteps,np,dt)
  !Calculate the velocity auto-corr function and determine D Using
  !Kubo's formula
  IMPLICIT NONE
  REAL*8,DIMENSION(nsteps,np,3)           :: vel_corr
  INTEGER,INTENT(IN)                    :: nsteps,np
  INTEGER*8                               :: i,j,k,l
  REAL*8,DIMENSION(nsteps)                :: v_c
  REAL*8                                  :: trapz,vc_temp,av_vel,dt

  OPEN(unit=4,file='vac.dat')
  v_c = 0.
  trapz = 0.

  do i=1,nsteps
    do j=1,nsteps-i
      do k=1,np 
    
        do l=1,3
          v_c(i) = v_c(i) + vel_corr(j,k,l)*vel_corr(i+j,k,l)
        end do

      end do !k
    end do !j
    v_c(i) = v_c(i) / DBLE(3*np*(j))
  end do !i
  
  do i=1,nsteps-1000
    WRITE(4,*) i*dt, v_c(i)
    trapz = trapz + dt*(v_c(i) + v_c(i+1)) / 2.
  end do !i

  PRINT *, trapz*dt

  CLOSE(4)

END SUBROUTINE Kubo

SUBROUTINE Struc_fac(pos,np,nsteps,box)
  IMPLICIT NONE
  REAL*8,DIMENSION(nsteps,np,3),INTENT(IN)         :: pos
  INTEGER,INTENT(IN)                             :: np,nsteps
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
        qv = (/ DBLE(g), DBLE(h), DBLE(i) /) *2*Pi/box
        sf=0.
        do j=1,nsteps
          sf1 = 0.
          sf2 = 0.
          !No need to include imaginary numbers in the structure factor formula
          do k=1,np
            sf1 = sf1 + cos(dot_product(qv,pos(j,k,:)))
            sf2 = sf2 + sin(dot_product(qv,pos(j,k,:)))
          end do !k
          sf = sf + sf1**2 + sf2**2
        end do !j
         total(c,1) = norm(qv)
         total(c,2) = sf / DBLE(np*nsteps)
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
  dw = (maxi - mini) / DBLE(nw)

  hist = 0.

  !Check each window for the points inside it
  do i=1,nw
    lower = mini + DBLE(i-1)*dw
    upper = mini + DBLE(i)*dw
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


