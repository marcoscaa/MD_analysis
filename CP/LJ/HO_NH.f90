PROGRAM HO_NH
  !Single harmonic oscilator with NH thermostat
  IMPLICIT NONE
  REAL*8,PARAMETER                      :: T=1.,dt=0.01,m=1. !m is the particle's mass
  INTEGER,PARAMETER                     :: nsteps=2000000,n_chains=2
  REAL*8                                :: x,v,AKIN
  REAL*8,DIMENSION(n_chains)            :: Q,x_NH,v_NH !DIMENSION = n_chains
  INTEGER                               :: step
  
  CALL Init(x,x_NH,v,v_NH,T,Q,n_chains,AKIN)
  do step=1,nsteps

    CALL NH_CHAIN(AKIN,x,v,x_NH,v_NH,Q,T,n_chains,dt)
    CALL Verlet(x,v,m,dt,AKIN)
    CALL NH_CHAIN(AKIN,x,v,x_NH,v_NH,Q,T,n_chains,dt)

    if (mod(step,1000) == 0) then
      CALL Print_Results(x,x_NH,v,v_NH,T,Q,step,n_chains)
    end if

  end do
  CALL Finish()

END PROGRAM HO_NH

SUBROUTINE Verlet(x,v,m,dt,AKIN)
  !Uses the Verlet Algorithm to propagate the equations of motion
  !with Nose-Hoover thermostat 
  IMPLICIT NONE
  REAL*8                                :: x,v,m,AKIN
  REAL*8,INTENT(IN)                     :: dt

  !F = -x
  v = v + dt/(2.*m) * (-x)
  x = x + dt*v
  !CALL FORCE
  v = v + dt/(2.*m) * (-x)
  AKIN = m*v*v/2.

END SUBROUTINE Verlet

SUBROUTINE NH_CHAIN(AKIN,x,v,x_NH,v_NH,Q,T,n_chains,dt)
  !Following the implementation of Martyna
  IMPLICIT NONE
  REAL*8,PARAMETER                        :: kb=1.
  INTEGER,PARAMETER                       :: n_ys=3,n_c=1,Np=1 !Martyna's paper
  REAL*8                                  :: x,v,dt,dts,SCAL,G,T,AKIN
  REAL*8,DIMENSION(n_chains)              :: x_NH,v_NH,Q
  REAL*8,DIMENSION(n_ys)                  :: w = (/ 1.351207192 , -1.702414384 , 1.351207192 /) !Martyna
  INTEGER                                 :: n_chains,i,j,k

  SCAL=1.

  DO i=1,n_c
    DO j=1,n_ys
      dts = w(j)*dt/DBLE(n_ys)

      IF (n_chains >= 2) THEN
        
        !Start from the end
        G = (Q(n_chains-1)*(v_NH(n_chains-1)**2) - kb*T)/Q(n_chains)
        v_NH(n_chains) = v_NH(n_chains) + G*(dts/4.)

        !until the beginning of the chain
        DO k=n_chains-1,1,-1
          v_NH(k) = v_NH(k) * exp(-dts*v_NH(k+1)/8.)
          if (k == 1 ) then
            G = (2.*AKIN - DBLE(Np)*kb*T) / Q(1)
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
            G = (2*AKIN - DBLE(Np)*kb*T) / Q(1)
          else 
            G = (Q(k-1)*(v_NH(k-1)**2) - kb*T) / Q(k)
          end if
          v_NH(k) = v_NH(k) + (dts*G/4.)
          v_NH(k) = v_NH(k) * exp(-dts*v_NH(k+1)/8.)
        END DO

        G = (Q(n_chains-1)*(v_NH(n_chains-1)**2) - kb*T)/Q(n_chains)
        v_NH(n_chains) = v_NH(n_chains) + G*(dts/4.)

      ELSE
        
        G = (2.*AKIN - DBLE(Np)*kb*T) / Q(1)
        v_NH(1) = v_NH(1) + (dts*G/4.)

        SCAL = SCAL * exp(-dts*v_NH(1)/2.)
        AKIN = AKIN * exp(-dts*v_NH(1))

        x_NH(1) = x_NH(1) + (dts*v_NH(1)/2.)

        G = (2.*AKIN - DBLE(Np)*kb*T) / Q(1)
        v_NH(1) = v_NH(1) + (dts*G/4.)

      END IF
    END DO !j
  END DO !i

  !Update particle velocities
  v = v*SCAL

END SUBROUTINE NH_CHAIN

SUBROUTINE Init(x,x_NH,v,v_NH,T,Q,n_chains,AKIN)
  !Set the initial parameters to run the simulation
  implicit NONE
  REAL*8                                :: x,v,T,AKIN
  REAL*8,DIMENSION(n_chains)            :: x_NH,v_NH,Q
  INTEGER                               :: n_chains

  OPEN(unit=1,file='pos_vel.dat')
  OPEN(unit=2,file='pos_vel_NH.dat')
  OPEN(unit=3,file='ener.dat')

  !particle in the origin with maximum kinetic energy
  x = 1.
  v = 1.!2.*(T**(0.5))-0.1
  AKIN = v*v/2.
  x_NH = 0. 
  v_NH = 0.
  Q = 1.

END SUBROUTINE Init

SUBROUTINE Print_Results(x,x_NH,v,v_NH,T,Q,step,n_chains)
  IMPLICIT NONE
  REAL*8                                :: x,v,T
  REAL*8,DIMENSION(n_chains)            :: x_NH,v_NH,Q
  INTEGER                               :: step,n_chains
  REAL*8                                :: Epot, Ekin, E_cons

  Epot = 0.5*x*x 
  Ekin = 0.5*v*v
  E_cons = Epot + Ekin + sum(Q*v_NH*v_NH/2.) + T*x_NH(1) + T*sum(x_NH(2:))

  WRITE(1,fmt='(I10, 2(3X,e15.8))') step, x, v 
  WRITE(2,fmt='(I10, 2(3X,e15.8))') step, x_NH(1), v_NH(1)
  WRITE(3,fmt='(I10, 4(3X,e15.8))') step, Epot, Ekin, Epot+Ekin, E_cons

END SUBROUTINE Print_Results

SUBROUTINE Finish()
  !Dummy subroutine to close the program
  IMPLICIT NONE

  CLOSE(1);CLOSE(2);CLOSE(3)

END SUBROUTINE Finish
