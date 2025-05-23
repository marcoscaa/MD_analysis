! ---------------------------------------------------------------------------------
! COMPUTES GLOBAL AND LOCAL BOND ORIENTATIONAL ORDER OF SUPERCOOLED LIQUID
! ---------------------------------------------------------------------------------

MODULE param
  
  IMPLICIT NONE
  SAVE
  
  INTEGER, PARAMETER :: npart = 2500, nsite = 3
  INTEGER, PARAMETER :: nframe = 100  
  INTEGER, PARAMETER :: max_nbhrs = 50
  INTEGER, PARAMETER :: nxtout = 25    

  INTEGER, PARAMETER :: l = 6          ! order of spherical harmonics
  INTEGER, PARAMETER :: n_ql_neigh = 16 ! no. of particles involved in order computation 
  INTEGER, PARAMETER :: flag_sorting = 1    ! flag for flag_sorting 'n_q6_neigh' neighbors from neighbor list

  INTEGER, PARAMETER :: flag_local = 1 ! flag for local OP computation
  INTEGER, PARAMETER :: Nbin  = 101
  DOUBLE PRECISION, PARAMETER :: ql_min = 0.0d0, ql_max = 1.0d0

  DOUBLE PRECISION, PARAMETER :: dt = 0.002*nxtout
  DOUBLE PRECISION,PARAMETER :: bcut = 5.0d0   ! cutoff for neighbor list
  DOUBLE PRECISION,PARAMETER :: pi = 3.141592653589d0

  DOUBLE PRECISION,PARAMETER :: av = 6.022141290d23
  DOUBLE PRECISION,PARAMETER :: mw = 18.015280d0 
 
END MODULE param

! ==================================================================================

PROGRAM main
  
  USE param
  
  IMPLICIT NONE
  
  INTEGER :: i, imol, isite, num(npart,nsite)
  INTEGER :: cframe, iframe, corr_frame
  INTEGER :: countf
  INTEGER :: bin 
 
  DOUBLE PRECISION,DIMENSION (npart) :: rx, ry, rz
  DOUBLE PRECISION :: box(3) 
  DOUBLE PRECISION :: rw(npart, nsite, 3)
  DOUBLE PRECISION :: Ql_g, Ql_l, Q_l(nframe), Q_g(nframe), Q_loc(npart)
  DOUBLE PRECISION :: sum_Ql, avg_Ql, sum_Qg, avg_Qg
  DOUBLE PRECISION :: dql, sum_qlocal(Nbin), avg_qlocal(Nbin), qlocal(Nbin)
  DOUBLE PRECISION :: n_dens, m_dens

  CHARACTER (LEN = 6) :: junk1
  CHARACTER (LEN = 3) :: junk2
  
  COMMON / coords / rx, ry, rz, box
  
  
  OPEN (UNIT = 1, FILE = "traj.xyz", ACTION = "READ")
  OPEN (UNIT = 11, FILE = "Q4_global_density.dat", ACTION = "WRITE")
  OPEN (UNIT = 12, FILE = "Q4_local.dat", ACTION = "WRITE")
  OPEN (UNIT = 13, FILE = "local_dist.dat", ACTION = "WRITE")

  Q_l = 0.0d0
  Q_g = 0.0d0
  sum_Qg  = 0.0d0
  sum_Ql  = 0.0d0

! bin size alomg order

  dql = (ql_max - ql_min) / (Nbin - 1)
  sum_qlocal = 0.0d0

  DO iframe = 1, nframe
     
     rw  = 0.0d0    
     box = 0.0d0

     READ(1,*)
     READ(1,*) box(:)
     DO imol = 1, npart
        DO isite = 1, nsite
           !READ(1,*) junk1, junk2, num(imol,isite), rw(imol,isite,:)
           READ(1,*) junk2, rw(imol,isite,:)
        END DO
     END DO
   
     n_dens = dble(npart) / (box(1)*box(2)*box(3))
     m_dens = (n_dens*mw) / (av*1.0d-24)   

     !    Extract coordinates of oxygen 
 
     DO imol = 1, npart
        rx(imol) = rw(imol,1,1)
        ry(imol) = rw(imol,1,2)
        rz(imol) = rw(imol,1,3)
     END DO
    
     !    Golbal Ql of the system     
 
     Ql_g = 0.0
     CALL Q_global(Ql_g)
     Q_g(iframe) = Ql_g
     sum_Qg = sum_Qg + Q_g(iframe)

     WRITE(11,*) (iframe-1)*dt, Q_g(iframe),  m_dens

     !   Local Ql of the system 

     IF (flag_local .EQ. 1) THEN
        Ql_l = 0.0d0
        Q_loc = 0.0d0
        CALL Q_local(Ql_l, Q_loc)    
        Q_l(iframe) = Ql_l
        sum_Ql    = sum_Ql + Q_l(iframe)   
        DO i = 1, npart
           bin = ((Q_loc(i) - ql_min) / dql) + 1
           sum_qlocal(bin) = sum_qlocal(bin) + 1
        END DO 
        WRITE(12,*) iframe*dt, Q_l(iframe)
     END IF

  END DO
  CLOSE(1)
  CLOSE(11)
  CLOSE(12)
  
  avg_Qg = sum_Qg / nframe 
  PRINT*, "global order = ", avg_Qg

  IF (flag_local .EQ. 1) THEN

     avg_Ql = sum_Ql / nframe
     PRINT*, "local order = ", avg_Ql

     !  Computation of local OP distribution

     DO bin = 1, Nbin
        qlocal(bin) = ql_min + (bin-1)*dql
        avg_qlocal(bin) = sum_qlocal(bin) / (nframe*npart*dql)
        WRITE (13,*) qlocal(bin), avg_qlocal(bin)
     END DO
     CLOSE(13)

  END IF 

END PROGRAM main

! ==================================================================================

SUBROUTINE Q_global(Ql_g)
  
  USE param
  IMPLICIT NONE
  
  INTEGER :: i, m 
  
  DOUBLE PRECISION :: avrg_rey(-l:l), avrg_imy(-l:l) 
  DOUBLE PRECISION :: rey(npart,-l:l), imy(npart,-l:l) 
  DOUBLE PRECISION :: Ql_g
  
  COMMON / spherical_harmonics / rey, imy
  
  CALL nn_list

  avrg_rey = 0.0d0
  avrg_imy = 0.0d0
 
  DO i = 1, npart
     CALL bond_order(i)
     DO m = -l, l
        avrg_rey(m) = avrg_rey(m) + rey(i,m)
        avrg_imy(m) = avrg_imy(m) + imy(i,m)
     END DO
  END DO
  
  Ql_g = 0.0d0
  
  DO m = -l, l
     avrg_rey(m) = avrg_rey(m) / npart
     avrg_imy(m) = avrg_imy(m) / npart
     Ql_g = Ql_g + (avrg_rey(m)**2 + avrg_imy(m)**2)
  END DO
  
  Ql_g = DSQRT(4*pi*Ql_g / (2*l+1))
  
  RETURN
END SUBROUTINE Q_global

! ===================================================================================

SUBROUTINE Q_local(Ql_l, Q_loc)

  USE param
  IMPLICIT NONE

  INTEGER :: i, m

  DOUBLE PRECISION :: avrg_ql, sum_loc
  DOUBLE PRECISION :: rey(npart,-l:l), imy(npart,-l:l)
  DOUBLE PRECISION :: Ql_l, Q_loc(npart)

  COMMON / spherical_harmonics / rey, imy

  CALL nn_list

  avrg_ql = 0.0d0

  DO i = 1, npart
     CALL bond_order(i)
     sum_loc = 0.0d0
     DO m = -l, l
        sum_loc = sum_loc + (rey(i,m)**2 + imy(i,m)**2)
     END DO
     Q_loc(i) =  DSQRT(4*pi*sum_loc / (2*l+1))
     avrg_ql = avrg_ql + Q_loc(i)
  END DO

  Ql_l = avrg_ql / npart

  RETURN
END SUBROUTINE Q_local

! ----------------------------------------------------------------------------------

SUBROUTINE bond_order(i)
  
  USE param
  IMPLICIT NONE
  
  INTEGER :: i, j, k, m
  INTEGER :: num_nbhrs(npart), nbhr_list(max_nbhrs, npart)
  
  DOUBLE PRECISION :: phi, cosmphi, sinmphi
  DOUBLE PRECISION :: rey(npart,-l:l), imy(npart,-l:l)
  DOUBLE PRECISION :: re, im
  DOUBLE PRECISION :: rxij, ryij, rzij
  DOUBLE PRECISION :: fact, legendre
  DOUBLE PRECISION :: rij_list(4,max_nbhrs,npart)

  EXTERNAL fact, legendre

  COMMON / nbr_list / num_nbhrs, nbhr_list 
  COMMON / nbr_dist / rij_list
  COMMON / spherical_harmonics / rey, imy
  
  DO m = -l, l
     rey(i,m) = 0
     imy(i,m) = 0
  END DO
  
  DO k = 1, num_nbhrs(i)
     j = nbhr_list(k, i)

     rxij = rij_list(1,k,i) / rij_list(4,k,i)
     ryij = rij_list(2,k,i) / rij_list(4,k,i)
     rzij = rij_list(3,k,i) / rij_list(4,k,i)
    
     phi  = ATAN2(ryij, rxij)

     DO m = 0, l
        cosmphi = DCOS(m*phi)
        sinmphi = DSIN(m*phi)
        
        re = DSQRT((2*l+1)*fact(l-m) / (4*pi*fact(l+m)))*legendre(m, rzij)*cosmphi
        im = DSQRT((2*l+1)*fact(l-m) / (4*pi*fact(l+m)))*legendre(m, rzij)*sinmphi
        rey(i,m) = rey(i,m) + re
        imy(i,m) = imy(i,m) + im
        
        IF (m .NE. 0) THEN
           rey(i,-m) = rey(i,-m) + (-1)**m * re
           imy(i,-m) = imy(i,-m) - (-1)**m * im
        END IF
        
     END DO 
  END DO                   
  
  DO m = -l, l
     rey(i,m) = rey(i,m) / num_nbhrs(i)
     imy(i,m) = imy(i,m) / num_nbhrs(i)
  END DO
  
  RETURN
END SUBROUTINE bond_order

! ----------------------------------------------------------------------
! Function to calculate legendre

FUNCTION legendre(m, x)

  USE param 
  IMPLICIT NONE
  
  INTEGER :: k, m, mtmp, m1
  
  DOUBLE PRECISION :: x, legendre
  DOUBLE PRECISION :: sum
  DOUBLE PRECISION :: fact
  
  EXTERNAL fact
  
  mtmp = ABS(m)
  m1   = INT((l-mtmp)/2)
  
  sum  = 0.0d0
  
  DO k = 0, m1
     sum = sum + (-1)**k*fact(2*(l-k))*x**(l-2*k-mtmp) / (2**l*fact(k)*fact(l-k)*fact(l-2*k-mtmp))
  END DO
  
  legendre = ((-1)**mtmp)*((1-x**2)**(1.0d0*mtmp/2))*sum
  
  RETURN
END FUNCTION legendre

! --------------------------------------------------
! Function to calculate factorial

FUNCTION fact(k)
  IMPLICIT NONE
  
  INTEGER :: i, k
  DOUBLE PRECISION :: fact
  
  IF (k .LE. 1) THEN
     fact = 1.0d0
     RETURN
  END IF
  
  fact = 1
  
  DO i = 1, k
     fact = fact * i
  END DO
  
  RETURN

END FUNCTION fact

! ====================================================================

SUBROUTINE nn_list 
  
  USE param
  
  IMPLICIT NONE 
  
  INTEGER :: imol, jmol, j, k 
  INTEGER :: num_nbhrs(npart), nbhr_list(max_nbhrs, npart)
  
  DOUBLE PRECISION, DIMENSION (npart) :: rx,  ry, rz
  DOUBLE PRECISION, DIMENSION (3) :: box
  DOUBLE PRECISION, DIMENSION (4,max_nbhrs,npart) :: rij_list
  DOUBLE PRECISION, DIMENSION (4) :: rij_temp

  DOUBLE PRECISION :: rxij, ryij, rzij, rij
  DOUBLE PRECISION :: nneigh_temp
  
  COMMON / nbr_list / num_nbhrs, nbhr_list
  COMMON / nbr_dist / rij_list
  COMMON / coords / rx, ry, rz, box
 
  num_nbhrs = 0
  nbhr_list = 0
  rij = 0.0d0
 
  ! neighbor list computation
  
  DO imol = 1, npart
     DO jmol = 1, npart
        
        IF (jmol .le. imol) CYCLE
        
        ! calculate the com separation vector
        rxij = rx(jmol) - rx(imol)
        ryij = ry(jmol) - ry(imol)
        rzij = rz(jmol) - rz(imol)
        
        ! apply the periodic boundary conditions
        rxij = rxij - dnint(rxij/box(1))*box(1)
        ryij = ryij - dnint(ryij/box(2))*box(2)
        rzij = rzij - dnint(rzij/box(3))*box(3)
        
        ! calculate the square of the com scalar separation distance
        rij = dsqrt(rxij**2 + ryij**2 + rzij**2)

        ! IF they are within the cutoff distance
        IF (rij .lt. bcut) THEN
           num_nbhrs(imol) = num_nbhrs(imol) + 1
           num_nbhrs(jmol) = num_nbhrs(jmol) + 1
           nbhr_list(num_nbhrs(imol), imol) = jmol
           nbhr_list(num_nbhrs(jmol), jmol) = imol

           ! store distance component
           
           rij_list(1, num_nbhrs(imol), imol) = rxij
           rij_list(2, num_nbhrs(imol), imol) = ryij
           rij_list(3, num_nbhrs(imol), imol) = rzij
           rij_list(4, num_nbhrs(imol), imol) = rij
          
           rij_list(:,num_nbhrs(jmol),jmol) = rij_list(:,num_nbhrs(imol),imol)
           
        END IF

     END DO

     ! check the minimum number of neighbors are present
     IF (num_nbhrs(imol) .lt. n_ql_neigh) THEN
        WRITE(*,'(a,i5,a,i3,a)') 'fatal error: ', imol, ' has less than ', n_ql_neigh,' nearest neighbors'
        STOP
     END IF
     
     ! check that the maximum number of neighbors is not exceeded
     IF (num_nbhrs(imol) .gt. max_nbhrs) THEN
        WRITE(*,'(a,i5,a)') 'fatal error: ', imol, ' has more than 200 nearest neighbors'
        STOP
     END IF
     
  END DO

  ! sort the neighbors using insertion sort IF using only the closest
  ! n_ql_neigh neighbors
  
  IF (flag_sorting .eq. 1) THEN

     IF (n_ql_neigh .NE. 0) THEN
        DO imol = 1, npart
           DO k = 2, num_nbhrs(imol)
              nneigh_temp = nbhr_list(k,imol)
              rij_temp(:) = rij_list(:,k,imol)
              DO j = k-1, 1, -1
                 IF (rij_list(4,j,imol) .le. rij_temp(4)) GOTO 10
                 nbhr_list(j + 1, imol) = nbhr_list(j,imol)
                 rij_list(:,j + 1,imol) = rij_list(:, j, imol)
              END DO
              jmol = 0
10            CONTINUE 
              nbhr_list(j + 1, imol) = nneigh_temp
              rij_list(:,j + 1, imol) = rij_temp(:)
           END DO
           num_nbhrs(imol) = n_ql_neigh
        END DO
     END IF
     
  END IF

  RETURN

END SUBROUTINE nn_list
