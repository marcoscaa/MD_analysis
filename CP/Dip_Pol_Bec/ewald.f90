!All units in atomic units

PROGRAM SFG
  IMPLICIT NONE
  INTEGER, PARAMETER         :: ind_at=1, natoms=1 
  DOUBLE PRECISION           :: pos(3,natoms),dipole(3,natoms)
  DOUBLE PRECISION           :: box(3,3) 

  pos(:,1) = (/ 0,0,0 /) 
  pos(:,2) = (/ 1,0,0 /) 

  dipole(:,1) = (/ 1.d0,0.d0,0.d0 /)
  dipole(:,2) = (/ 1,0,0 /)

  box=0.d0
  box(1,1)=8.d0;box(2,2)=8.d0;box(3,3)=8.d0

  CALL INDUCED_FIELD_EWALD ( ind_at,pos,dipole,box,natoms )

  CALL LATTICE_SUM ( ind_at,pos,dipole,box,natoms )

END PROGRAM SFG

SUBROUTINE INDUCED_FIELD_EWALD ( ind_at,pos,dipole,box,natoms )
  !Calculate the re-irradiated field produced by the change in 
  !polarization after we apply the external electric field
  IMPLICIT NONE
  INTEGER, INTENT(IN)          :: ind_at, natoms
  DOUBLE PRECISION, INTENT(IN) :: pos(3,natoms), box(3,3), dipole(3,natoms)
  DOUBLE PRECISION             :: LR_field(3)
  DOUBLE PRECISION             :: SR_field(3)
  DOUBLE PRECISION             :: sum_dip(3)
  DOUBLE PRECISION             :: boxrec(3), Vcell 
  DOUBLE PRECISION             :: Kvec(3), K2 
  DOUBLE PRECISION             :: r2, fact
  DOUBLE PRECISION             :: rLvec(3), rL,rL2
  DOUBLE PRECISION             :: f1, f2
  DOUBLE PRECISION, PARAMETER  :: Pi=3.14159265359,alpha=2.5
  DOUBLE PRECISION             :: eps(3)
  COMPLEX,PARAMETER            :: i=(0,1)
  INTEGER                      :: iat, ipol, ieps
  INTEGER                      :: ik1, ik2, ik3
  INTEGER,PARAMETER            :: nk=100, nL=20
  
  Vcell = 1.d0
  DO ipol=1,3
    boxrec(ipol) =  2*Pi/box(ipol,ipol)
    Vcell = Vcell * box(ipol,ipol)
  END DO

  LR_field = 0.d0
  SR_field = 0.d0

  sum_dip = 0.d0
  DO iat=1,natoms
    sum_dip = sum_dip + dipole(:,iat)
  END DO

  eps = 0.d0

  ik2=0;ik3=0
  !Fourier part of Ewald sum
  DO ik1=-nk,nk 
    DO ik2=-nk,nk
      DO ik3=-nk,nk
        Kvec = (/ ik1,ik2,ik3 /) * boxrec
        K2 = sum(Kvec*Kvec)
        IF (K2<1.D-8) CYCLE

        f2 = 4.d0*(Pi/K2) * exp(-K2/(4*alpha)) 
        DO iat=1,natoms
          !f1 = f2 * cos(sum(Kvec*pos(:,iat))) 
          f1 = f2 * REAL(exp(i*sum(Kvec*pos(:,iat))))

          eps = eps - &
            !f1 * Kvec * sum(Kvec*dipole(:,iat)) 
            f2 * REAL(Kvec*exp(i*sum(Kvec*pos(:,iat))) & 
               * sum(Kvec*dipole(:,iat)))

        END DO

        eps = eps + &
           f2 * Kvec * sum(Kvec*dipole(:,ind_at))

      END DO
    END DO
  END DO

  print *, "###Reciprocal-space###"
  print *, eps/Vcell 
  
  !K=0 term 
  eps = eps - (4.d0*Pi/3.d0) * ( sum_dip - dipole(:,ind_at) ) 
  !eps = eps - (2.d0*Pi/3.d0) * sum_dip  
  print *, "###K=0###"
  !print *, - ((2.d0*Pi/3.d0) * sum_dip)/Vcell
  print *, - (4.d0*Pi/3.d0) * ( sum_dip - dipole(:,ind_at) )/Vcell 


  LR_field = eps/Vcell

  !Remove self-interaction
  !LR_field = LR_field + (2.d0/3.d0) * (alpha**3/Pi)**(0.5d0) * dipole(:,ind_at)
  print *, "###Self-interaction###"
  print *, (2.d0/3.d0) * (alpha**3/Pi)**(0.5d0) * dipole(:,ind_at)
  
  !Short range part of Ewald sum
  DO ik1=-nL,nL
    DO ik2=-nL,nL
      DO ik3=-nL,nL

        DO iat=1,natoms

          rLvec = pos(:,iat) + MATMUL((/ ik1, ik2, ik3 /),box)
          rL2 = sum(rLvec*rLvec)
          rL = sqrt(rL2)
          IF (rL.eq.0) CYCLE

          f1 = 2.d0*sqrt(alpha/Pi) * rL * exp( -alpha*rL2 ) &
             + erfc(sqrt(alpha)*rL)
          f2 = 4.d0*alpha*sqrt(alpha/Pi) * exp( -alpha*rL2 ) / rL2
          !f1 = 3.d0*erfc(sqrt(alpha)*rL)/(rL2*rL) &
          !   + 6.d0*sqrt(alpha/Pi) * Exp(-alpha*rL2)/rL2 &
          !   + 4.d0*alpha*sqrt(alpha/Pi) * Exp(-alpha*rL2)
          !f1 = f1 / rL2
          !     
          !f2 = erfc(sqrt(alpha)*rL)/rL &
          !   + 2.d0*sqrt(alpha/Pi)*exp(-alpha*rL2)
          !f2 = f2 / rL2

          SR_field = SR_field & 
                   + f1*3.d0*rLvec * sum(rLvec*dipole(:,iat)) & 
                   / (rL2*rL2*rL) &
                   - f1*dipole(:,iat) / (rL2*rL) &
                   + f2* rLvec * sum( rLvec*dipole(:,iat) )
          !SR_field(:) = SR_field(:) &
          !                 + f1 * rLvec * sum( rLvec*dipole(:,iat) ) & 
          !                 - f2 * dipole(:,iat) 
         
        END DO

      END DO
    END DO
  END DO

  print *, "###Long-Range###"
  print *, LR_field
  print *, "###Short-Range###"
  print *, SR_field
  print *, "###Total###"
  print *, LR_field + SR_field
    
END SUBROUTINE INDUCED_FIELD_EWALD

SUBROUTINE LATTICE_SUM ( ind_at,pos,dipole,box,natoms )
  !Calculate the re-irradiated field produced by the change in 
  !polarization after we apply the external electric field
  IMPLICIT NONE
  INTEGER, INTENT(IN)          :: ind_at, natoms
  DOUBLE PRECISION, INTENT(IN) :: pos(3,natoms), box(3,3), dipole(3,natoms)
  DOUBLE PRECISION             :: rLvec(3), rL,rL2
  DOUBLE PRECISION             :: field(3)
  DOUBLE PRECISION, PARAMETER  :: Pi=3.14159265359
  INTEGER                      :: iat, ipol, ieps
  INTEGER                      :: ik1, ik2, ik3
  INTEGER,PARAMETER            :: nL=200

  field = 0.d0
  ik3=0;ik2=0

  DO ik1=-nL,nL
    DO ik2=-nL,nL
      DO ik3=-nL,nL

        DO iat=1,natoms

          rLvec = pos(:,iat) + MATMUL((/ ik1, ik2, ik3 /),box)
          rL2 = sum(rLvec*rLvec)
          rL = sqrt(rL2)
          IF (rL .lt. 1.D-8) CYCLE

          field = field + 3.d0 * sum(dipole(:,iat)*rLvec) &
                        * rLvec / (rL2*rL2) &
                        - dipole(:,iat) / (rL2*rL)
         
        END DO

     END DO
   END DO
  END DO

  print *, "###Lattice-Sum###"
  print *, field
    
END SUBROUTINE LATTICE_SUM
