!Computes the PT CV. This is not a general code!!!
!Only works if atype is ordered with O1 ... ON H1 ... H2N

MODULE OP
  REAL*8 :: lambda
END MODULE OP

PROGRAM Coord_number
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  
  DO frame = 1,nframes
    CALL READ_RAW_POS_BOX
    CALL READ_RAW_FORCE
    CALL COMPUTE_OP_DERIVATIVE(frame)
  END DO

  CLOSE(1); CLOSE(2)

END PROGRAM Coord_number

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nwater, pos, atype, &
                         nframes, nequil, force
  USE OP, ONLY : lambda
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nwater, nframes, nequil, lambda
  CLOSE(1)

  !Read atype from type.raw
  ALLOCATE(atype(natoms))
  OPEN(unit = 1,file = 'type.raw')
  READ(1, *) (atype(i), i=1,natoms)
  CLOSE(1)

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(force(3,natoms))

  OPEN(unit = 1,file = 'coord.raw')
  OPEN(unit = 2,file = 'box.raw')
  OPEN(unit = 5,file = 'force.raw')
  OPEN(unit = 3,file = 'S_w.dat')

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_OP(step)
  !Compute center of mass along the z direction
  !Here we use the particle with index 1 as a reference
  USE parameters, ONLY : pos, nwater, natoms
  USE OP, ONLY: lambda
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: step
  REAL*8                     :: charge(nwater), compute_charge
  REAL*8                     :: dist_mat(natoms-nwater,nwater), sd
  REAL*8                     :: Dist, Total_charge
  REAL*8, PARAMETER          :: alpha=0.0001 !Smoothness of total charge
  INTEGER                    :: iw,jw

  CALL COMPUTE_DIST_MATRIX(dist_mat)

  IF (lambda.lt.0) then
    CALL compute_charge_non_smooth(dist_mat,charge)
  else
    DO iw = 1, nwater
      charge(iw) = compute_charge(iw,dist_mat)
    END DO
  end if

  sd = 0.d0

  DO iw = 1,nwater
!$omp parallel do reduction(+:sd)
    DO jw = iw+1,nwater
      sd = sd - charge(iw)*charge(jw)*Dist(iw,jw)
    END DO
!$end omp parallel do
  END DO

  IF (lambda.lt.0) then

    Total_charge = sum(abs(charge))

  ELSE
    Total_charge=0.d0

!$omp parallel do reduction(+:Total_charge)
    DO iw = 1,nwater
      Total_charge = Total_charge + sqrt(charge(iw)**2+alpha)
    END DO
!$end omp parallel do

    Total_charge = Total_charge - nwater*sqrt(alpha)

  ENDIF

  WRITE(3,*) step, sd, Total_charge

END SUBROUTINE COMPUTE_OP

SUBROUTINE COMPUTE_OP2(step)
  !Compute center of mass along the z direction
  !Here we use the particle with index 1 as a reference
  USE parameters, ONLY : pos, nwater, natoms
  USE OP, ONLY: lambda
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: step
  REAL*8                     :: charge(nwater), compute_charge,Dist
  REAL*8                     :: dist_mat(natoms-nwater,nwater), s(3)
  INTEGER                    :: iw,jw

  CALL COMPUTE_DIST_MATRIX(dist_mat)

  IF (lambda.lt.0) then
    CALL compute_charge_non_smooth(dist_mat,charge)
  else
    DO iw = 1, nwater
      charge(iw) = compute_charge(iw,dist_mat)
    END DO
  end if

  s = 0.d0

!$omp parallel do reduction(+:s)
  DO iw = 1,nwater
    s = s + abs(charge(iw))*pos(:,iw)
  END DO
!$end omp parallel do

  WRITE(3,fmt="(I10, 3(3X,F12.8))") step, s

END SUBROUTINE COMPUTE_OP2

REAL*8 FUNCTION compute_charge(ind, dist_mat)
  USE parameters, only: nwater, natoms
  USE op, ONLY : lambda
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  INTEGER                    :: iw, iH
  REAL*8                     :: denominator, coord
  REAL*8                     :: dih, dist_mat(natoms-nwater,nwater)

  coord = 0.d0
  DO iH = nwater+1,natoms
    denominator = 0.d0
    dih=dist_mat(iH-nwater,ind)
    IF (dih<3.0) THEN 
!$omp parallel do reduction(+:denominator)
      DO iw = 1,nwater
        denominator = denominator + exp(-lambda*dist_mat(iH-nwater,iw))
      END DO
!$end omp parallel do
    coord = coord + exp(-lambda*dih)/denominator
    END IF
  END DO

  compute_charge = coord - 2.d0

END FUNCTION compute_charge

subroutine compute_charge_non_smooth(dist_mat, charge)
  use parameters, only: nwater, natoms
  implicit none
  integer                    :: iw, ih, min_ind
  real*8                     :: charge(nwater)
  real*8                     :: mindist, dist_mat(natoms-nwater,nwater)

  charge = 0.d0
  do ih = nwater+1,natoms
    mindist=100.d0
    do iw = 1,nwater
      if (dist_mat(ih-nwater,iw)<mindist) then
        mindist=dist_mat(ih-nwater,iw)
        min_ind=iw
      endif
    end do
    charge(min_ind) = charge(min_ind) + 1.d0
  end do

  charge = charge - 2.d0

end subroutine compute_charge_non_smooth

SUBROUTINE COMPUTE_DIST_MATRIX(mat)
  USE parameters, only : natoms, nwater
  IMPLICIT NONE
  REAL*8                     :: mat(natoms-nwater,nwater)
  REAL*8                     :: Dist
  INTEGER                    :: iw, ih

  DO ih=nwater+1,natoms
    DO iw=1,nwater
      mat(ih-nwater,iw) = Dist(ih,iw)
    END DO
  END DO

END SUBROUTINE COMPUTE_DIST_MATRIX

SUBROUTINE COMPUTE_DIST_MATRIX_VECTOR_FULL(mat)
  USE parameters, only : natoms, pos
  IMPLICIT NONE
  REAL*8                     :: mat(4,natoms,natoms)
  REAL*8                     :: dij(4)
  INTEGER                    :: iat,jat

  DO iat=1,natoms
    DO jat=iat,natoms
      CALL DISTANCE_VECTOR(pos(:,jat),pos(:,iat),dij) !pos(iat)-pos(jat)
      mat(:,iat,jat)=dij
      mat(4,jat,iat)=dij(4)
      mat(1:3,jat,iat)=-dij(1:3)
    END DO
  END DO

END SUBROUTINE COMPUTE_DIST_MATRIX_VECTOR_FULL

SUBROUTINE COMPUTE_OP_DERIVATIVE(step)
  !Compute center of mass along the z direction
  !Here we use the particle with index 1 as a reference
  USE parameters, ONLY : pos, nwater, natoms, force
  USE OP, ONLY: lambda
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: step
  REAL*8                     :: charge(nwater), compute_charge
  REAL*8                     :: distmat(4,natoms,natoms), sd
  REAL*8                     :: dH_dS, dS_dR(3,natoms), dij(4)
  REAL*8                     :: dQ_dR(3,nwater,natoms), den 
  INTEGER                    :: iw,jw,iat, iH

  CALL COMPUTE_DIST_MATRIX_VECTOR_FULL(distmat)

  DO iw = 1, nwater
    charge(iw) = compute_charge(iw,distmat(4,nwater+1:natoms,1:nwater))
  END DO

  dQ_dR=0.d0

  DO iw=1,nwater

    !Do derivative for Oxygen first
    DO iat=1,nwater
      DO iH=nwater+1,natoms 
        den = 0.d0
        DO jw=1,nwater
          den = den + exp(-lambda*distmat(4,iH,jw))
        END DO
        if (iat==iw) then
          dQ_dR(:,iw,iat) = dQ_dR(:,iw,iat) - lambda*distmat(1:3,iat,iH)/distmat(4,iat,iH) * & 
                                              exp(-lambda*distmat(4,iat,iH))/den
        end if
        dQ_dR(:,iw,iat) = dQ_dR(:,iw,iat) + lambda*distmat(1:3,iat,iH)/distmat(4,iat,iH) * &
                                            exp(-lambda*distmat(4,iat,iH)) * &
                                            exp(-lambda*distmat(4,iw,iH))/den/den
      END DO
    END DO

    !Do derivative for hydrogen 
    DO iat=nwater+1,natoms
      den = 0.d0
      DO jw=1,nwater
        den = den + exp(-lambda*distmat(4,iat,jw))
      END DO
      dQ_dR(:,iw,iat) = dQ_dR(:,iw,iat) +lambda*distmat(1:3,iw,iat)/distmat(4,iw,iat) * &
                                                exp(-lambda*distmat(4,iw,iat))/den
      DO jw=1,nwater
        dQ_dR(:,iw,iat) = dQ_dR(:,iw,iat) -lambda*distmat(1:3,jw,iat)/distmat(4,jw,iat) * &
                                                exp(-lambda*distmat(4,jw,iat))* & 
                                                exp(-lambda*distmat(4,iw,iat))/den/den
      END DO
    END DO
  END DO

  dS_dR=0.d0

  DO iat=1,natoms
    DO iw=1,nwater
      DO jw=iw+1,nwater
        dS_dR(:,iat) = dS_dR(:,iat) - &
                      ( dQ_dR(:,iw,iat)*charge(jw) + &
                      charge(iw)*dQ_dR(:,jw,iat) ) * distmat(4,iw,jw) 
      END DO
      if ((iat.le.nwater).and.(iat.ne.iw)) then
        dS_dR(:,iat) = dS_dR(:,iat) - charge(iw)*charge(iat)*distmat(1:3,iw,iat)/distmat(4,iw,iat)
      end if
    END DO
  END DO
  
  dH_dS = 0.d0

  DO iat = 1,natoms
    dH_dS = dH_dS - sum(force(:,iat)/dS_dR(:,iat))
  END DO

  sd = 0.d0

  DO iw = 1,nwater
    DO jw = iw+1,nwater
      sd = sd - charge(iw)*charge(jw)*distmat(4,iw,jw)
    END DO
  END DO

  WRITE(3,*) step, sd, dH_dS

END SUBROUTINE COMPUTE_OP_DERIVATIVE
