!Computes the PT CV. This is not a general code!!!
!Only works if atype is ordered with O1 ... ON H1 ... H2N

MODULE OP
  REAL*8 :: lambda
END MODULE OP

PROGRAM Coord_number
  USE parameters, ONLY : nframes, box
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_POS_CEL
    CALL COMPUTE_OP(frame)
  END DO

  CLOSE(1); CLOSE(2)

END PROGRAM Coord_number

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nwater, pos, atype, &
                         nframes, nequil
  USE OP, ONLY : lambda
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file, out_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)
  CALL getarg(3, out_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nwater, nframes, nequil, lambda
  CLOSE(1)

  !Allocate module arrays
  ALLOCATE(atype(natoms))
  ALLOCATE(pos(3,natoms))

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = out_file)

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
  REAL*8                     :: Dist
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

  WRITE(2,*) step, sd, sum(abs(charge))

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

  WRITE(2,fmt="(I10, 3(3X,F12.8))") step, s

END SUBROUTINE COMPUTE_OP2

real*8 function compute_charge(ind, dist_mat)
  use parameters, only: nwater, natoms
  use op, only : lambda
  implicit none
  integer, intent(in)        :: ind
  integer                    :: iw, ih
  real*8                     :: denominator, coord
  real*8                     :: dih, dist_mat(natoms-nwater,nwater)

  coord = 0.d0
  do ih = nwater+1,natoms
    denominator = 0.d0
    dih=dist_mat(ih-nwater,ind)
    if (dih<3.0) then 
!$omp parallel do reduction(+:denominator)
      do iw = 1,nwater
        denominator = denominator + exp(-dble(lambda*dist_mat(ih-nwater,iw)))
      end do
!$end omp parallel do
    coord = coord + exp(-lambda*dih)/denominator
    end if
  end do

  compute_charge = coord - 2.d0

end function compute_charge

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
