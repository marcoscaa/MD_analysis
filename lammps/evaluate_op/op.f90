!Computes the PT CV. This is not a general code!!!
!Only works if atype is ordered with O1 ... ON H1 ... H2N

MODULE OP
  REAL*8 :: lambda
END MODULE OP

PROGRAM Coord_number
  USE parameters, ONLY : nframes, nequil
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED_CNT
    !CALL READ_XYZ_IPI
    IF (frame>nequil) THEN
      !CALL FIND_INDEX_PROTON_DEFECT(frame)
      CALL COMPUTE_OP2(frame)
    END IF
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
  !ALLOCATE(pos(3,natoms))

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

  WRITE(2,*) step, sd, Total_charge

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

  WRITE(2,fmt="(I10, 3(3X,F16.8))") step, s

END SUBROUTINE COMPUTE_OP2

SUBROUTINE FIND_INDEX_PROTON_DEFECT(step)
  USE parameters, ONLY : nwater,natoms
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: step
  REAL*8                     :: charge(nwater)
  REAL*8                     :: dist_mat(natoms-nwater,nwater)
  INTEGER                    :: iw,index_defect

  CALL COMPUTE_DIST_MATRIX(dist_mat)

  CALL compute_charge_non_smooth(dist_mat,charge)

  index_defect = 0

  DO iw = 1,nwater
    IF (abs(charge(iw)).eq.1) THEN
       index_defect=iw 
   END IF
  END DO

  WRITE(2,fmt="(I10, 1X, I4)") step, index_defect

END SUBROUTINE FIND_INDEX_PROTON_DEFECT

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
  USE parameters, only : natoms, nwater, pos
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

SUBROUTINE READ_ATOM_REDUCED_CNT 
  !Read LAMMPS atom file
  USE parameters, ONLY : pos, atype, box, natoms
  IMPLICIT NONE
  INTEGER                    :: iat, ind, iwater, nat, atyp, ipol
  REAL*8                     :: box_tmp(2), box_min(3)
  REAL*8, ALLOCATABLE        :: pos_(:,:)
 
  DO iat=1,3
    READ(1,*)
  END DO

  READ(1,*) nat
  READ(1,*)

  allocate(pos_(3,nat))

  !Read Box
  DO iat=1,3
    READ(1,*) box_tmp(1), box_tmp(2)
    box(iat) = box_tmp(2)-box_tmp(1)
    box_min(iat)=box_tmp(1)
  END DO

  READ(1,*)

  iwater=0
  DO iat = 1, nat
    READ(1,*) ind, atype(ind), pos_(1,ind), pos_(2,ind), pos_(3,ind)
    if ((atype(ind)==2).or.(atype(ind)==3)) iwater=iwater+1
  END DO

  IF (.not. allocated(pos)) THEN
      natoms=iwater
      allocate(pos(3,natoms))
  END IF

  iwater=0
  DO atyp = 2,3
    DO iat = 1, nat
      if (atype(iat)==atyp) then
          iwater=iwater+1
          pos(:,iwater) = pos_(:,iat) * box + box_min
      end if
    END DO
  END DO

  deallocate(pos_)

  !!print *, iwater, natoms
  !print *, pos
  !STOP

END SUBROUTINE READ_ATOM_REDUCED_CNT
