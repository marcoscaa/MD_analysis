!Computes the PT CV. This is not a general code!!!
!Only works if H and O atoms belong to water molecules (or ions) only. 

PROGRAM Water_Ion 
  USE, intrinsic :: iso_fortran_env, Only : iostat_end
  IMPLICIT NONE
  INTEGER :: frame, iostat

  CALL INITIALIZE
  
  frame=1
  DO 
    CALL READ_EXTXYZ_IO (iostat)
    IF (iostat == iostat_end ) THEN
      EXIT
    END IF
    CALL REALLOCATE_POS
    CALL COMPUTE_OP(frame)
    frame=frame+1
  END DO

  CLOSE(1); CLOSE(2)

END PROGRAM Water_Ion

SUBROUTINE INITIALIZE
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = "pos_water_ion.dat")
  WRITE(2,fmt="(A24)") "# Step Pos_x Pos_y Pos_z"

END SUBROUTINE INITIALIZE

SUBROUTINE REALLOCATE_POS
  USE parameters, only : pos, nwater, nH, natoms, atype
  IMPLICIT NONE
  INTEGER                    :: iw, iat
  REAL*8                     :: pos_(3,natoms)

  iw=0
  DO iat=1,natoms
    IF (TRIM(atype(iat))=='O') THEN
        iw = iw+1
        pos_(:,iw) = pos(:,iat)
    END IF
  END DO

  nwater=iw

  DO iat=1,natoms
    IF (TRIM(atype(iat))=='H') THEN
        iw = iw+1
        pos_(:,iw) = pos(:,iat)
    END IF
  END DO

  nH=iw-nwater

  DEALLOCATE(pos)
  ALLOCATE(pos(3,nwater+nH))

  DO iw=1,nwater+nH
    pos(:,iw) = pos_(:,iw)
  END DO

END SUBROUTINE REALLOCATE_POS

SUBROUTINE COMPUTE_OP(step)
  !Compute center of mass along the z direction
  !Here we use the particle with index 1 as a reference
  USE parameters, ONLY : pos, nwater, nH
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: step
  REAL*8                     :: charge(nwater), compute_charge,Dist
  REAL*8                     :: dist_mat(nH,nwater), s(3)
  INTEGER                    :: iw,jw

  CALL COMPUTE_DIST_MATRIX(dist_mat)
  CALL compute_charge_non_smooth(dist_mat,charge)

  s = 0.d0

!$omp parallel do reduction(+:s)
  DO iw = 1,nwater
    s = s + abs(charge(iw))*pos(:,iw)
  END DO
!$end omp parallel do

  WRITE(2,fmt="(I10, 3(3X,F12.8))") step, s

END SUBROUTINE COMPUTE_OP

subroutine compute_charge_non_smooth(dist_mat, charge)
  use parameters, only: nwater, nH
  implicit none
  integer                    :: iw, ih, min_ind
  real*8                     :: charge(nwater)
  real*8                     :: mindist, dist_mat(nH,nwater)

  charge = 0.d0
  do ih = 1,nH
    mindist=100.d0
    do iw = 1,nwater
      if (dist_mat(ih,iw)<mindist) then
        mindist=dist_mat(ih,iw)
        min_ind=iw
      endif
    end do
    charge(min_ind) = charge(min_ind) + 1.d0
  end do

  charge = charge - 2.d0

end subroutine compute_charge_non_smooth

SUBROUTINE COMPUTE_DIST_MATRIX(mat)
  USE parameters, only : natoms, nwater, atype, nH
  IMPLICIT NONE
  REAL*8                     :: mat(nH,nwater)
  REAL*8                     :: Dist
  INTEGER                    :: iw, ih, indH, indO

  indH = 0
  DO ih=nwater+1,nwater+nH
    DO iw=1,nwater
        mat(ih-nwater,iw) = Dist(ih,iw)
    END DO
  END DO

END SUBROUTINE COMPUTE_DIST_MATRIX
