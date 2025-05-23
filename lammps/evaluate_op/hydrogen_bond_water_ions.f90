!Computes the PT CV. This is not a general code!!!
!Only works if atype is ordered with O1 ... ON H1 ... H2N

MODULE OP
  REAL*8 :: lambda
  REAL*8, ALLOCATABLE :: charge(:)
  INTEGER,ALLOCATABLE :: min_ind(:)
END MODULE OP

PROGRAM Coord_number
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER :: frame
  REAL*8 :: S

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL COMPUTE_OP(S)
    CALL PRINT_HYDROGEN_BONDS(S,frame)
  END DO

  CLOSE(1); CLOSE(2)

END PROGRAM Coord_number

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nwater, pos, atype, &
                         nframes, nequil
  USE OP, ONLY : lambda, min_ind, charge
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
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(min_ind(2*nwater))
  ALLOCATE(charge(nwater))

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = out_file)

  WRITE(2,*) "# nhb_accepted(h3O) nhb_accepted(oh) nhb_donated(h3O) nhb_donated(oh)"

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_OP(sd)
  !Compute center of mass along the z direction
  !Here we use the particle with index 1 as a reference
  USE parameters, ONLY : pos, nwater, natoms
  USE OP, ONLY: lambda, charge
  IMPLICIT NONE
  REAL*8                     :: compute_charge
  REAL*8                     :: dist_mat(natoms-nwater,nwater), sd
  REAL*8                     :: Dist, Total_charge
  REAL*8, PARAMETER          :: alpha=0.0001 !Smoothness of total charge
  INTEGER                    :: iw,jw

  CALL COMPUTE_DIST_MATRIX(dist_mat)

  CALL compute_charge_non_smooth(dist_mat,charge)

  sd = 0.d0

  DO iw = 1,nwater
!$omp parallel do reduction(+:sd)
    DO jw = iw+1,nwater
      sd = sd - charge(iw)*charge(jw)*Dist(iw,jw)
    END DO
!$end omp parallel do
  END DO

  Total_charge = sum(abs(charge))

END SUBROUTINE COMPUTE_OP

subroutine compute_charge_non_smooth(dist_mat, charge)
  use parameters, only: nwater, natoms
  use OP, only: min_ind
  implicit none
  integer                    :: iw, ih
  real*8                     :: charge(nwater)
  real*8                     :: mindist, dist_mat(natoms-nwater,nwater)

  min_ind=0
  charge = 0.d0
  do ih = nwater+1,natoms
    mindist=100.d0
    do iw = 1,nwater
      if (dist_mat(ih-nwater,iw)<mindist) then
        mindist=dist_mat(ih-nwater,iw)
        min_ind(ih-nwater)=iw
      endif
    end do
    charge(min_ind(ih-nwater)) = charge(min_ind(ih-nwater)) + 1.d0
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

SUBROUTINE PRINT_HYDROGEN_BONDS(S,step)
  USE parameters, only: natoms, nwater, pos
  USE op, only: charge 
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: S
  INTEGER, INTENT(IN) :: step
  INTEGER :: ih, iw
  INTEGER :: nhb_accepted(2), nhb_donated(2)

  DO iw=1,nwater
    if(charge(iw)==1) then
      CALL NUMBER_HBOND(iw,nhb_accepted(1),nhb_donated(1))
    else if (charge(iw)==-1) then
      CALL NUMBER_HBOND(iw,nhb_accepted(2),nhb_donated(2))
    end if
  END DO

  if (any(charge==1)) then
    WRITE(2,fmt="(I8,3X,F12.8,3X,4(I3,3X))"), step, S, nhb_accepted, nhb_donated
  end if

END SUBROUTINE PRINT_HYDROGEN_BONDS

SUBROUTINE NUMBER_HBOND(indO,nhb_accepted,nhb_donated)
  USE parameters, only: pos, natoms, nwater
  USE OP, only: min_ind 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: indO
  INTEGER :: nhb_accepted, nhb_donated
  INTEGER :: iO, iH
  REAL*8, PARAMETER :: max_dOO = 3.5d0, max_angle=-0.8660 
  REAL*8 :: dOO, Dist, ang_hoh, Angle

  nhb_accepted=0
  nhb_donated =0

  DO iO=1,nwater
    dOO=Dist(indO,iO)
    IF ((iO.ne.indO).and.(dOO < max_dOO)) THEN
      DO iH=nwater+1,natoms
        IF (min_ind(iH-nwater).eq.iO) THEN
          ang_hoh=Angle(pos(:,iH),pos(:,iO),pos(:,indO))
          IF (ang_hoh<max_angle) nhb_accepted = nhb_accepted + 1
        ELSE IF (min_ind(iH-nwater).eq.indO) THEN
          ang_hoh=Angle(pos(:,iH),pos(:,iO),pos(:,indO))
          IF (ang_hoh<max_angle) nhb_donated = nhb_donated + 1
        END IF 
      END DO
    END IF
  END DO

END SUBROUTINE NUMBER_HBOND

