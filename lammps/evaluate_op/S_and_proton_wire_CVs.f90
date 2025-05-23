!Computes the PT CV. This is not a general code!!!
!Only works if atype is ordered with O1 ... ON H1 ... H2N

MODULE OP
  INTEGER,ALLOCATABLE :: min_ind(:,:)
  INTEGER             :: hind_wire
  REAL*8, ALLOCATABLE :: charge(:)
  REAL*8, ALLOCATABLE :: dist_mat(:,:) 
  REAL*8              :: total_charge, sd
END MODULE OP

PROGRAM Coord_number
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL COMPUTE_DIST_MATRIX
    CALL compute_charge_non_smooth
    CALL COMPUTE_OP
    CALL PROTON_WIRE_LENGTH
    CALL PRINT_CV(frame)
  END DO

  CLOSE(1); CLOSE(2)

END PROGRAM Coord_number

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nwater, pos, atype, &
                         nframes, nequil
  USE OP
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file, out_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)
  CALL getarg(3, out_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nwater, nframes, nequil
  CLOSE(1)

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(min_ind(2*nwater,2))
  ALLOCATE(charge(nwater))
  ALLOCATE(dist_mat(2*nwater,nwater))

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = out_file)

END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_OP
  USE parameters, ONLY : nwater
  USE OP, ONLY: charge, dist_mat, total_charge, sd
  IMPLICIT NONE
  REAL*8                     :: Dist
  INTEGER                    :: iw,jw

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

subroutine compute_charge_non_smooth
  use parameters, only: nwater, natoms
  use OP, only: min_ind, dist_mat, charge
  implicit none
  integer                    :: iw, ih
  integer                    :: index_mindist_oxygen

  charge = 0.d0
  do ih = 1,natoms-nwater
    min_ind(ih,1)=index_mindist_oxygen(ih,0)
    charge(min_ind(ih,1)) = charge(min_ind(ih,1)) + 1.d0
  end do

  charge = charge - 2.d0

end subroutine compute_charge_non_smooth

SUBROUTINE PROTON_WIRE_LENGTH
  USE parameters, only: natoms,nwater
  USE OP
  IMPLICIT NONE
  INTEGER                    :: ih,iw,minind
  INTEGER                    :: index_mindist_oxygen
  REAL*8                     :: minsum, sumh, chargeprod

  minsum=100.d0
  minind=0
  DO ih=1,natoms-nwater
    min_ind(ih,2)=index_mindist_oxygen(ih,min_ind(ih,1))
    chargeprod=charge(min_ind(ih,1))*charge(min_ind(ih,2))
    if ((total_charge<1).or.(chargeprod.eq.-1)) then
      sumh=dist_mat(ih,min_ind(ih,1))+dist_mat(ih,min_ind(ih,2))
      if (sumh<minsum) then
        minsum=sumh
        minind=ih
      end if
    end if
  END DO
 
  hind_wire=minind

  !print *, min_ind(minind,1), min_ind(minind,2), sumh

END SUBROUTINE PROTON_WIRE_LENGTH

SUBROUTINE COMPUTE_DIST_MATRIX
  USE parameters, only : natoms, nwater
  USE OP, ONLY :dist_mat
  IMPLICIT NONE
  REAL*8                     :: Dist
  INTEGER                    :: iw, ih

  DO ih=nwater+1,natoms
    DO iw=1,nwater
      dist_mat(ih-nwater,iw) = Dist(ih,iw)
    END DO
  END DO

END SUBROUTINE COMPUTE_DIST_MATRIX

SUBROUTINE PRINT_CV(step)
  USE OP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: step
  REAL*8              :: Dist

  WRITE(2,fmt='(I10,7(3X,F12.8),2(3X,I4))') step, sd, total_charge, &
        charge(min_ind(hind_wire,1)), charge(min_ind(hind_wire,2)), &
        dist_mat(hind_wire,min_ind(hind_wire,1)), &
        dist_mat(hind_wire,min_ind(hind_wire,2)), &
        Dist(min_ind(hind_wire,1),min_ind(hind_wire,2)), &
        min_ind(hind_wire,1), min_ind(hind_wire,2)
  
END SUBROUTINE PRINT_CV

INTEGER FUNCTION index_mindist_oxygen(ih,iw_exclude) 
  USE parameters, only: natoms,nwater
  USE OP, only: dist_mat
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ih, iw_exclude
  INTEGER :: iw, minind
  REAL*8 :: mindist
  
  mindist=100.d0
  do iw = 1,nwater
    if ((dist_mat(ih,iw)<mindist).and.(iw/=iw_exclude)) then
      !if ((total_charge>0).and.(charge(iw).ne.0)) then
        mindist=dist_mat(ih,iw)
        minind=iw
      !end if
    endif
  end do
  
  index_mindist_oxygen=minind

END FUNCTION index_mindist_oxygen
