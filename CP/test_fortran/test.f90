program test

complex*8 :: a(3,3)
real*8    :: b(3)
integer   :: i,j

do i =1,3
  b(i) = dble(i+2)
  do j=1,3
    a(i,j) = CMPLX(dble(i),dble(j))  
  end do
end do

print *, b
print *, a(:,1)
print *, sum(a(:,1) * b(:))
end program test
