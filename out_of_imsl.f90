module out_of_imsl
contains


!-----------------------------------------------inverse
subroutine solve_inverse(a,a_inverse,n)
	real(8)::a(n,n),a_inverse(n,n)
	integer::n
    real(8)::E(n,n)
    integer::i
	E=0
		do i=1,n
			E(i,i)=1
		end do
	call mateq(a,E,a_inverse,n,n)
	return
end subroutine solve_inverse


subroutine mateq(a,b,x,n,m)
implicit real*8(a-z)
integer::i,k,n,m
integer::id_max
real(8)::a(n,n),b(n,m),x(n,m)
real(8)::aup(n,n),bup(n,m)
real(8)::ab(n,n+m)
real(8)::vtemp1(n+m),vtemp2(n+m)
real(8)::vtmp(n),xtmp(n)
ab(1:n,1:n)=a
ab(:,n+1:n+m)=b
do k=1,n-1
	elmax=dabs(ab(k,k))
	id_max=k
	do i=k+1,n
	if(dabs(ab(i,k))>elmax) then
		elmax=ab(i,k)
        id_max=i
	end if
	end do
vtemp1=ab(k,:)
vtemp2=ab(id_max,:)
ab(k,:)=vtemp2
ab(id_max,:)=vtemp1
do i=k+1,n
	temp=ab(i,k)/ab(k,k)
	ab(i,:)=ab(i,:)-temp*ab(k,:)
	end do
end do
aup(:,:)=ab(1:n,1:n)
do i=1,m
vtmp=ab(:,n+i)
call uptri(aup,vtmp,xtmp,n)
x(:,i)=xtmp
end do
end subroutine mateq

subroutine uptri(a,b,x,n)
implicit real*8(a-z)
integer::i,j,n
real(8)::a(n,n),b(n),x(n)
x(n)=b(n)/a(n,n)
do i=n-1,1,-1
	x(i)=b(i)
	do j=i+1,n
		x(i)=x(i)-a(i,j)*x(j)
	end do
		x(i)=x(i)/a(i,i)
	end do
end subroutine uptri
!-----------------------------------------------inverse
!-----------------------------------------------diag
subroutine diag_extract(a,diag_vector,n)
integer::n
real(8)::a(n,n)
real(8)::diag_vector(n)
integer::i
	do i=1,n
		diag_vector(i)=a(i,i)
	end do
return
end subroutine diag_extract

subroutine diag_enchase(b,diag_vector,n)
integer::n
real(8)::b(n,n)
real(8)::diag_vector(n)
integer::i
	do i=1,n
		b(i,i)=diag_vector(i)
	end do
return
end subroutine diag_enchase
!----------------------------------------------diag
!----------------------------------------------QR_decompose
subroutine QR_decompose(a,Q,R,m,n)
implicit real*8(a-z)
integer::m,n
real(8)::a(m,n)
real(8)::Q(m,m),R(m,n)
real(8)::H0(m,m),H1(m,m),H2(m,m),qt(m,m)
real(8)::A1(m,n),A2(m,n),u(m)

integer::i,j,k
A1=a
H1=0d0
do j=1,m
	H1(j,j)=1d0
end do
do k=1,n
	h0=0d0
do i=1,m
	h0(i,i)=1d0
end do
	s=0d0
do i=k,m
	s=s+a1(i,k)*a1(i,k)
end do
	s=dsqrt(s)
	u=0d0
if(a1(k,k)>=0) then
	u(k)=a1(k,k)+s
else
	u(k)=a1(k,k)-s
end if

do i=k+1,m
	u(i)=a1(i,k)
end do
du=0
do i=k,m
du=du+u(i)*u(i)
end do
do i=k,m
	do j=k,m
		H0(i,j)=-2d0*u(i)*u(j)/du
		if(i==j) then
		H0(i,j)=1d0+H0(i,j)
		end if
	end do
end do
A2=matmul(H0,A1)
A1=A2
H1=matmul(H1,H0)
Q=H1
R=A1
end do
return	
end subroutine QR_decompose

subroutine eig_calculate(a,n,eig_value,tol)
implicit real*8(a-z)
integer::n
real(8)::a(n,n)
real(8)::eig_value(n)
real(8)::tol
integer::i,j,k
real(8)::a1(n,n),Q(n,n),R(n,n)
a1=a

do i=1,200
	call QR_decompose(a1,Q,R,n,n)
	a1=matmul(R,Q)
	do k=1,n
		ds=0d0
		ds=ds+A1(k,k)**2
	end do
	do j=1,n
		eig_value(j)=A1(j,j)
	end do
	if(ds<tol) exit
end do 


return
end subroutine eig_calculate
!----------------------------------------------------------
!----------------------------------------------------------det fucntion

recursive function det(a,n)
implicit real(a-z)
integer::n,i
real::a(n,n)
real::a1(n-1,n-1),a2(n-1,n-1),a3(n-1,n-1),a4(n-1,n-1)
	if(n==3)then
	det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
	else if(n==4)then 
		a1(:,:)=a(2:4,2:4)

		a2(:,1)=a(2:4,1)
        a2(:,2:3)=a(2:4,3:4)

		a3(:,1:2)=a(2:4,1:2)	
		a3(:,3)=a(2:4,4)	

		a4(:,:)=a(2:4,1:3)

	det=a(1,1)*det(a1,3)-a(1,2)*det(a2,3)+a(1,3)*det(a3,3)-a(1,4)*det(a4,3)
end if
return
end function det




end module




