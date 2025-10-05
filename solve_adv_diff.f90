! -----------------------------------------------------------------------------
! Project: Fractal Growth DBM
! File: fsolve_adv_diff.f90
! Author: Fabien Chauvet
! License: MIT (see LICENSE file for details)
! Copyright (c) 2025 Fabien Chauvet
! -----------------------------------------------------------------------------


subroutine solve_adv_diff(Ld,AA,cin,pl,n,m,i0,iff,j0,jff,erreur_conv,c)

!implicit none
double precision, intent(in) :: Ld
integer, intent(in) :: n
integer, intent(in) :: m
integer, intent(in) :: i0
integer, intent(in) :: iff
integer, intent(in) :: j0
integer, intent(in) :: jff
double precision, intent(in) :: erreur_conv
double precision, dimension (0:n-1,0:m-1), intent(in) :: cin
!f2py depend(n,m), cin 
double precision, dimension (0:n-1,0:m-1), intent(in) :: pl
!f2py depend(n,m), pl 
double precision, dimension (0:n-1,0:m-1), intent(in) :: AA
!f2py depend(n,m), AA
double precision, dimension (0:n-1,0:m-1), intent(out) :: c
double precision :: residu
integer :: k
double precision, dimension (0:n-1,0:m-1) :: cc 
double precision, dimension (0:iff-i0,0:jff-j0) :: c_sud, c_nord, c_est, c_ouest
double precision, dimension (0:iff-i0) :: c_sud_0, c_nord_0, c_est_0, c_ouest_0, c_sud_1, c_nord_1, c_est_1, c_ouest_1

residu = 1
k=0
c=cin
cc=c

do while (residu>erreur_conv)
	k=k+1
	cc=c

	! periodic boundary condition on lateral walls	
	c_sud_0=c(i0+1:iff+1,0)
	c_nord_0=c(i0-1:iff-1,0)
	c_est_0=c(i0:iff,1)
	c_ouest_0=c(i0:iff,m-1)

	c_sud_1=c(i0+1:iff+1,m-1)
	c_nord_1=c(i0-1:iff-1,m-1)
	c_est_1=c(i0:iff,0)
	c_ouest_1=c(i0:iff,m-2)

	c(i0:iff,0)=AA(i0:iff,0)*((c_sud_0+c_nord_0+ &
				& c_est_0+c_ouest_0)/4+ &
				& (1/(8*Ld))*(c_nord_0-c_sud_0))

	c(i0:iff,m-1)=AA(i0:iff,m-1)*((c_sud_1+c_nord_1+ &
				& c_est_1+c_ouest_1)/4+ &
				& (1/(8*Ld))*(c_nord_1-c_sud_1))

	!interior of the domain
	c_sud=c(i0+1:iff+1,j0:jff)
	c_nord=c(i0-1:iff-1,j0:jff)
	c_est=c(i0:iff,j0+1:jff+1)
	c_ouest=c(i0:iff,j0-1:jff-1)

	c(i0:iff+1-1,j0:jff+1-1)=AA(i0:iff+1-1,j0:jff+1-1)*((c_sud+c_nord+c_est+c_ouest)/4+(c_nord-c_sud)/(8*Ld))
!	c(i0:iff+1,j0:jff+1)=AA(i0:iff+1,j0:jff+1)*((c(i0+1:iff+2,j0:jff+1)+ &
!				& c(i0-1:iff,j0:jff+1) + &
!				& c(i0:iff+1,j0+1:jff+2) + &
!				& c(i0:iff+1,j0-1:jff))/4 + &
!				& (c(i0-1:iff,j0:jff+1)-c(i0+1:iff+2,j0:jff+1))/(8*Ld))

	c=c*AA

	residu = maxval(((c-cc)**2*pl))**0.5
!	print*, residu


end do

print*, 'N iter=', k

end subroutine solve_adv_diff
