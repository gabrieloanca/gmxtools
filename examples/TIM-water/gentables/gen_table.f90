!fortran code to generate gromacs tables
!these tables contain the functions f(x), -f'(x), g(x), -g'(x), h(x), and -h'(x) from eq:
!V = q1*q1/(4*pi*e0)*f(r) + C*g(r) + A*h(r)
!where normally f(r) =  1/r,       -f'(r) =  1/(r**2)   - Coulomb
!               g(r) = -1/(r**6),  -g'(r) = -6/(r**7)   - dipersion
!               h(r) =  1/(r**12), -h'(r) =  12/(r**13) - short range repulsion
!            or h(r) =  exp(-b*r), -h'(r) = b/exp(-b*r) - soft-core repulsion
!
!
! For a guide of how to use table, see the following link:
! https://www.gromacs.org/Documentation_of_outdated_versions/How-tos/Tabulated_Potentials
! and the pdf there in:
! https://www.gromacs.org/@api/deki/files/94/=gromacs_nb.pdf
!
! For a guide of how to compile and run fortran code, see:
! https://gcc.gnu.org/wiki/GFortranGettingStarted
! To compile this code:
!   gfortran gen_table.f90 -o gen_table.out
! Then run it as:
!   ./gen_table.out > table_r1_hr.xvg
! then change the values at r=0 (the infinities) to some arbitrary values.
!
!! It is recommended that the spacing between the adjacent rs in your table should equal
!! 0.002 nm or 0.0005 nm for the single or double precision versions of gromacs respectively.
!
!
program gen_table
implicit none
real, parameter :: delr=0.002, rcut=4.0
real :: r
integer :: nbins, j

nbins=int((rcut+1)/delr) + 1

do j=0,nbins
    r=delr*j
    !!! NB tables
    !write(6,*) r, 1/r, 1/(r*r), -1/(r**6), -6/(r**7), exp(-25.0*r), 25.0*exp(-25.0*r)
    ! Coul and LJ6-12
    !write(6,*) r, 1/r, 1/(r*r), -1/(r**6), -6/(r**7), 1/(r**12), 12/(r**13)
    ! LJ12
    !write(6,*) r, -1/(r**6), -6/(r**7)
    !!! Bond tables
    !soft-core, b=25.0
    write(6,*) r, exp(-28.0*r), 28.0*exp(-28.0*r)
    !! electrostatic
    !write(6,*) r, 1/r, 1/(r*r)
    !! Morse, a=20.0, b0=0.1
    !! f'(x) = 2 * beta * x0 * [1-exp(-beta*(x-x0))] * exp(-beta*(x-x0))
    !write(6,*) r, (1-exp(-20.0*(r-0.01)))**2, -4.0*exp(-20.0*(r-0.1))*(1-exp(-20.0*(r-0.1)))
    !! Harmonic
    !write(6,*) r, (r-0.09572)**2, -2*(r-0.09572)
    !write(6,*) r, (r-0.1)**2, -2*(r-0.1)
end do

end program

