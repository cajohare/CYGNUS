program runCYGNUS_Exposure
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax,i1,i2,i,n_ex
	double precision :: m_min,m_max,sigma_min,sigma_max,Vol,Time,ex_min,ex_max,m
	double precision,dimension(:),allocatable :: ex_vals,DL
	character(len=100) :: filename
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)


	!------------------Optimisation-----------------!
	!-----Set binning ------------------------------!
	  nE_bins = 5 ! Number of energy bins		  		!
	  nT_bins = 1	 ! Number of time bins					  !
	  nside = 16  ! Order of pixelation (2,4 or 8)		!
	!----- Mass range for limits -------------------!
	  nm = 200 !  Number of mass points							!
	  m_min = 0.5d0 ! Min mass											!
	  m_max = 10000.0d0 ! Max mass										!
	!----- Cross section range for limits ----------!
	  ns = 500 ! resolution of cs scan							!
	  sigma_min = 1.0d-50 ! min. expected cs				!
	  sigma_max = 1.0d-41 ! max expected cs					!
	!-----------------------------------------------!
	!-----------------------------------------------!

	n_ex = 20
	ex_min = 0.01*3.0*(0.16/1000.0d0) ! Convert m^3-years into ton-years
	ex_max = 1.0e10*3.0*(0.16/1000.0d0) ! Convert m^3-years into ton-years
	allocate(ex_vals(n_ex))
	allocate(DL(n_ex))


	lat = Boulby(1)
	long = Boulby(2)
	E_max = 200.0d0
	E_th = 3.0d0
	nucleus = Fluorine
	readout = 2
	energy_on = .true. ! energy info is currently turned on for best limits
	energyres_on = .false.
	searchmode = .false.
	efficiency_on = .true.


	m = 10.0

	filename = '../data/CYGNUS1000-Exposure-Ideal.txt'
	angres_on = .false.
	headtail_on = .false.
	call LoadReadout(readout,	fn_end)
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL)
	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(3*(0.16/1000.0))
	write(1000,*) DL
	close(1000)

	filename = '../data/CYGNUS1000-Exposure-angres.txt'
	angres_on = .true.
	headtail_on = .false.
	call LoadReadout(readout,	fn_end)
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL)
	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(3*(0.16/1000.0))
	write(1000,*) DL
	close(1000)

	filename = '../data/CYGNUS1000-Exposure-angres-HT.txt'
	angres_on = .true.
	headtail_on = .true.
	call LoadReadout(readout,	fn_end)
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL)
	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(3*(0.16/1000.0))
	write(1000,*) DL
	close(1000)

	nside = 0
	filename = '../data/CYGNUS1000-Exposure-Nondirectional.txt'
	angres_on = .true.
	headtail_on = .true.
	call LoadReadout(readout,	fn_end)
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL)
	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(3*(0.16/1000.0))
	write(1000,*) DL
	close(1000)

end program runCYGNUS_Exposure
