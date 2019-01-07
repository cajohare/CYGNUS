program runCYGNUS_thresholds
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax
	double precision :: m_min,m_max,sigma_min,sigma_max,Vol,Time
	double precision,dimension(:),allocatable :: m_vals,DLHe,DLF
	character(len=100) :: filename
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)


	!------------------Optimisation-----------------!
	!-----Set binning ------------------------------!
	  nE_bins = 20 ! Number of energy bins		  		!
	  nT_bins = 1	 ! Number of time bins					  !
	  nside = 4  ! Order of pixelation (2,4 or 8)		!
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
	allocate(m_vals(nm))
	allocate(DLF(nm))
	allocate(DLHe(nm))

	VolTime = 1000.0*3.0d0
	Exposure = VolTime*(0.16/1000.0d0) ! Convert m^3-years into ton-years
	lat = Boulby(1)
	long = Boulby(2)
	E_max = 200.0d0

	readout = 2
	energy_on = .true. ! energy info is currently turned on for best limits
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .true.
	energyres_on = .false.
	searchmode = .false.
	call LoadReadout(readout,	fn_end)



	write(*,*) '----------------------------------------------------'
	filename = '../data/CYGNUS1000-1keV3keV.txt'
	open(unit=123,file=trim(filename))

	nucleus = Helium
	E_th = 1.8d0
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)

	nucleus = Fluorine
	E_th = 3.0d0
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)

	write(123,*) m_vals
	write(123,*) DLF
	write(123,*) DLHe
	close(123)
	write(*,*) 'writing to: ',trim(filename)
	write(*,*) '----------------------------------------------------'



	write(*,*) '----------------------------------------------------'
	filename = '../data/CYGNUS1000-3keV5keV.txt'
	open(unit=123,file=trim(filename))

	nucleus = Helium
	E_th = 3.0d0
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)

	nucleus = Fluorine
	E_th = 5.0d0
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)

	write(123,*) m_vals
	write(123,*) DLF
	write(123,*) DLHe
	close(123)
	write(*,*) 'writing to: ',trim(filename)
	write(*,*) '----------------------------------------------------'


	write(*,*) '----------------------------------------------------'
	filename = '../data/CYGNUS1000-6keV10keV.txt'
	open(unit=123,file=trim(filename))

	nucleus = Helium
	E_th = 6.0d0
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)

	nucleus = Fluorine
	E_th = 10.0d0
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)

	write(123,*) m_vals
	write(123,*) DLF
	write(123,*) DLHe
	close(123)
	write(*,*) 'writing to: ',trim(filename)
	write(*,*) '----------------------------------------------------'


end program runCYGNUS_thresholds
