program runCYGNUS_Ideal
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax
	double precision :: m_min,m_max,sigma_min,sigma_max,Vol,Time
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

	!-----------------------------------------------!
	! CYGNUS 1000 x 3 years of running
	! 30 degrees angular res.
	nside = 0
	VolTime = 1000.0*3.0d0
	filename = '../data/CYGNUS1000-Nondirectional.txt'
	readout = 2
	energy_on = .true. ! energy info is currently turned on for best limits
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .false.
	energyres_on = .false.
	searchmode = .false.
	call LoadReadout(readout,	fn_end)
	call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)

stop
	nside = 4

	!-----------------------------------------------!
	! CYGNUS 1000 x 3 years of running
	! 30 degrees angular res.
	VolTime = 1000.0*3.0d0
	filename = '../data/CYGNUS1000-30deg.txt'
	readout = 2
	energy_on = .true. ! energy info is currently turned on for best limits
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .false.
	energyres_on = .false.
	searchmode = .false.
	call LoadReadout(readout,	fn_end)
	call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)

	! 30 degrees angular res. + head-tail
	filename = '../data/CYGNUS1000-30deg-HT.txt'
	headtail_on = .true.
	call LoadReadout(readout,	fn_end)
	call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)



stop


	readout = 0
	energy_on = .true. ! energy info is currently turned on for best limits
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .true.
	energyres_on = .false.
	searchmode = .false.
	call LoadReadout(readout,	fn_end)


	! NORMAL OPERATION
	!-----------------------------------------------!
	! CYGNUS 1000 x 3 years of running
	VolTime = 1000.0*3.0d0
	filename = '../data/CYGNUS1000-Final.txt'
	call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)

	! CYGNUS 100k x 3 years of running
	VolTime = 100000.0*3.0d0
	filename = '../data/CYGNUS100k-Final.txt'
	call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	!-----------------------------------------------!


	! SEARCH MODE OPERATION
	searchmode = .true.
	nside = 0
	!-----------------------------------------------!
	VolTime = 1000.0*3.0d0
	filename = '../data/CYGNUS1000-SearchMode.txt'
	call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)

	! CYGNUS 100k x 3 years of running
	VolTime = 100000.0*3.0d0
	filename = '../data/CYGNUS100k-SearchMode.txt'
	call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	!-----------------------------------------------!


end program runCYGNUS_Ideal
