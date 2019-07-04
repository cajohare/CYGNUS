program runCYGNUS_Ideal
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax
	double precision :: m_min,m_max,sigma_min,sigma_max,Vol,Time,E_th_F,E_th_He
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
	  nm = 50 !  Number of mass points							!
	  m_min = 0.2d0 ! Min mass											!
	  m_max = 10000.0d0 ! Max mass										!
	!----- Cross section range for limits ----------!
	  ns = 100 ! resolution of cs scan							!
	  sigma_min = 1.0d-50 ! min. expected cs				!
	  sigma_max = 1.0d-39 ! max expected cs					!
	!-----------------------------------------------!
	!-----------------------------------------------!

	! ! NON-DIRECTIONAL
	! nside = 0
	! VolTime = 1000.0*3.0d0
	! filename = '../data/CYGNUS1000-Nondirectional.txt'
	! readout = 2
	! energy_on = .true. ! energy info is currently turned on for best limits
	! angres_on = .true.
	! efficiency_on = .true.
	! headtail_on = .false.
	! energyres_on = .false.
	! searchmode = .false.
	! call LoadReadout(readout,	fn_end)
	! call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	!
	! ! JUST ANGULAR RES.
	! nside = 4
	! VolTime = 1000.0*3.0d0
	! filename = '../data/CYGNUS1000-30deg.txt'
	! readout = 2
	! energy_on = .true. ! energy info is currently turned on for best limits
	! angres_on = .true.
	! efficiency_on = .true.
	! headtail_on = .false.
	! energyres_on = .false.
	! searchmode = .false.
	! call LoadReadout(readout,	fn_end)
	! call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	!
	! ! ANGULAR RES + HEAD-TAIL
	! filename = '../data/CYGNUS1000-30deg-HT.txt'
	! headtail_on = .true.
	! call LoadReadout(readout,	fn_end)
	! call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)

	! PREDRIFT (NORMAL MODE)
	readout = 1
	energy_on = .true.
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .true.
	energyres_on = .true.
	searchmode = .false.
	E_th_F = 2.5
	E_th_He = 1.8
	Time = 6.0 ! Years

	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d1*Time
	filename = '../data/CYGNUS10-predrift_Final.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d2*Time
	filename = '../data/CYGNUS100-predrift_Final.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d3*Time
	filename = '../data/CYGNUS1000-predrift_Final.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d4*Time
	filename = '../data/CYGNUS10k-predrift_Final.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d5*Time
	filename = '../data/CYGNUS100k-predrift_Final.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	



	! PREDRIFT SEARCH MODE OPERATION
	searchmode = .true.
	nside = 0
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d1*Time
	filename = '../data/CYGNUS10-SearchMode-predrift.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d2*Time
	filename = '../data/CYGNUS100-SearchMode-predrift.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d3*Time
	filename = '../data/CYGNUS1000-SearchMode-predrift.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d4*Time
	filename = '../data/CYGNUS10k-SearchMode-predrift.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d5*Time
	filename = '../data/CYGNUS100k-SearchMode-predrift.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	
	
	! NONDIRECTIONAL COMPARISON
	searchmode = .false.

	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d1*Time
	filename = '../data/CYGNUS10-NonDirectional_Final.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d2*Time
	filename = '../data/CYGNUS100-NonDirectional_Final.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d3*Time
	filename = '../data/CYGNUS1000-NonDirectional_Final.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d4*Time
	filename = '../data/CYGNUS10k-NonDirectional_Final.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	call LoadReadout(readout,	fn_end)
	VolTime = 1.0d5*Time
	filename = '../data/CYGNUS100k-NonDirectional_Final.txt'
	call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
	
	


end program runCYGNUS_Ideal
