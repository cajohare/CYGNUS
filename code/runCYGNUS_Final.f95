program runCYGNUS_Ideal
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer,parameter :: n_expts = 6
	integer :: i,jj,ns,nm,readout_selection,loopmin,loopmax
	double precision :: m_min,m_max,sigma_min,sigma_max,Time
	double precision :: E_th_F(n_expts),E_th_He(n_expts),Vol(n_expts)
	logical :: saveNwimp
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


	! POSTDRIFT (NORMAL MODE)
	readout = 2
	energy_on = .true.
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .true.
	energyres_on = .true.
	searchmode = .false.
	saveNwimp = .false.
	Time = 6.0 ! Years
	


	Vol = (/1.0,10.0,100.0,1000.0,10000.0,100000.0/)
	E_th_He = (/0.25,1.0,3.0,8.0,8.0,8.0/)
	E_th_F = (/0.25,1.0,3.0,8.0,8.0,8.0/)

	i = 1
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS1-Postdrift_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)
	
	i = 2
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS10-Postdrift_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)
	
	i = 3
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS100-Postdrift_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)
	
	i = 4
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS1000-Postdrift_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)
	
	i = 5
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS10k-Postdrift_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)

	i = 6
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS100k-Postdrift_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)
	
	
	
	
	!!!!!!!!!!!
		
	
	! NONDIRECTIONAL COMPARISON
	nside = 0
	
	

	i = 1
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS1-Nondirectional_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)
	
	i = 2
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS10-Nondirectional_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)
	
	i = 3
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS100-Nondirectional_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)
	
	i = 4
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS1000-Nondirectional_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)
	
	i = 5
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS10k-Nondirectional_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)

	i = 6
	call LoadReadout(readout,	fn_end)
	VolTime = Vol(i)*Time
	filename = '../data/CYGNUS100k-Nondirectional_Final.txt'
	call CYGNUSLimit(E_th_F(i),E_th_He(i),m_min,m_max,nm,sigma_min,sigma_max,ns,filename,saveNwimp)
	
	


end program runCYGNUS_Ideal
