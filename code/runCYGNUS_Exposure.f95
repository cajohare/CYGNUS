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
	double precision,dimension(:),allocatable :: ex_vals,DL_F,Nbg_F,Nsig_F
	character(len=100) :: filename
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)


	!------------------Optimisation-----------------!
	!-----Set binning ------------------------------!
	  nE_bins = 20 ! Number of energy bins		  		
	  nT_bins = 1	 ! Number of time bins					 
	  nside = 4  ! Order of pixelation (2,4 or 8)		
	!----- Cross section range for limits ----------!
	  ns = 1000 ! resolution of cs scan							
	  sigma_min = 1.0d-51 ! min. expected cs				
	  sigma_max = 1.0d-38 ! max expected cs					
	!-----------------------------------------------!
	  n_ex = 20
	  ex_min = 0.001*6.0*((5.0/20.0)*0.16/1000.0d0) ! Convert m^3-years into ton-years
	  ex_max = 1.0e10*6.0*((5.0/20.0)*0.16/1000.0d0) ! Convert m^3-years into ton-years
	!-----------------------------------------------!

	allocate(ex_vals(n_ex))
	allocate(DL_F(n_ex))
	allocate(Nsig_F(n_ex))
	allocate(Nbg_F(n_ex))
	E_th = 2.5d0
	nucleus = Fluorine
	lat = Boulby(1)
	long = Boulby(2)
	energy_on = .true. ! energy info is currently turned on for best limits
	searchmode = .false.
	ex_vals = logspace(ex_min,ex_max,n_ex)

	write(*,*) '----------------------------------------------------'
	readout = 2
	call LoadReadout(readout,	fn_end)
	write(*,*) 'Readout = ',fn_end
	efficiency_on = .true.
	energyres_on = .true.
	angres_on = .true.
	headtail_on = .true.

	filename = '../data/CYGNUS-Exposure-Postdrift.txt'
	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(6.0*(5.0/20.0)*(0.16/1000.0)) ! Volume (F)

	m = 10.0
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
	write(1000,*) DL_F
	write(1000,*) Nsig_F
	write(1000,*) Nbg_F

	m = 100.0
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
	write(1000,*) DL_F
	write(1000,*) Nsig_F
	write(1000,*) Nbg_F
	close(1000)
	
	
	
	write(*,*) '----------------------------------------------------'
	readout = 2
	nside = 0
	call LoadReadout(readout,	fn_end)
	write(*,*) 'Readout = ',fn_end
	efficiency_on = .true.
	energyres_on = .true.
	angres_on = .false.
	headtail_on = .false.

	filename = '../data/CYGNUS-Exposure-Nondirectional.txt'
	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(6.0*(5.0/20.0)*(0.16/1000.0)) ! Volume (F)

	m = 10.0
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
	write(1000,*) DL_F
	write(1000,*) Nsig_F
	write(1000,*) Nbg_F

	m = 100.0
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
	write(1000,*) DL_F
	write(1000,*) Nsig_F
	write(1000,*) Nbg_F
	close(1000)
	
	
	call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start

end program runCYGNUS_Exposure
