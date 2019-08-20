program runCYGNUS_Exposure
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax,i1,i2,i,n_ex,nthvals
	double precision :: m_min,m_max,sigma_min,sigma_max,Vol,Time,ex_min,ex_max,m,m1,m2
	double precision,dimension(:),allocatable :: ex_vals,DL_F,Nbg_F,Nsig_F,E_th_vals
	character(len=100) :: filename
	character(len=100) :: fn_end
	character(len=100) :: fname
	

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)


	!------------------Optimisation-----------------!
	!-----Set binning ------------------------------!
	  nE_bins = 20 ! Number of energy bins		  			
	!----- Cross section range for limits ----------!
	  ns = 1000 ! resolution of cs scan							
	  sigma_min = 1.0d-51 ! min. expected cs				
	  sigma_max = 1.0d-38 ! max expected cs					
	!-----------------------------------------------!
	  n_ex = 10
	  ex_min = 1.0*6.0*((5.0/20.0)*0.16/1000.0d0) ! Convert m^3-years into ton-years
	  ex_max = 1.0e10*6.0*((5.0/20.0)*0.16/1000.0d0) ! Convert m^3-years into ton-years
	!-----------------------------------------------!

	allocate(ex_vals(n_ex))
	allocate(DL_F(n_ex))
	allocate(Nsig_F(n_ex))
	allocate(Nbg_F(n_ex))
	nucleus = Fluorine
	lat = Boulby(1)
	long = Boulby(2)
	searchmode = .false.
	ex_vals = logspace(ex_min,ex_max,n_ex)


	m = 9.0
	nthvals = 9
	allocate(E_th_vals(nthvals))
	E_th_vals = (/0.25d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0/)
	
	
!
! 	write(*,*) '----------------------------------------------------'
! 	filename = '../data/CYGNUS-Exposure9GeV-Predrift.txt'
! 	open(unit=1000,file=trim(filename))
! 	write(1000,*) 0.0,0.0,ex_vals/(6.0*(5.0/20.0)*(0.16/1000.0)) ! Volume (F)
! 	readout = 1
! 	nside = 4
!   	nT_bins = 10
! 	call LoadReadout(readout,	fn_end)
! 	write(*,*) 'Readout = ',fn_end
! 	efficiency_on = .true.
! 	energyres_on = .true.
! 	angres_on = .true.
! 	headtail_on = .true.
! 	energy_on = .true. ! energy info is currently turned on for best limits
! 	do i = 1,nthvals
! 		E_th = E_th_vals(i)
! 		call LoadReadout(readout,	fn_end)
! 		call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
! 		write(1000,*) Nbg_F(1),Nbg_F(n_ex),DL_F
! 	end do
! 	close(1000)
! 	write(*,*) '----------------------------------------------------'
!
	
	
	write(*,*) '----------------------------------------------------'
	filename = '../data/CYGNUS-Exposure9GeV-Postdrift.txt'
	open(unit=1000,file=trim(filename))
	write(1000,*) 0.0,0.0,ex_vals/(6.0*(5.0/20.0)*(0.16/1000.0)) ! Volume (F)
	readout = 2
	nside = 4
  	nT_bins = 10
	call LoadReadout(readout,	fn_end)
	write(*,*) 'Readout = ',fn_end
	efficiency_on = .true.
	energyres_on = .true.
	angres_on = .true.
	headtail_on = .true.
	energy_on = .true. ! energy info is currently turned on for best limits
	do i = 1,nthvals
		E_th = E_th_vals(i)
		call LoadReadout(readout,	fn_end)
		call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
		write(1000,*) Nbg_F(1),Nbg_F(n_ex),DL_F
	end do
	close(1000)
	write(*,*) '----------------------------------------------------'
	
	
	
	
	write(*,*) '----------------------------------------------------'
	filename = '../data/CYGNUS-Exposure9GeV-Predrift-Nondirectional.txt'
	open(unit=1000,file=trim(filename))
	write(1000,*) 0.0,0.0,ex_vals/(6.0*(5.0/20.0)*(0.16/1000.0)) ! Volume (F)
	readout = 1
	nside = 0
  	nT_bins = 10
	call LoadReadout(readout,	fn_end)
	write(*,*) 'Readout = ',fn_end
	efficiency_on = .true.
	energyres_on = .true.
	angres_on = .true.
	headtail_on = .true.
	energy_on = .true. ! energy info is currently turned on for best limits
	do i = 1,nthvals
		E_th = E_th_vals(i)
		call LoadReadout(readout,	fn_end)
		call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
		write(1000,*) Nbg_F(1),Nbg_F(n_ex),DL_F
	end do
	close(1000)
	write(*,*) '----------------------------------------------------'
	
	
	write(*,*) '----------------------------------------------------'
	filename = '../data/CYGNUS-Exposure9GeV-Predrift-NoEnergy.txt'
	open(unit=1000,file=trim(filename))
	write(1000,*) 0.0,0.0,ex_vals/(6.0*(5.0/20.0)*(0.16/1000.0)) ! Volume (F)
	readout = 1
	nside = 0
  	nT_bins = 10
	call LoadReadout(readout,	fn_end)
	write(*,*) 'Readout = ',fn_end
	efficiency_on = .true.
	energyres_on = .true.
	angres_on = .true.
	headtail_on = .true.
	energy_on = .false. ! energy info is currently turned on for best limits
	do i = 1,nthvals
		E_th = E_th_vals(i)
		call LoadReadout(readout,	fn_end)
		call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
		write(1000,*) Nbg_F(1),Nbg_F(n_ex),DL_F
	end do
	close(1000)
	write(*,*) '----------------------------------------------------'
	
	call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start



	stop

	write(*,*) '----------------------------------------------------'
	readout = 2
	nside = 0
  	nT_bins = 50	 ! Number of time bins					 
	call LoadReadout(readout,	fn_end)
	write(*,*) 'Readout = ',fn_end
	efficiency_on = .true.
	energyres_on = .true.
	angres_on = .false.
	headtail_on = .false.
	energy_on = .false.
	filename = '../data/CYGNUS-Exposure-NoEnergy.txt'
	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(6.0*(5.0/20.0)*(0.16/1000.0)) ! Volume (F)
	m = m1
	nucleus = Fluorine
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
	write(1000,*) DL_F
	m = m2
	nucleus = Helium
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
	write(1000,*) DL_F
	write(1000,*) Nbg_F
	close(1000)
					 
	
	write(*,*) '----------------------------------------------------'
	filename = '../data/CYGNUS-Exposure9GeV-Nondirectional.txt'
	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(6.0*(5.0/20.0)*(0.16/1000.0)) ! Volume (F)
	readout = 2
	nside = 0
  	nT_bins = 1
	call LoadReadout(readout,	fn_end)
	write(*,*) 'Readout = ',fn_end
	efficiency_on = .true.
	energyres_on = .true.
	angres_on = .true.
	headtail_on = .true.
	energy_on = .true. ! energy info is currently turned on for best limits

	do i = 1,nthvals
		E_th = E_th_vals(i)
		call LoadReadout(readout,	fn_end)
		call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
		write(1000,*) DL_F
	end do
	close(1000)
	write(*,*) '----------------------------------------------------'
	
	
	stop
	
	
	
	write(*,*) '----------------------------------------------------'
	readout = 2
	nside = 0
	nT_bins = 50	 
	call LoadReadout(readout,	fn_end)
	write(*,*) 'Readout = ',fn_end
	efficiency_on = .true.
	energyres_on = .true.
	angres_on = .false.
	headtail_on = .false.
	energy_on = .true. ! energy info is currently turned on for best limits
	filename = '../data/CYGNUS-Exposure-Nondirectional.txt'
	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(6.0*(5.0/20.0)*(0.16/1000.0)) ! Volume (F)
	m = m1
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
	write(1000,*) DL_F
	m = m2
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)
	write(1000,*) DL_F
	write(1000,*) Nbg_F
	close(1000)
	
	
	
	

	
	
	
	
	call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start

end program runCYGNUS_Exposure
