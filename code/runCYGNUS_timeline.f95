program runCYGNUS_thresholds
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer,parameter :: n_expts = 6
	integer :: i,j,jj,ns,nm,readout_selection,loopmin,loopmax,n_ex
	double precision :: m_min,m_max,sigma_min,sigma_max,ex_min,ex_max,Time
	double precision,dimension(:),allocatable :: m_vals,DLHe,DLF,DL
	double precision :: E_th_F(n_expts),E_th_He(n_expts),Vol(n_expts)
	character(len=100) :: filename
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)

	
	
	!------------------Optimisation-------------!
	!-----Set binning --------------------------!
	nE_bins = 20 ! Number of energy bins	
	nT_bins = 1	 ! Number of time bins		
	nside = 4  ! Order of pixelation (2,4 or 8)	
	!----- Mass range for limits ---------------!
	nm = 200 !  Number of mass points	
	m_min = 0.1d0 ! Min mass									
	m_max = 10000.0d0 ! Max mass							
	!----- Cross section range for limits ------!
	ns = 500 ! resolution of cs scan			
	sigma_min = 1.0d-50 ! min. expected cs				
	sigma_max = 1.0d-38 ! max expected cs	


	allocate(m_vals(nm))
	allocate(DLF(nm))
	allocate(DLHe(nm))
	allocate(DL(nm))


	Vol = (/1.0,10.0,100.0,1000.0,10000.0,100000.0/)
	E_th_He = (/0.25,0.5,1.0,1.0,1.0,2.0/)
	E_th_F = (/0.5,0.75,1.0,1.5,2.0,3.0/)
	Time = 6.0
	lat = Boulby(1)
	long = Boulby(2)
	

	nside = 4
	readout = 1 ! Predrift
	energy_on = .true. ! energy info is currently turned on for best limits
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .true.
	energyres_on = .false.
	searchmode = .false.
	call LoadReadout(readout,	fn_end)
	
	
	filename = '../data/CYGNUS10-timeline.txt'
	open(unit=123,file=trim(filename))
	write(123,*) m_vals
	do j = 1,n_expts
		write(*,*) '----------------------------------------------------'
		
		Exposure = Vol(j)*Time*(0.16/1000.0d0)*(755.0/740.0)
		nucleus = Helium
		E_th = E_th_He(j)
		call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)
		
		Exposure = Vol(j)*Time*(0.16/1000.0d0)*(5.0/20.0)
		nucleus = Fluorine
		E_th = E_th_F(j)
		call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)

		do i = 1,nm
			if (DLF(i)==0.0) then
				DL(i) = DLHe(i)
			elseif (DLHe(i)==0.0) then
				DL(i) = DLF(i)
			else
				DL(i) = min(DLF(i),DLHe(i))
			end if
		end do
		write(123,*) Vol,DL
		
	end do
	close(123)


end program runCYGNUS_thresholds
