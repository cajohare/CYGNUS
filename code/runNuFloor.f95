program runNuFloor
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax,i1,i2,i,n_ex,n_Eth
	double precision :: m_min,m_max,sigma_min,sigma_max,Vol,Time,ex_min,ex_max
	double precision,dimension(:),allocatable :: m_vals,ex_vals,Eth_vals
	double precision,dimension(:,:),allocatable :: DL
	double precision :: Eth_max,Eth_min
	character(len=100) :: filename
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)

	! Set all directional stuff to something, not relevant for this file
	lat = Boulby(1)
	long = Boulby(2)
	energy_on = .true. ! energy info is currently turned on for best limits
	searchmode = .false.
	angres_on = .false.
	efficiency_on = .false.
	headtail_on = .false.
	energyres_on = .false.
	call LoadReadout(9,	fn_end)


	!------------------Optimisation-----------------!
	!-----Set binning ------------------------------!
	  nE_bins = 20 ! Number of energy bins		  		!
	  nT_bins = 1	 ! Number of time bins					  !
	  nside = 0  ! Order of pixelation (2,4 or 8)		!
	!----- Mass range for limits -------------------!
	  nm = 100 !  Number of mass points							!
	  m_min = 0.5d0 ! Min mass											!
	  m_max = 100000.0d0 ! Max mass										!
	!----- Cross section range for limits ----------!
	  ns = 200 ! resolution of cs scan							!
	  sigma_min = 1.0d-50 ! min. expected cs				!
	  sigma_max = 1.0d-43 ! max expected cs					!
	!--------- Exposure range for limits -----------!
		n_ex = 200																			!
		ex_min = 0.001																!
		ex_max = 1.0e13													!
	!--------- Exposure range for limits -----------!
	!	n_Eth = 20																			!
	!	Eth_min = 0.001																!
	!	Eth_max = 50.0																!
	!-----------------------------------------------!
	!-----------------------------------------------!
	allocate(m_vals(nm))
	allocate(ex_vals(n_ex))
	allocate(DL(nm,n_ex))
	m_vals = logspace(m_min,m_max,nm)

	E_th = 10.0
	E_max = 200.0d0
	nucleus = Fluorine
	call GetLimits_MassExposure2(m_min,m_max,nm,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,DL)
	filename = 'nuTest-F.txt'
	open(unit=1000,file=trim(filename))
	write(1000,*) 0.0,logspace(ex_min,ex_max,n_ex)
	do i = 1,nm
		write(1000,*) m_vals(i),DL(i,:)
	end do
	close(1000)


	E_th = 1.8
	E_max = 200.0d0
	nucleus = Helium
	call GetLimits_MassExposure(m_min,m_max,nm,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,DL)
	filename = 'nuTest-He.txt'
	open(unit=1000,file=trim(filename))
	write(1000,*) 0.0,logspace(ex_min,ex_max,n_ex)
	do i = 1,nm
		write(1000,*) m_vals(i),DL(i,:)
	end do
	close(1000)

	call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start
end program runNuFloor