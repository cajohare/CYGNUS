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
	double precision,dimension(:),allocatable :: ex_vals,DL_F,DL_He,Nbg_F,Nsig_F,Nbg_He,Nsig_He
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
	  m_max = 1000.0d0 ! Max mass										!
	!----- Cross section range for limits ----------!
	  ns = 1000 ! resolution of cs scan							!
	  sigma_min = 1.0d-51 ! min. expected cs				!
	  sigma_max = 1.0d-38 ! max expected cs					!
	!-----------------------------------------------!
	!-----------------------------------------------!

	n_ex = 20
	ex_min = 0.001*6.0*(0.16/1000.0d0) ! Convert m^3-years into ton-years
	ex_max = 1.0e10*6.0*(0.16/1000.0d0) ! Convert m^3-years into ton-years
	allocate(ex_vals(n_ex))
	allocate(DL_He(n_ex))
	allocate(DL_F(n_ex))
	allocate(Nsig_F(n_ex))
	allocate(Nbg_F(n_ex))
	allocate(Nsig_He(n_ex))
	allocate(Nbg_He(n_ex))
	lat = Boulby(1)
	long = Boulby(2)
	E_max = 200.0d0

	! CYGNUS 1000 x 3 years of running (can be changed for whatever preference)
	Vol = 1000.0
	Time = 6.0d0
	VolTime = Vol*Time
	energy_on = .true. ! energy info is currently turned on for best limits
	searchmode = .false.

	write(*,*) '====================================================='
	write(*,*) 'Starting program: CYGNUS-1000'
	write(*,*) 'Select redout technology (0-10):'
	write(*,*) '0 --- Idealised detector (3D)'
	write(*,*) '1 --- Pixel grid readout (3D)'
	write(*,*) '2 --- Pre drift grid readout (3D)'
	write(*,*) '3 --- Post drift grid readout (3D)'
	write(*,*) '4 --- Strip readout (3D)'
	write(*,*) '5 --- Optical readout (2D)'
	write(*,*) '6 --- Wire readout (2D)'
	write(*,*) '7 --- Pad readout (1D)'
	write(*,*) '8 --- Planar readout (1D)'
	write(*,*) '9 --- Nondirectional (0D)'
	write(*,*) '10 --- RUN ALL READOUTS'
	write(*,*) '11 --- RUN ALL PIXEL READOUTS'
	readout_selection = -1
	do while ((readout_selection.gt.11).or.(readout_selection.lt.0))
		write(*,*) 'Select 0-11'
		write(*,*) '>...'
		read(*,*) readout_selection
	end do

	if (readout_selection.eq.10) then
		loopmax = 9
		loopmin = 0
		write(*,*) 'Running all readouts'
	elseif (readout_selection.eq.11) then
		loopmax = 2
		loopmin = 0
	else
		loopmin = readout_selection
		loopmax = readout_selection
	end if

	write(*,*) '----------------------------------------------------'
	do readout = loopmin,loopmax
		call LoadReadout(readout,	fn_end)
		write(*,*) 'Readout = ',fn_end

		if (readout.eq.0) then
			angres_on = .false.
			efficiency_on = .false.
			headtail_on = .false.
			energyres_on = .false.
		elseif (readout.eq.9) then
			angres_on = .false.
			efficiency_on = .true.
			headtail_on = .false.
			energyres_on = .false.
		else
			angres_on = .true.
			efficiency_on = .true.
			headtail_on = .true.
			energyres_on = .false.
		end if

		if (readout.eq.9) then
			nside = 0
		end if


		filename = '../data/CYGNUS-Exposure-'//trim(fn_end)//'.txt'


		E_th = 3.0d0
		nucleus = Fluorine
		m = 9.0
		call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_F,Nsig_F,Nbg_F)

		E_th = 1.8d0
		nucleus = Helium
		m = 1000.0
		call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL_He,Nsig_He,Nbg_He)

		open(unit=1000,file=trim(filename))
		write(1000,*) ex_vals/(6.0*(0.16/1000.0))
		write(1000,*) DL_F
		write(1000,*) Nsig_F
		write(1000,*) Nbg_F
		write(1000,*) DL_He
		write(1000,*) Nsig_He
		write(1000,*) Nbg_He
		close(1000)

		call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start
	end do

end program runCYGNUS_Exposure
