program runCYGNUS1000
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax
	double precision :: m_min,m_max,sigma_min,sigma_max
	character(len=100) :: filename
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)
	
	
	!------------------Optimisation-----------------!
	!-----Set binning ------------------------------!	
	  nE_bins = 20 ! Number of energy bins			!
	  nT_bins = 1 ! Number of time bins				!
	  nside = 4  ! Order of pixelation (2,4 or 8)	!
	!----- Mass range for limits -------------------!	
	  nm = 100 !  Number of mass points				!
	  m_min = 0.5d0 ! Min mass						!
	  m_max = 1000.0d0 ! Max mass					!
	!----- Cross section range for limits ----------!	
	  ns = 500 ! resolution of cs scan				!
	  sigma_min = 1.0d-49 ! min. expected cs		!
	  sigma_max = 1.0d-40 ! max expected cs			!
	!-----------------------------------------------!
	!-----------------------------------------------!
	
	
	
	! CYGNUS 1000 x 3 years of running (can be changed for whatever preference)
	VolTime = 1000.0d2*3.0d0
	energy_on = .true. ! energy info is currently turned on for best limits

	write(*,*) '====================================================='
	write(*,*) 'Starting program: CYGNUS-1000'
	write(*,*) 'Select redout technology (0-8):'
	write(*,*) '0 --- Idealised detector (3D)'
	write(*,*) '1 --- Pixel grid readout (3D)'
	write(*,*) '2 --- Strip readout (3D)'
	write(*,*) '3 --- Optical readout (2D)'
	write(*,*) '4 --- Wire readout (2D)'
	write(*,*) '5 --- Pad readout (1D)'
	write(*,*) '6 --- Planar readout (1D)'
	write(*,*) '7 --- Nondirectional (0D)'
	write(*,*) '8 --- RUN ALL READOUTS'
	readout_selection = -1
	do while ((readout_selection.gt.8).or.(readout_selection.lt.0)) 
		write(*,*) 'Select 1-8'
		write(*,*) '>...'
		read(*,*) readout_selection
	end do

	if (readout_selection.eq.8) then
		loopmax = 7
		loopmin = 0
		write(*,*) 'Running all readouts'
	else
		loopmin = readout_selection
		loopmax = readout_selection
	end if

	write(*,*) '----------------------------------------------------'
	do readout = loopmin,loopmax
		call LoadReadout(readout,	fn_end)
		write(*,*) 'Readout = ',fn_end

		if ((readout.eq.0).or.((readout.eq.7))) then
			angres_on = .false.
			efficiency_on = .false.
			headtail_on = .false.
			energyres_on = .false.
		else
			angres_on = .true.
			efficiency_on = .true.
			headtail_on = .true.
			energyres_on = .false.
		end if

		if (readout.eq.7) then
			nside = 0
		end if

		filename = '../data/CYGNUS1000-'//trim(fn_end)//'.txt'
		call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
		call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start
	end do

  end program runCYGNUS1000
