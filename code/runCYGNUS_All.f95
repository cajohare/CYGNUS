program runCYGNUS_All
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
	  nE_bins = 50 ! Number of energy bins		  		
	  nT_bins = 1	 ! Number of time bins					  
	  nside = 4  ! Order of pixelation (2,4 or 8)		
	!----- Mass range for limits -------------------!
	  nm = 50 !  Number of mass points							
	  m_min = 0.5d0 ! Min mass											
	  m_max = 1000.0d0 ! Max mass										
	!----- Cross section range for limits ----------!
	  ns = 100 ! resolution of cs scan							
	  sigma_min = 1.0d-49 ! min. expected cs				
	  sigma_max = 1.0d-40 ! max expected cs					
	!-----------------------------------------------!
	!-----------------------------------------------!



	! CYGNUS 1000 x 3 years of running (can be changed for whatever preference)
	E_th_F = 2.5
	E_th_He = 1.8
	Vol = 1000.0
	Time = 6.0d0
	VolTime = Vol*Time
	energy_on = .true. ! energy info is currently turned on for best limits
	searchmode = .false.
	
	write(*,*) '====================================================='
	write(*,*) 'Select run (0-10):'
	write(*,*) '0 --- Pixel Ideal'
	write(*,*) '1 --- Pixel pre-drift'
	write(*,*) '2 --- Pixel post-drift'
	write(*,*) '3 --- Gem readout'
	write(*,*) '4 --- Wire readout'
	write(*,*) '5 --- Pad readout'
	write(*,*) '6 --- Strip readout'
	write(*,*) '7 --- Pixel readout'
	write(*,*) '8 --- RUN ALL READOUTS'
	
	readout_selection = -1
	do while ((readout_selection.gt.8).or.(readout_selection.lt.0))
		write(*,*) 'Select 0-8'
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
		efficiency_on = .true.
		energyres_on = .true.
		
		if (readout.eq.0) then ! IDEAL
			angres_on = .false.
			headtail_on = .false.
		else
			angres_on = .true.
			headtail_on = .true.
		end if


		if (Vol.eq.1000.0) then
			filename = '../data/CYGNUS1000-'//trim(fn_end)//'.txt'
		elseif (Vol.eq.100000.0) then
			filename = '../data/CYGNUS100k-'//trim(fn_end)//'.txt'
		elseif (Vol.eq.10000.0) then
			filename = '../data/CYGNUS10k-'//trim(fn_end)//'.txt'
		end if

		call CYGNUSLimit(E_th_F,E_th_He,m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
		call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start
	end do

end program runCYGNUS_All
