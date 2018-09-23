program runLimits
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: ns,nm
	double precision :: m_min,m_max,sigma_min,sigma_max
	character(len=100) :: filename

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)

	!----------------------------------------------------------!
	! Set binning (CONTROLS ACCURACY/EFFICIENCY OF CODE)
	  nE_bins = 20
	  nT_bins = 1
	  nside = 4  
	! Mass range for limits
	  nm = 100
	  m_min = 0.5d0
	  m_max = 1000.0d0
	! Cross section range for limits
	  ns = 500
	  sigma_min = 1.0d-49
	  sigma_max = 1.0d-40
	!----------------------------------------------------------!


	! CYGNUS 1000
	  ! Idealised
	  VolTime = 1000.0*3.0d0
	  sig_E = 0.0d0
	  energy_on = 1
	  angres_on = 0
	  eff_on = 0
	  headtail_on = 0
	  readout = 1
	  filename = '../data/CYGNUS1000-ideal.txt'
	  call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)

	  ! With angular resolution
	  angres_on = 1
	  eff_on = 0
	  headtail_on = 0
	  readout = 1
	  filename = '../data/CYGNUS100k-angres.txt'
	  call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)



	call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start
end program runLimits
