program runCYGNUS_Final2
	use params
	use WIMPFuncs
	use util
	use LabFuncs
	use like
	use NeutrinoFuncs
	implicit none

	integer :: jj,ns,nm,readout_selection,loopmin,loopmax,n_ex
	double precision :: m_min,m_max,sigma_min,sigma_max,ex_min,ex_max,Vol,Time,m,m1,m2,m3
	double precision,dimension(:),allocatable :: m_vals,DLHe,DLF
	double precision :: E_th_F_edisc,E_th_He_edisc,E_th_F,E_th_He
	double precision,dimension(:),allocatable :: ex_vals,DL1_Ex,DL2_Ex,DL3_Ex,Nbg1,Nsig1,Nbg2,Nsig2,Nbg3,Nsig3,DLHe_1wimp,DLF_1wimp
	logical :: saveNwimp
	character(len=100) :: filename
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)

	!------------------Optimisation-----------------!
	!-----Set binning ------------------------------!
	  nE_bins = 20 ! Number of energy bins
	  nT_bins = 1	 ! Number of time bins
	!----- Mass range for limits -------------------!
	  nm = 50 !  Number of mass points
	  m_min = 0.1d0 ! Min mass
	  m_max = 1000.0d0 ! Max mass
	!----- Cross section range for limits ----------!
	  ns = 100 ! resolution of cs scan
	  sigma_min = 1.0d-50 ! min. expected cs
	  sigma_max = 1.0d-36 ! max expected cs


	allocate(m_vals(nm))
	allocate(DLF(nm))
	allocate(DLHe(nm))

	Time = 6.0
	VolTime = 1000.0*Time
	lat = Boulby(1)
	long = Boulby(2)

	energy_on = .true. ! energy info is currently turned on for best limits
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .true.
	energyres_on = .true.
	searchmode = .false.
  readout = 1 ! Predrift

  !Vol = 1000.0
  !filename = '../data/CYGNUS2000-Final.txt'

  !Vol = 10000.0
  !filename = '../data/CYGNUS20k-Final.txt'

  Vol = 100000.0
  filename = '../data/CYGNUS200k-Final.txt'


	open(unit=123,file=trim(filename))
	m_vals = logspace(m_min,m_max,nm)
	write(123,*) m_vals


  E_th_F = 0.25
  E_th_He = 0.25
  !=======!=======!=======!=======!=======!=======!=======!=======!=======!=======
  ! Directional mode
  nside = 4
  call LoadReadout(readout,	fn_end)

  Exposure = Vol*6*(0.16/1000.0d0)*(755.0/740.0)
	nucleus = Helium
	E_th = E_th_He
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)
	Exposure = Vol*6*(0.16/1000.0d0)*(20.0/20.0)
	nucleus = Fluorine
	E_th = E_th_F
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)
  write(123,*) DLHe
  write(123,*) DLF
  !=======!=======!=======!=======!=======!=======!=======!=======!=======!=======

  !=======!=======!=======!=======!=======!=======!=======!=======!=======!=======
  ! Search mode
  nside = 0
  call LoadReadout(readout,	fn_end)
  Exposure = Vol*6*(0.16/1000.0d0)*(160.0/740.0)
	nucleus = Helium
	E_th = E_th_He
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)
	Exposure = Vol*6*(0.16/1000.0d0)*(600.0/20.0)
	nucleus = Fluorine
	E_th = E_th_F
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)
  write(123,*) DLHe
  write(123,*) DLF
  !=======!=======!=======!=======!=======!=======!=======!=======!=======!=======








  E_th_F = 8.0
  E_th_He = 8.0
  !=======!=======!=======!=======!=======!=======!=======!=======!=======!=======
  ! Directional mode
  nside = 4
  call LoadReadout(readout,	fn_end)

  Exposure = Vol*6*(0.16/1000.0d0)*(755.0/740.0)
	nucleus = Helium
	E_th = E_th_He
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)
	Exposure = Vol*6*(0.16/1000.0d0)*(20.0/20.0)
	nucleus = Fluorine
	E_th = E_th_F
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)
  write(123,*) DLHe
  write(123,*) DLF
  !=======!=======!=======!=======!=======!=======!=======!=======!=======!=======

  !=======!=======!=======!=======!=======!=======!=======!=======!=======!=======
  ! Search mode
  nside = 0
  call LoadReadout(readout,	fn_end)
  Exposure = Vol*6*(0.16/1000.0d0)*(160.0/740.0)
	nucleus = Helium
	E_th = E_th_He
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)
	Exposure = Vol*6*(0.16/1000.0d0)*(600.0/20.0)
	nucleus = Fluorine
	E_th = E_th_F
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)
  write(123,*) DLHe
  write(123,*) DLF
  !=======!=======!=======!=======!=======!=======!=======!=======!=======!=======



end program runCYGNUS_Final2
