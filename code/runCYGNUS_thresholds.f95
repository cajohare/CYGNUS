program runCYGNUS_thresholds
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
	double precision :: E_th_F_edisc,E_th_He_edisc,E_th_F_1e,E_th_He_1e,E_th_F_predrift,E_th_He_predrift
	double precision,dimension(:),allocatable :: ex_vals,DL1_Ex,DL2_Ex,DL3_Ex,Nbg1,Nsig1,Nbg2,Nsig2,Nbg3,Nsig3,DLHe_1wimp,DLF_1wimp
	character(len=100) :: filename
	character(len=100) :: fn_end

	! Seed rand and random
	call itime(mytime)
	call cpu_time(clock_start)


! 	write(*,*) '====================================================='
! 	write(*,*) 'Starting program: CYGNUS-1000'
! 	write(*,*) 'Select redout technology (0-10):'
! 	write(*,*) '0 --- Idealised detector (3D)'
! 	write(*,*) '1 --- Pixel grid readout (3D)'
! 	write(*,*) '2 --- Pre drift grid readout (3D)'
! 	write(*,*) '3 --- Post drift grid readout (3D)'
! 	write(*,*) '4 --- Strip readout (3D)'
! 	write(*,*) '5 --- Optical readout (2D)'
! 	write(*,*) '6 --- Wire readout (2D)'
! 	write(*,*) '7 --- Pad readout (1D)'
! 	write(*,*) '8 --- Planar readout (1D)'
! 	write(*,*) '9 --- Nondirectional (0D)'
! 	write(*,*) '10 --- RUN ALL READOUTS'
! 	write(*,*) '11 --- RUN ALL PIXEL READOUTS'

	!------------------Optimisation-----------------!
	!-----Set binning ------------------------------!
	  nE_bins = 20 ! Number of energy bins		  		!
	  nT_bins = 1	 ! Number of time bins					  !
	  nside = 4  ! Order of pixelation (2,4 or 8)		!
	!----- Mass range for limits -------------------!
	  nm = 200 !  Number of mass points							!
	  m_min = 0.1d0 ! Min mass											!
	  m_max = 10000.0d0 ! Max mass										!
	!----- Cross section range for limits ----------!
	  ns = 500 ! resolution of cs scan							!
	  sigma_min = 1.0d-50 ! min. expected cs				!
	  sigma_max = 1.0d-38 ! max expected cs					!

  	n_ex = 20
  	ex_min = 0.1*(0.16/1000.0d0) ! Convert m^3-years into ton-years
  	ex_max = 1.0e10*(0.16/1000.0d0) ! Convert m^3-years into ton-years
	!-----------------------------------------------!
	!-----------------------------------------------!
	allocate(m_vals(nm))
	allocate(DLF(nm))
	allocate(DLHe(nm))
	allocate(DLF_1wimp(nm))
	allocate(DLHe_1wimp(nm))

	allocate(ex_vals(n_ex))
	allocate(DL1_ex(n_ex))
	allocate(DL2_ex(n_ex))
	allocate(DL3_ex(n_ex))
	allocate(Nsig1(n_ex))
	allocate(Nbg1(n_ex))
	allocate(Nsig2(n_ex))
	allocate(Nbg2(n_ex))
	allocate(Nsig3(n_ex))
	allocate(Nbg3(n_ex))

	Time = 6.0
	VolTime = 10000.0*Time
	lat = Boulby(1)
	long = Boulby(2)
	E_max = 200.0d0
	ex_min = ex_min*Time
	ex_max = ex_max*Time

	E_th_F_1e = 0.25
	E_th_He_1e = 0.25

	E_th_F_edisc = 8.0
	E_th_He_edisc = 8.0

	E_th_F_predrift = 0.1
	E_th_He_predrift = 0.1

	m1 = 2.0
	m2 = 10.0
	m3 = 100.0

	! DIRECTIONAL LIMITS


	nside = 4
	readout = 2
	energy_on = .true. ! energy info is currently turned on for best limits
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .true.
	energyres_on = .false.
	searchmode = .false.
	call LoadReadout(readout,	fn_end)


	write(*,*) '----------------------------------------------------'
	Exposure = VolTime*(0.16/1000.0d0)
	filename = '../data/CYGNUS10k-predrift.txt'
	open(unit=123,file=trim(filename))

	nucleus = Helium
	E_th = E_th_He_predrift
	call NwimpEvents(1.0d0,m_min,m_max,nm,m_vals,DLHe_1wimp)
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)

	nucleus = Fluorine
	E_th = E_th_F_predrift
	call NwimpEvents(1.0d0,m_min,m_max,nm,m_vals,DLF_1wimp)
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)

	write(123,*) m_vals
	write(123,*) DLF
	write(123,*) DLHe
	write(123,*) DLF_1wimp
	write(123,*) DLHe_1wimp
	close(123)
	write(*,*) 'writing to: ',trim(filename)
	write(*,*) '----------------------------------------------------'


	write(*,*) '----------------------------------------------------'
	filename = '../data/CYGNUS10k-Exposure-predrift.txt'
	E_th = E_th_He_predrift
	nucleus = Helium
	m = m1
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL1_Ex,Nsig1,Nbg1)

	E_th = E_th_F_predrift
	nucleus = Fluorine
	m = m2
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL2_Ex,Nsig2,Nbg2)

	E_th = E_th_F_predrift
	nucleus = Fluorine
	m = m3
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL3_Ex,Nsig3,Nbg3)


	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(Time*(0.16/1000.0))
	write(1000,*) DL1_Ex
	write(1000,*) Nsig1
	write(1000,*) Nbg1
	write(1000,*) DL2_Ex
	write(1000,*) Nsig2
	write(1000,*) Nbg2
	write(1000,*) DL3_Ex
	write(1000,*) Nsig3
	write(1000,*) Nbg3
	close(1000)
	write(*,*) '----------------------------------------------------'






	! NON-DIRECTIONAL
	nside = 0
	readout = 2
	energy_on = .true. ! energy info is currently turned on for best limits
	angres_on = .true.
	efficiency_on = .true.
	headtail_on = .true.
	energyres_on = .false.
	searchmode = .false.
	call LoadReadout(readout,	fn_end)


	write(*,*) '----------------------------------------------------'
	Exposure = VolTime*(0.16/1000.0d0)
	filename = '../data/CYGNUS10k-1e.txt'
	open(unit=123,file=trim(filename))

	nucleus = Helium
	E_th = E_th_He_1e
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)
	call NwimpEvents(1.0d0,m_min,m_max,nm,m_vals,DLHe_1wimp)

	nucleus = Fluorine
	E_th = E_th_F_1e
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)
	call NwimpEvents(1.0d0,m_min,m_max,nm,m_vals,DLF_1wimp)

	write(123,*) m_vals
	write(123,*) DLF
	write(123,*) DLHe
	write(123,*) DLF_1wimp
	write(123,*) DLHe_1wimp
	close(123)
	write(*,*) 'writing to: ',trim(filename)
	write(*,*) '----------------------------------------------------'


	write(*,*) '----------------------------------------------------'
	Exposure = VolTime*(0.16/1000.0d0)
	filename = '../data/CYGNUS10k-edisc.txt'
	open(unit=123,file=trim(filename))

	nucleus = Helium
	E_th = E_th_He_edisc
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLHe)
  call NwimpEvents(1.0d0,m_min,m_max,nm,m_vals,DLHe_1wimp)

	nucleus = Fluorine
	E_th = E_th_F_edisc
	call GetLimits(m_min,m_max,nm,sigma_min,sigma_max,ns,	m_vals,DLF)
	call NwimpEvents(1.0d0,m_min,m_max,nm,m_vals,DLF_1wimp)

	write(123,*) m_vals
	write(123,*) DLF
	write(123,*) DLHe
	write(123,*) DLF_1wimp
	write(123,*) DLHe_1wimp
	close(123)
	write(*,*) 'writing to: ',trim(filename)
	write(*,*) '----------------------------------------------------'




	write(*,*) '----------------------------------------------------'
	filename = '../data/CYGNUS10k-Exposure-1e.txt'
	E_th = E_th_He_1e
	nucleus = Helium
	m = m1
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL1_Ex,Nsig1,Nbg1)

	E_th = E_th_F_1e
	nucleus = Fluorine
	m = m2
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL2_Ex,Nsig2,Nbg2)

	E_th = E_th_F_1e
	nucleus = Fluorine
	m = m3
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL3_Ex,Nsig3,Nbg3)

	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(Time*(0.16/1000.0))
	write(1000,*) DL1_Ex
	write(1000,*) Nsig1
	write(1000,*) Nbg1
	write(1000,*) DL2_Ex
	write(1000,*) Nsig2
	write(1000,*) Nbg2
	write(1000,*) DL3_Ex
	write(1000,*) Nsig3
	write(1000,*) Nbg3
	close(1000)
	write(*,*) '----------------------------------------------------'






	write(*,*) '----------------------------------------------------'
	filename = '../data/CYGNUS10k-Exposure-edisc.txt'
	E_th = E_th_He_edisc
	nucleus = Helium
	m = m1
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL1_Ex,Nsig1,Nbg1)

	E_th = E_th_F_edisc
	nucleus = Fluorine
	m = m2
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL2_Ex,Nsig2,Nbg2)

	E_th = E_th_F_edisc
	nucleus = Fluorine
	m = m3
	call GetLimits_Exposure(m,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,ex_vals,DL3_Ex,Nsig3,Nbg3)

	open(unit=1000,file=trim(filename))
	write(1000,*) ex_vals/(Time*(0.16/1000.0))
	write(1000,*) DL1_Ex
	write(1000,*) Nsig1
	write(1000,*) Nbg1
	write(1000,*) DL2_Ex
	write(1000,*) Nsig2
	write(1000,*) Nbg2
	write(1000,*) DL3_Ex
	write(1000,*) Nsig3
	write(1000,*) Nbg3
	close(1000)
	write(*,*) '----------------------------------------------------'





end program runCYGNUS_thresholds
