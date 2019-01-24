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
	double precision,dimension(:),allocatable :: m_vals,ex_vals,Eth_vals,DL_ex
	double precision,dimension(:,:),allocatable :: DL,Nsig,Nbg
	double precision :: Eth_max,Eth_min,EB8
	character(len=100) :: filename1,filename2,filename3
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

	!-----Set binning ------------------------------!
	nE_bins = 200 ! Number of energy bins		  		!
	nT_bins = 1	 ! Number of time bins					  !
	nside = 0  ! Order of pixelation (2,4 or 8)		!

	nm = 300
	ns = 300
	n_ex = 2000
	allocate(m_vals(nm))
	allocate(ex_vals(n_ex))
	allocate(DL(nm,ns))
	allocate(DL_ex(nm))
	allocate(Nsig(nm,ns))
	allocate(Nbg(nm,ns))

	nucleus = Fluorine

	E_th = 1.0d-3
	E_max = 200.0d0
	m_min = 0.1
	m_max = 1.0d4
	sigma_min = 1.0d-50
	sigma_max = 1.0d-43
	ex_min = 1.0e-3
	ex_max = 1.0e19
	call GetLimits_MassExposure2(m_min,m_max,nm,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,DL,Nsig,Nbg)
	m_vals = logspace(m_min,m_max,nm)
	filename1 = 'nuTest-F.txt'
	open(unit=1000,file=trim(filename1))
	write(1000,*) 0.0,logspace(sigma_min,sigma_max,ns)
	do i = 1,nm
		write(1000,*) m_vals(i),DL(i,:)
		write(2000,*) Nsig(i,:)
		write(3000,*) Nbg(i,:)
	end do
	close(1000)
	call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start


	nucleus = Helium
	sigma_min = 1.0d-50
	sigma_max = 1.0d-42
	call GetLimits_MassExposure2(m_min,m_max,nm,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,DL,Nsig,Nbg)
	m_vals = logspace(m_min,m_max,nm)
	filename1 = 'nuTest-He.txt'
	open(unit=1000,file=trim(filename1))
	write(1000,*) 0.0,logspace(sigma_min,sigma_max,ns)
	do i = 1,nm
		write(1000,*) m_vals(i),DL(i,:)
		write(2000,*) Nsig(i,:)
		write(3000,*) Nbg(i,:)
	end do
	close(1000)
	call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start



stop

	! High mass Fluorine
	EB8 = MaxNuRecoilEnergy(9)
	E_th = EB8+0.01
	E_max = 200.0d0
	m_min = 1.0
	m_max = 1000.0
	sigma_min = 1.0d-51
	sigma_max = 1.0d-45
	ex_min = 0.1
	ex_max = 1.0e10
	call GetLimits_ExposureGradient(m_min,m_max,nm,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,DL_ex)
	filename1 = 'nuTest-F-High-Ex.txt'
	m_vals = logspace(m_min,m_max,nm)
	open(unit=1000,file=trim(filename1))
	write(1000,*) m_vals
	write(1000,*) DL_ex
	! call GetLimits_MassExposure2(m_min,m_max,nm,ex_min,ex_max,n_ex,sigma_min,sigma_max,ns,Nsig,Nbg)
	! m_vals = logspace(m_min,m_max,nm)
	! filename1 = 'nuTest-F-High.txt'
	! open(unit=1000,file=trim(filename1))
	! write(1000,*) 0.0,logspace(sigma_min,sigma_max,ns)
	! do i = 1,nm
	! 	write(1000,*) m_vals(i),DL(i,:)
	! end do
	! close(1000)
	! call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start


end program runNuFloor
