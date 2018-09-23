program runMollweide
  use params
  use NeutrinoFuncs
  use WIMPFuncs
  use util
  use LabFuncs
  implicit none

  double precision,dimension(:,:),allocatable :: skymap_wimp,skymap_nu
  double precision,dimension(:),allocatable :: xvals,yvals,Evals,dRdE_nu,dRdE_wimp
  double precision :: th,ph,q_rec(3),JD,x,y
  integer :: expt,np,ti,i,j,k,si,choice
  character(len=100) :: filename_wimp,filename_nu
  call itime(mytime)
  call cpu_time(clock_start)

  
  ! EXPT FIXED PARAMS
  E_th = 3.0d0
  E_max = 100.0d0
  sig_E = 0.0d0
  angres_on = 0
  eff_on = 0
  headtail_on = 0
  Exposure = 0.16d0*3
  
  !---------------- SETUP MOLLWEIDE-----------------------!
  np = 500
  nE_bins = 50
  nT_bins = 1
  nside = 0
  allocate(Evals(nE_bins))
  allocate(dRdE_nu(nE_bins))
  allocate(dRdE_wimp(nE_bins))
  allocate(xvals(2*np))
  allocate(yvals(np))
  allocate(skymap_nu(np,2*np))
  allocate(skymap_wimp(np,2*np))
  xvals = linspace(-2*sqrt(2.0d0),2*sqrt(2.0d0),2*np)
  yvals = linspace(-sqrt(2.0d0),sqrt(2.0d0),np)
  !----------------------------------------------------------!
  
  
  
  ! DO MOLLWEIDE
  do choice = 1,2
	  if (choice.eq.1) then
		  nucleus = Helium
		  E_th = 1.2d0
		  filename_wimp = '../data/Moll-He-wimp.txt'
		  filename_nu = '../data/Moll-He-nu.txt'
		  Evals = logspace(E_th,E_max,nE_bins)  
	   elseif (choice.eq.2) then
		   nucleus = Fluorine
		   E_th = 3.0d0
		   filename_wimp = '../data/Moll-F-wimp.txt'
		   filename_nu = '../data/Moll-F-nu.txt'
		   Evals = logspace(E_th,E_max,nE_bins)
	   end if
	   
	   ! LOAD DATA
	   call GetNuFluxes ! Backgrounds	
	   call PreAllocate ! Data		
	   call SHM
	   JD = JulianDay(9,6,2018,18.0d0)-Jan1 ! September 6th
	   v_lab_all(1,:) = LabVelocity(JD)
	   T_bin_centers(1) = JD 
	   ! FIX WIMP
	   m_chi = 6.0d0
	   sigma_p = 1.0d-45
		  
	  skymap_wimp = 0.0d0
	  skymap_nu = 0.0d0
	  do i = 1,np
	  	do j = 1,2*np
			y = yvals(i)
			x = xvals(j)
		 	if (((x/(2*sqrt(2.0d0)))**2.0d0 + (y/sqrt(2.0d0))**2.0d0).lt.1.0d0) then
		    	th = asin(y/sqrt(2.0d0))
		    	ph =  pi - pi*x/(2*sqrt(2.0d0)*sqrt(1.0d0-y**2.0d0/2.0d0))
		    	q_rec = -1.0*(/cos(th)*cos(ph),cos(th)*sin(ph),sin(th)/)
				dRdE_nu = 0.0d0
		    	do k = 1,nE_bins
					do si = 1,n_bg-2
						dRdE_nu(k) = dRdE_nu(k)+NeutrinoRecoilSpectrum_Solar(Evals(k)*q_rec,1,E_nu_all(:,si),Flux_all(:,si))
					end do
					dRdE_wimp(k) = WIMPRate_Direction(Evals(k)*q_rec,1)
		    	end do
		    	skymap_nu(i,j) = sum(0.5*(Evals(2:nE_bins)-Evals(1:nE_bins-1))*(dRdE_nu(2:nE_bins)+dRdE_nu(1:nE_bins-1)))
		    	skymap_wimp(i,j) = sum(0.5*(Evals(2:nE_bins)-Evals(1:nE_bins-1))*(dRdE_wimp(2:nE_bins)+dRdE_wimp(1:nE_bins-1)))
		  	end if
	    end do
	    write(*,*) i,'/',np
	  end do

	  ! SAVE
	  write(*,*) 'Writing to: ',trim(filename_wimp), '...'
	  write(*,*) 'Writing to: ',trim(filename_nu), '...'
	  open(unit=111,file=trim(filename_wimp))
	  open(unit=222,file=trim(filename_nu))
	  do i = 1,np
	    write(111,*) skymap_wimp(i,:)
	    write(222,*) skymap_nu(i,:)
	  end do
	  close(111)
	  close(222)
	  
	  call Unallocate
	end do

  call cpu_time(clock_stop)
  write(*,*) 'Time elapsed = ',clock_stop-clock_start
end program runMollweide
