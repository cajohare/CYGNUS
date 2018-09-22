program runLimits
  use params
  use WIMPFuncs
  use util
  use LabFuncs
  use like
  use NeutrinoFuncs
  implicit none

  integer :: ns,nm,i,n_ex,pix
  double precision :: m_min,m_max,sigma_min,sigma_max,exposure_min,exposure_max,Ec,x_pix_db(3)
  character(len=100) :: filename

  ! Seed rand and random
  call itime(mytime)
  call cpu_time(clock_start)
   
  !---------------- CODE SPEED/EFFICIENCY/ACCURACY------------!
  ! SET BINNING (CONTROLS EFFICIENCY OF CODE)
  nE_bins = 20
  nT_bins = 1 ! Min 1
  nside = 4  
  nm = 100 ! mass limits
  m_min = 0.5d0
  m_max = 1000.0d0
  ns = 500 ! cross section limits
  sigma_min = 1.0d-49
  sigma_max = 1.0d-40
  !----------------------------------------------------------!
  
  VolTime = 100.0d3*3.0d0
  sig_E = 0.0d0
  energy_on = 1
  angres_on = 0
  eff_on = 0
  headtail_on = 0
  readout = 1
  filename = 'data/CYGNUS100k-ideal.txt'
  call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
  
  
  angres_on = 1
  eff_on = 0
  headtail_on = 0
  readout = 1
  filename = 'data/CYGNUS100k-angres.txt'
  call CYGNUSLimit(m_min,m_max,nm,sigma_min,sigma_max,ns,filename)
  
  
  
  call cpu_time(clock_stop); write(*,*) 'Time elapsed = ',clock_stop-clock_start
end program runLimits
