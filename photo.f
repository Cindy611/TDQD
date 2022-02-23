      PROGRAM PHOTO                                                     
C                                                                       
C.....Modified 11 Sept 1997 .. Gabriel Balint-Kurti
C.....This program computes the total energy absorption cross section   
C.....for the photodissociation of a diatomic molecule using            
C.....Time-Dependent Quantum Dynamics [ E.J. Heller, J. Chem. Phys. 68, 
C.....2066 (1978).]. The excitation is assumed to be from a ground     
C.....state to a single repulsive excited state. A wavepacket is formed 
C.....by the multiplying an initial vibrational wavefunction by a       
C.....transition dipole functionmoment. This is then evolved in time    
C.....using a Chebychev polynomial expansion of the evolution operator  
C.....[R. Kosloff, J. Phys. Chem., 92, 2087(1988)]. The program         
C.....calculates the total absorption cross section as a function of    
C.....photon energy for the transition by Fourier Transform analysis of 
C.....the Autocorrelation Function [ G.G. Balint-Kurti, R.N. Dixon and  
C.....C. Clay Marston, J. Chem. Soc. Faraday Trans., 86, 1741(1990)].   
C                                                                       
C.....An initial vibrational wavefunction (plus the corresponding       
C.....eigenvalue) is computed from a bound potential energy curve which 
C.....is supplied by the user. The excited (repulsive) state potential  
C.....energy curve and a transition dipole moment function must also be 
C.....supplied.                                                         
C.....The ground state wavefunction is calculated using the Fourier     
C.....Grid Hamiltonian method (C. C. Marston and G. G. Balint-Kurti,    
C.....J. Chem. Phys., 91, 3571 (1989); G. G. Balint-Kurti, C. L.        
C.....Ward and C. C. Marston, Comput. Phys. Commun., 67, 285 (1991)).   
C.....This requires input of the ground state potential governing the   
C.....nuclear motion plus the reduced mass of the system. A specific    
C.....vibrational state wavefunction and corresponding eigenvalue       
C.....is selected from the output. For the propagation the excited             
C.....state potential is required along with the dipole moment for         
C.....the spectroscopic transition to form the initial wavepacket. 
C                                                                       
C.....The Fourier Transform methods used require the absorption         
C.....(damping) of the wavefunction at the end of the grid. This is     
C.....achieved using a complex damping potential [ A. Vibok G.G.        
C.....Balint-Kurti 96, 7615(1992); A. Vibok and G. G. Balint-Kurti J.   
C.....Phys. Chem. (submitted)]. The current program is for the          
C.....photodissociation of HCl or DCl, with an additional facility to   
C.....choose the initial vibrational state for either molecule.         
C                                                                       
C********************************************************************
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
C.....NX    : number of grid points.                                       
C.....NRECT : number of time steps.                                     
C.....Active parameters below for HCl.                                 
      PARAMETER (NX=256,NRECT=100)  
      PARAMETER (NSIZE=4096)                                            
C.....The parameters for DCl are:                                      
C      PARAMETER (NX=256,NRECT=80)                                      
C                                                                       
C.....Declare arrays.                                                   
C.....C00(NX) contains initial wavefunction.                            
C.....C(NX)   contains wavefunction at current time t.                    
      COMPLEX*16 C00(NX),C(NX)                                          
      COMPLEX*16 PHIJ(NX)                                               
      COMPLEX*16 X1A(NX),CR(NX)                                         
      COMPLEX*16 CCD(NSIZE)        
      COMPLEX*16 CZERO                                                  
C                                                                       
C.....Arrays for Fast Fourier Transform routines.                       
      COMPLEX*16 FF1(NX)
      COMPLEX*16 CD
C                                                                       
C.....VVN(NX) contains potential evaluated on the grid.                 
      DIMENSION VVN(NX)                                                 
      DIMENSION TIMEFS(NRECT+1),AUTOC(NRECT+1)                          
      DIMENSION TVN(NX)
C                                                                       
C.....Temporary storage.                                                
      DIMENSION SCRV1(NX),SCRV2(NX),SCRV3(NX),SCRV4(NX)                 
      DIMENSION SCRM1(NX,NX),SCRM2(NX,NX)                               
      DIMENSION XPL(NX)                                                 
C                                                                       
C.....Variable data input.                                              
C.....ZMA : atomic mass of atom A.                                      
C.....ZMB : atomic mass of atom B.                                      
C.....RO  : equilibrium bond separation         
C
C.....IVS     : initial vibrational ground state from which                 
C.....          photodissociation occurs.                                     
C.....RLENGTH : grid length (expressed as a multiple of the equilibrium 
C.....          separation).                                                  
C.....RMIN    : starting point of grid (expressed as a multiple of the     
C.....          equilibrium separation).                                      
C.....RDAMP   : starting point of damping region (expressed as a          
C.....          multiple of the equilibrium separation).                      
C.....TTIME   : total time for the propagation cycle.                     
C.....Active parameters below for HCl.                                 
      DATA ZMA,ZMB/1836.9822D0,64621.91975D0/                           
C.... The parameters for DCl are:                                       
C     DATA ZMA,ZMB/3670.5381D0,64621.91975D0/                           
      DATA RO/2.408558D0/                                               
      DATA IVS/0/                                                       
      DATA RLENGTH/5.0D0/                                               
      DATA RMIN/0.21D0/                                                 
      DATA RDAMP/4.0D0/                                                 
C.....Active parameters below for HCl.                                 
      DATA TTIME/500.0D0/                                               
C.... The parameters for DCl are:                                       
C      DATA TTIME/640.0D0/                                              
C                                                                       
C.....IWRIT : parameter for varying screen output.                      
C..... 0 for restricted output                                         
C..... 1 for more detailed output                                      
      DATA IWRIT/1/                                                     
C.....ICHEX : parameter for calculation of energy at each time step.    
C..... 0 for no calculation                                            
C..... 1 for calculation of Complex Hamiltonian Expectation            
C.....   value (CHEX). Data written to file.                             
      DATA ICHEX/1/                                                     
C.....INORM : parameter for evaluation of wavepacket norm at            
C.....        each time step.                                                   
C..... 0 for no calculation                                            
C..... 1 for calculation of <C|C>. Data written to file.               
      DATA INORM/1/                                                     
C.....IAUTOC : parameter for evaluation of autocorrelation function     
C.....         at end of propagation.                                            
C..... 0 for no calculation                                            
C..... 1 for calculation. Data written to file.                        
      DATA IAUTOC/1/                                                    
C.....IUNITA : file destination of absorption cross section data.       
C.....IUNITC : file destination of (any) CHEX data.                     
C.....IUNITF : file destination of (any) autocorrelation data.          
C.....IUNITN : file destination of (any) wavepacket norm data.          
      DATA IUNITA,IUNITC,IUNITF,IUNITN/7,8,9,10/                        
C                                                                       
C.....Constants.                                                        
      DATA CZERO/(0.0D0,0.0D0)/                                         
      DATA PI/3.141592653589793D0/                                      
C.....Conversion factors.                                               
C.....HTEV : HARTREES  --> eV                                            
C.....ATB  : ANGSTROMS --> BOHR                                          
      DATA HTEV/27.211648D0/                                            
C      DATA ATB/1.88973D0/                 
C                                                                       
C.....Grid dimensions.                                                  
C.....Calculate length of grid (ZLZ).  
      OPEN(6, FILE='photo.out', ACCESS='SEQUENTIAL')
      ZLZ=RLENGTH*RO                                                    
      WRITE(6,*) ' '                                                   
      WRITE(6,*) 'All quantities used in the program are in au
     & (unless otherwise stated).'      
      WRITE(6,*) 'The energy absorption cross sections are output in
     & Angstroms squared'
      WRITE(6,*) 'and the photon energy is given in cm-1'      
      WRITE(6,*) ' '                                                   
      WRITE(6,'('' Number of grid points in 1D grid (NX) = '',I3)') NX 
      WRITE(6,'('' Range of grid (ZLZ) = '',F11.8)') ZLZ               
C.....Specify start of grid.                                            
      X0=RMIN*RO                                            
      WRITE(6,'('' Starting point for grid (X0) = '',F11.8)') X0       
      WRITE(6,'('' End point of grid = '',F11.8)') (X0+ZLZ)            
C.....Calculate spatial increment DZ.                                   
      DZ=ZLZ/DFLOAT(NX)                                                 
      WRITE(6,'('' Spatial increment (DZ) =  '',F11.8)') DZ            
C.....Fill out XPL(NX) array with spatial coordinates on grid.          
      ZZ=X0                                                             
       DO 002 I=1,NX                                                    
          XPL(I)=ZZ                                                     
          ZZ=ZZ+DZ                                                      
002    CONTINUE                                                         
C                                                                       
C.....Calculate ground state vibrational wavefunction.                  
      CALL GROUND(NX,ZMA,ZMB,ZMU,ZLZ,IVS,SCRV1,SCRM1,SCRV2,             
     &            SCRV3,SCRV4,SCRM2,VMAX,RO,X0,DZ,CZERO,C00,IWRIT)     
C.....Save vibrational state Eigenvalue.                                
      EIG=SCRV2(IVS+1)                                                  
       IF(IWRIT.GT.0) THEN                                              
          WRITE(6,'('' Eigenvalue of ground vibrational state (EIG) =
     & '',F11.8)') EIG                                      
       END IF                                                        
C                                                                       
C.....Define damping region.                                            
      NDAMP=INT((RDAMP*RO-X0)*DFLOAT(NX)/ZLZ)                           
      WRITE(6,'('' Damping to start at grid point (NDAMP) = '',I3)')
     & NDAMP                                                      
      WRITE(6,'('' Damping region starts at '',F11.8)') (NDAMP*DZ+X0)  
C.....Calculate minimum momentum on grid.                               
      ZKMIN=-PI/DZ                                                      
       IF(IWRIT.EQ.1) THEN                                              
          WRITE(6,'('' Minimum momentum on grid (ZKMIN) = '', 
     & F13.8)') ZKMIN                                
       END IF                                                           
C.....Calculate increment in momentum space (DKZ).                      
      DKZ=2.0D0*PI/ZLZ                                                  
       IF(IWRIT.EQ.1) THEN                                              
          WRITE(6,'('' Increment in momentum space (DKZ) = '',
     & F12.8)') DZ                                             
       END IF                                                           
C                                                                       
C.....Set maximum kinetic energy.                                       
      TMAX=(ZKMIN**2)/(2.0D0*ZMU)                                       
       IF(IWRIT.EQ.1) THEN                                              
          WRITE(6,'('' Maximum kinetic energy on grid (TMAX) = '',    
     & F12.8)') TMAX                                
       END IF                                                           
C.....Set VMAX to 10eV.                                                 
      VMAX=10.D0/HTEV                                                   
C.....Set VMIN to zero.                                                 
      VMIN=0.0D0                                                        
      DELTAV=VMAX-VMIN                                                  
      DELTAE=TMAX+DELTAV                                                
C.....E0  is the zero of energy for energy scaling (i.e. the centre      
C.....    of the range).                                                
C.....ESC is the energy scaling parameter ( so that scaled energy       
C.....    lies between -1 and 1).                                       
      ESC=0.5D0*DELTAE                                                  
      E0= VMIN + ESC                                                    
C.....Calculate and store kinetic energy array.                         
       DO 003 I=1,NX                                                   
          ZKA=ZKMIN+DFLOAT(I-1)*DKZ                                     
          SCRV1(I)=(ZKA**2)/(2.0D0*ZMU)                                 
          SCRV1(I)=SCRV1(I)/ESC                                         
003    CONTINUE                                                        
C                                                                       
C.....Kinetic energy array is now in SCRV1(I). Re-order for FFT         
C.....and place in TVN(I).                                              
      CALL UTGR(SCRV1,TVN,NX)                                           
C                                                                       
C.....Multiply initial wavefunction by dipole moment operator           
C.....to obtain wavepacket on the repulsive state.                      
C      OPEN(30, FILE='initial_wavepacket.txt')
       DO 004 I=1,NX 
          CALL DIPMOM(XPL(I),DIP)                                      
          C00(I)=C00(I)*DIP
C          WRITE(30,'(2F25.20)') C00(I)                                 
004    CONTINUE                                                       
C                                                                       
C.....Create potential grid for repulsive state.                        
C.....Zero out potential and scale diagonal elements. 
      OPEN(22, FILE='excited_PES.txt')                   
       DO 006 I=1,NX                                                    
          VVN(I)=0.0D0                                                  
          CALL ESTATE(XPL(I),SCRV2(I))                                  
          IF(SCRV2(I).GT.VMAX) SCRV2(I)=VMAX                           
          VVN(I)=(SCRV2(I)-E0)/ESC                                     
          WRITE(22,*) XPL(I), SCRV2(I)
006    CONTINUE                                                         
C                                                                       
      CALL PROP(IWRIT,PHIJ,C00,C,X1A,CR,TVN,VVN,NX,ESC,
     &          NDAMP,FF1,TTIME,E0,CD,CCD,TIMEFS,DTIME,NRECT, 
     &          ICHEX,INORM,IUNITC,IUNITN)                           
C                                                                       
C.....Calculate absolute value of autocorrelation function.             
C.....Need to check first and last values to determine whether          
C.....total time of propagation is sufficient.                          
       IF(IAUTOC.EQ.1) THEN                                             
          WRITE(IUNITF,*) 'Absolute value of Autocorrelation function' 
          WRITE(IUNITF,*) 'Time elapsed vs absolute value'             
          WRITE(IUNITF,*) ' '                                           
          DO 060 I=1,NRECT+1                                           
             CALL AUTOCRN(CCD(I),AUTOC(I))                             
             WRITE(IUNITF,*) TIMEFS(I),AUTOC(I)                         
060       CONTINUE                                                     
          TEMP1=0.00001D0*AUTOC(1)                                      
          IF(AUTOC(NRECT+1).GT.TEMP1) THEN                              
             WRITE(6,*) 'Final value of Autocorrelation function
     & too high'                             
             WRITE(6,*) 'Insufficient total propagation time'        
             WRITE(6,*) 'Increase value of TTIME'                     
             STOP                                                    
          END IF                                                      
       ELSE                                                             
          CALL AUTOCRN(CCD(1),TEMP2)                                  
          CALL AUTOCRN(CCD(NRECT+1),TEMP3)                            
          TEMP4=0.00001D0*TEMP2                                       
          IF(TEMP3.GT.TEMP4) THEN                                    
             WRITE(6,*) 'Final value of Autocorrelation function
     & too high'             
             WRITE(6,*) 'Insufficient total propagation time'       
             WRITE(6,*) 'Increase value of TTIME'                   
             STOP                                                      
          END IF                                                      
       END IF                                                           
C                                                                       
C.....Calculate and plot total absorption cross section.                
      CALL ABSSPEC(DTIME,NRECT,CCD,EIG,IUNITA,IWRIT)
C                                                                       
      STOP                                                              
      END PROGRAM                                               
C                                                                       
C--------------------------------------------------------------------
C               
      SUBROUTINE PROP(IWRIT,PHIJ,C00,C,X1A,CR,TVN,VVN,NX,ESC,
     &                NDAMP,FF1,TTIME,E0,CD,CCD,TIMEFS,DTIME,NRECT,
     &                ICHEX,INORM,IUNITC,IUNITN)                   
C                                                                       
C.....This performs the basic time propagation of the wavepacket.       
C.....X1A is the wavepacket on the excited dissociative surface.        
C     
C********************************************************************
C                                                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMPLEX*16 C00(NX),C(NX)                                          
      COMPLEX*16 PHIJ(NX)                                               
      COMPLEX*16 X1A(NX),CR(NX)                                         
      COMPLEX*16 CD                                                     
      COMPLEX*16 CCD(NRECT+1)                                           
      COMPLEX*16 EYE                                                    
      COMPLEX*16 CHEX,CHEX0                                             
      COMPLEX*16 CPF                                                    
      COMPLEX*16 CZERO                                                  
      DIMENSION VVN(NX),CBESS(0:255)                                    
      DIMENSION TVN(NX)
      DIMENSION TIMEFS(NRECT+1)                                         
C                                                                       
C.....Temporary arrays needed for Fast Fourier Transform routine FOUR1. 
      COMPLEX*16 FF1(NX)
C                                                                       
C.....Conversion factors and constants.                                 
C      DATA HTEV/27.211648D0/       
      CZERO=DCMPLX(0.D0,0.D0)                                           
      EYE=DCMPLX(0.D0,1.D0)                                             
C                                                                       
C.....Initialise wavefunction |C>.                                      
       DO 008 I=1,NX                                                    
          C(I)=C00(I)                                                
008    CONTINUE                                                         
C                                                                       
C.....NRECT : number of time steps to be computed.                      
      WRITE(6,'('' Number of time steps for propagation (NRECT) = '',I3
     & )') NRECT                                                 
C                                                                       
      WRITE(6,'('' Total time over which wavepacket is propagated'',   
     &'' (TTIME) = '',F10.5)') TTIME                                    
C.....Calculate length of each time step.                               
      DTIME=TTIME/DFLOAT(NRECT)                                         
      WRITE(6,'('' Length of each time step (DTIME)= '',F11.8)') DTIME  
C                                                                       
C.....Compute initial autocorrelation function.                         
      CALL ACEV(C00,C,CD,NX)                                            
C.....Set start time to zero.                                           
      TIME0=0.D0                                                        
C.....Save initial autocorrelation and corresponding time.              
      TIMEFS(1)=TIME0                                                   
      CCD(1)=CD                                                         
       IF(IWRIT.GT.0) THEN                                            
          WRITE(6,'('' Initial autocorrelation = '',F11.8,
     & '',''F11.8)') CD                                               
       END IF                                                        
C                                                                       
C.....Pre-calculate Bessel functions.                                   
C.....Calculate argument.                                               
      ALPHA=ESC*DTIME                                                   
       IF(IWRIT.EQ.1) THEN                                              
          WRITE(6,'('' Argument for Bessel functions (ALPHA) = '',
     & F11.8)') ALPHA                                                  
       END IF                                                           
C.....Calculate NTERMS : number of terms in the Chebychev polynomial    
C.....                   expansion.                                                       
      CBESS(0)=BESSJ0(ALPHA)                                            
      CBESS(1)=BESSJ1(ALPHA)                                            
      N=1                                                               
300    N=N+1                                                            
       IF(N.GT.256) THEN                                               
          WRITE(6,*) 'NTERMS greater than 256'                          
          WRITE(6,*) 'Dimensions of CBESS(N) must be increased'         
          STOP                                                         
       END IF                                                          
      CBESS(N)=BESSJ(CBESS(0),CBESS(1),N,ALPHA)                       
       IF(DABS(CBESS(N)).GT.1.0D-15) GO TO 300                        
      NTERMS=N                                                          
       IF(IWRIT.EQ.1) THEN                                             
          WRITE(6,'('' Number of terms in Chebychev expansion (NTERMS)
     & = '',I3)') NTERMS                                         
       END IF                                                            
      C0= CBESS(0)                                                      
      C1=2.0D0*CBESS(1)                                                 
C                                                                       
C.....Calculate phase factor.                                           
C.....E0 is Kosloff's VMIN+DELTA(E)/2, midpoint of energy range.        
      CPF=CDEXP(-EYE*E0*DTIME)                                          
      E0SC=E0/ESC                                                       
       IF(IWRIT.EQ.1) THEN                                              
          WRITE(6,'('' Midpoint of energy range (E0) =  '',F11.8)') E0
          WRITE(6,'('' Scaling factor (ESC) =  '',F11.8)') ESC          
          WRITE(6,'('' E0SC (E0/ESC) =  '',F11.8)') E0SC                
          WRITE(6,'('' Phase factor (CPF) = '',F11.8,'',''F11.8)') CPF 
       END IF                                                           
C                                                                       
C.....Compute PHIJ> = I H |C>                                           
C.....Calculate damping factor (DAMPF) for complex damping potential.   
      DAMPF=1.0D0/DFLOAT(NX-NDAMP)                                      
      CALL HCAFI(EYE,TVN,VVN,C,NX,PHIJ,ESC,FF1,DAMPF,NDAMP,1)  
C                                                                       
C.....Compute expectation value of energy (CHEX).                       
      CALL HEX2FI(EYE,C,PHIJ,CHEX0,NX)                                  
      CALL WAVENORM(C,NX,RNORM)
      CHEX=ESC*(E0SC*RNORM+CHEX0)                                       
       IF(IWRIT.GT.0) THEN                                         
          WRITE(6,*)' '                                              
          WRITE(6,'('' Energy before time propagation = '',F11.8,    
     & '',''F11.8)') CHEX                                              
       END IF                                                          
C                                                                       
       IF(INORM.EQ.1) THEN                                              
          WRITE(IUNITN,*) 'Wavepacket norm data'                    
          WRITE(IUNITN,*) 'Time step index vs norm'            
          WRITE(IUNITN,*) ' '                              
       END IF                                                           
C                                                                       
C          ************************************************             
C          ********** MAIN TIME PROPAGATION LOOP **********             
C          ************************************************             
C                                                                       
C.....Loop over NREC time steps.                                        
C.....The large arrays C, PHIJ, CR and X1A are needed in this loop.     
C.....X1A is the wavefunction at time t which is being accumulated      
C.....in the Chebychev expansion.                                       
C.....The basic recursion formula is:                                   
C.....|PHIJ> = -2.*I*H*|C> + |CR>                                      
C.....|X1A> = |X1A> + BESSJ(N)*|PHIJ>                                  
C       
      OPEN(31, FILE='dynamics.txt')                               
      JFLAG=0                                                           
       DO 001 NREC=1,NRECT                                           
          WRITE(31,*)
C                                                  
          TIMET=DFLOAT(NREC)*DTIME                                   
          DO 010 I=1,NX                                             
             X1A(I)=(C0*C(I)-C1*PHIJ(I))                             
             CR(I)=C(I)                                              
             C(I)=-PHIJ(I)                                           
010       CONTINUE                                                   
C                                                                       
C.....Implement Chebychev recursion formula.                            
          DO 014 N=2,NTERMS                                          
             CI=CBESS(N)*2.0D0                                       
             CALL HCAFI(EYE,TVN,VVN,C,NX,PHIJ,ESC,FF1,DAMPF,NDAMP,1)  
C                                                                       
             DO 012 I=1,NX                                           
                PHIJ(I)=-2.0D0*PHIJ(I)                                 
                PHIJ(I)=PHIJ(I)+CR(I)                                
                CR(I)=C(I)                                            
                C(I)=PHIJ(I)                                          
                X1A(I)=X1A(I)+CI*PHIJ(I)                              
012          CONTINUE                                                 
014       CONTINUE                                                    
C                                                                       
C.....We now have converged wavefunction (|X1A>) at time t.             
C.....Multiply by phase factor (CPF).                                   
          DO 015 I=1,NX                                               
             X1A(I)=X1A(I)*CPF                                        
015       CONTINUE                                                     
C                                                                       
C.....Compute autocorrelation function and save.                        
          CALL ACEV(C00,X1A,CD,NX)                                  
          TIMEFS(NREC+1)=TIMET                                        
          CCD(NREC+1)=CD                                             
C                                                                       
C.....Place |C> = |X1A> for next stage of propagation.                  
          DO 017 I=1,NX                                               
             C(I)=X1A(I)
             aaa = 0.50579718d0 + dfloat(I-1)*12.04279/dfloat(256)
             bbb = sqrt(dimag(X1A(I))**2+real(X1A(I))**2)           
             IF(aaa.LT.10.d0) WRITE(31,'(3F25.20)') TIMET, aaa, bbb
017       CONTINUE
C                                                                       
C.....Check norm of wavefunction to see that it is fully absorbed       
C.....by damping function.                                              
          IF(INORM.EQ.1) THEN                                         
             CALL WAVENORM(C,NX,RNORM)                                
             WRITE(IUNITN,*) NREC,RNORM                               
          END IF                                                      
C                                                                       
C.....Compute |PHIJ> = I*H*|C>                                          
          CALL HCAFI(EYE,TVN,VVN,C,NX,PHIJ,ESC,FF1,DAMPF,NDAMP,1)  
C                                                                       
C.....Evaluate expectation value of Hamiltonian:                        
C.....CHEX = <X1A|H|X1A> / <X1A|X1A>                                    
          IF(ICHEX.EQ.1) THEN                                       
             CALL HEX2FI(EYE,C,PHIJ,CHEX0,NX)                          
             CHEX=ESC*(CHEX0+E0SC*RNORM)                              
             WRITE(IUNITC,*)' Convergence energy (CHEX) = ',CHEX     
          END IF                                                       
C                                                                       
001    CONTINUE                                                       
C                                                                       
C         **************************************************            
C         ********** END OF TIME PROPOGATION LOOP **********            
C         **************************************************            
C                                                                       
      RETURN                                                            
      END SUBROUTINE PROP                                           
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE UTGR(VI,VO,N)                                          
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION VI(N),VO(N)                                             
      N2=N/2                                                            
      NS=N2                                                             
       DO 019 I=1,N2                                                    
          VO(I)=VI(I+NS)                                           
          VO(I+NS)=VI(I)                                          
019    CONTINUE                                                         
      RETURN                                                            
      END SUBROUTINE UTGR                                          
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE GROUND(NX,ZMA,ZMB,ZMU,ZLZ,IVS,SCRV1,SCRM1,SCRV2,       
     &                  SCRV3,SCRV4,SCRM2,VMAX,RO,X0,DZ,CZERO,C00,IWRIT)
C.....Establishes a wavefunction for propagation on the excited state.  
C.....Fragment masses defined by ZMA and ZMB. IVS is the vibrational    
C.....state selected.                                                   
C
C********************************************************************
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
C.....Declare arrays.                                                   
C.....C00 contains initial wavefunction.                                
      COMPLEX*16 C00(NX)                                                
      COMPLEX*16 CZERO                                                  
C                                                                       
C.....Temporary storage.                                                
      DIMENSION SCRV1(NX),SCRV2(NX),SCRV3(NX),SCRV4(NX)                 
      DIMENSION SCRM1(NX,NX),SCRM2(NX,NX)                               
C     DATA CZERO/(0.0D0,0.0D0)/                                         
C                                                                       
C.....Calculate reduced mass.                                           
      ZMU=ZMA*ZMB/(ZMA+ZMB)                                             
      WRITE(6,'('' Reduced mass of AB (ZMU) = '',F14.8)') ZMU          
      WRITE(6,'('' Vibrational state selected (IVS) = '',I2)') IVS     
C                                                                       
C.....Create potential grid for ground state.                           
      OPEN(21, FILE='ground_PES.txt')
      XX=X0                                                             
       DO 020 I=1,NX                                                    
          CALL GSTATE(XX,SCRV4(I))
          IF(SCRV4(I).GT.VMAX) SCRV4(I)=VMAX
          WRITE(21,*) XX,SCRV4(I)
          XX=XX+DZ                                                   
020    CONTINUE                                                         
C                                                                       
C.....Calculate initial wavefunction.                                   
C.....Use points up to 5.0*RO for evaluation of bound vib state.        
      NEND=INT((5.0D0*RO-X0)*DFLOAT(NX)/ZLZ)                            
       IF(NEND.GT.NX) THEN                                              
          WRITE(6,*)' Increase number of grid points for FGH1D'       
          STOP                                                     
       END IF                                                           
C.....Check that NEND is even.                                          
      ITEST=MOD(NEND,2)                                                 
       IF(ITEST.NE.0) THEN                                              
          NEND=NEND+1                           
       END IF                                                           
       IF(IWRIT.GT.0) THEN                                           
          WRITE(6,'( '' Number of points for wavefunction evaluation
     & (NEND) =  '', I3)') NEND                             
       END IF                                                        
C.....Calculate corresponding length for NEND.                          
      ZFGH=ZLZ*DFLOAT(NEND)/DFLOAT(NX)                                  
       IF(IWRIT.GT.0) THEN                                              
          WRITE(6,'('' Corresponding grid length (ZFGH) =  '',F11.8)')
     & ZFGH                                      
       END IF                                                      
C                                                                       
      CALL FGH1D(NX,NEND,ZMU,ZFGH,SCRV4,SCRM1,SCRV2,SCRV3,              
     &           SCRV1,SCRM2,IWRIT)                                 
C.....SCRM1 contains the initial wavefunction.                          
C.....SCRV2 contains the eigenvalues for the vibrational states.        
C                                                                       
C.....Zero out large separation part of the eigenvector arrays.         
       IF(NX.GT.NEND) THEN                                           
          DO 022 I=1,NX                                               
             DO 021 J=(NEND+1),NX                                    
                SCRM1(J,I)=0.0D0                                     
021          CONTINUE                                                 
022       CONTINUE                                                    
       END IF                                                           
C                                                                       
C.....Zero out C00 and fill with selected ground vibrational state.     
       DO 023 I=1,NX                                                   
          C00(I)=CZERO                                                  
023    CONTINUE                                                        
       DO 025 I=1,NEND                                                 
          C00(I)=DCMPLX(SCRM1(I,IVS+1),0.0D0)                           
025    CONTINUE                                                        
C.....Zero out wavefunction at edges of grid.                           
      C00(1)=CZERO                                                     
      C00(2)=CZERO                                                     
      C00(NX-1)=CZERO                                                  
      C00(NX)=CZERO                                                    
C                                                                       
C.....C00 now contains the initial wavefunction.                        
C                                                                       
      RETURN                                                            
      END SUBROUTINE GROUND                              
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE ESTATE(XAA,YYA)                                        
C.....This subroutine calculates the A1Pi HCl potential energy curve.   
C.....(Ref: Simon C. Givertz, Ph.D thesis (1986))                       
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
C.....Data for potential energy calculation.                            
      DATA A/40.57/                                                     
      DATA B/-2.995/                                                    
      DATA C/2.801/                                                     
      DATA D/2.243D-13/                                                 
      DATA ALPHA/1.811/                                                 
C                                                                       
C.....Evaluate potential energy on grid.                                
      EX=DEXP(-ALPHA*XAA)                                             
      YYA=A*(1.0+B/XAA+C/XAA**2+D/XAA**3)*EX                          
       IF(XAA.GT.6.5D0) YYA=0.0D0                                      
C                                                                       
      END SUBROUTINE ESTATE                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE GSTATE(XAA,X1S)                                        
C.....Subroutine to calculate the X1Sigma HCl (ground state) attractive 
C.....potential using the truncated polynomial presented by J.F.Olgilvie
C.....Proc.R.Soc.Lond. A 378 287-300 (1981)                             
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
C                                                                       
C.....Set up array for coefficients and fill.                           
      DIMENSION HHCL(9)                                                 
      DATA HHCL /0.961914D0,-1.362999D0,0.86675D0,-0.49804D0,           
     & 0.1727D0,0.2687D0,-1.977D0,2.78D0,4.89D0/                     
C                                                                       
C.....RE is the equilibrium bond length in HCl.                         
      DATA RE /2.408558D0/                                              
C                                                                       
C.....Calculate potential on the grid.                                  
      OMEGA=2.0D00*(XAA-RE)/(XAA+RE)                                   
      HCL=0.0D00                                                       
       DO 028 K=2,9                                                    
          HCL=HCL+HHCL(K)*(OMEGA**(K-1))                                
028    CONTINUE                                                        
      HCL=HHCL(1)*OMEGA*OMEGA*(HCL+1.0D00)                             
      HCL=HCL-0.169695D+00                                             
      X1S=HCL                                                          
C                                                                       
C.....Real potential rises to zero after the well.                      
       IF(XAA.GE.4.0D0) THEN                                            
          RX=XAA-4.0D0                                                
          X1S=-4.85D-2*((1.0D0-DTANH(RX))**1.5D0)                    
       END IF                                                           
       IF(XAA.GE.6.5D0) X1S=0.0D0                                       
C                                                                       
      RETURN                                                            
      END SUBROUTINE GSTATE                                         
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE HCAFI(EYE,TVN,VVN,C00,NX,HCAV,ESC, 
     &                 FF1,DAMPF,NDAMP,IDAMPTYPE)            
C.....This subroutine acts with the Hamiltonian on a wavefunction       
C.....represented on a two dimensional grid and then multiplies by i    
C.....|HCAV> = I * H * |C00>                                          
C.....VVN contains the potential evaluated on the grid.               
C.....TVN is a vector and contains the kinetic energy at the grid     
C.....points in momentum space.                                       
C.....IDAMPTYPE is a damping option:                                    
C.....IDAMPTYPE=0  no damping.                                     
C.....IDAMPTYPE=1  damping with complex potential.                 
C     
C********************************************************************
C                                                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMPLEX*16 EYE,C00(NX),HCAV(NX)                                   
      COMPLEX*16 FF1(NX)
      COMPLEX*16 CZERO                                                  
      DIMENSION VVN(NX),TVN(NX)                                         
      DATA CZERO/(0.0D0,0.0D0)/                                         
C                                                                       
C.....Data for complex damping potential.                               
      DATA ZMIN/-0.3346D0/                                              
C                                                                       
C.....Evaluate action of Hamiltonian on each of the wavepackets.        
C.....Array FF1 will contain the NX Fourier Transforms of array C00.    
      ISGN=+1                                                          
       DO 061 I=1,NX                                                   
          FF1(I)=C00(I)                                                 
061    CONTINUE                                                        
      CALL FOUR1(FF1,NX,ISGN)                                          
       DO 029 I=1,NX                                                   
          FF1(I)=FF1(I)*TVN(I)                                          
029    CONTINUE                                                        
      ISGN=-1                                                          
      CALL FOUR1(FF1,NX,ISGN)                                          
C                                                                       
C.....Evaluate potential energy term and apply damping function.        
       IF((IDAMPTYPE.NE.0).AND.(IDAMPTYPE.NE.1)) THEN                
          WRITE(6,*) 'DAMPTYPE variable not allowed'                  
          WRITE(6,*) 'Change call to HCAFI'                         
          STOP                                                      
       END IF                                                         
C                                                                       
C.....Apply complex damping potential.                                  
       DO 032 I=1,NX                                                
          IF((IDAMPTYPE.EQ.1).AND.(I.GE.NDAMP)) THEN               
             DF=DFLOAT(I-NDAMP)*DAMPF                                
             DF=DF**3                                                  
             DF=2.0D0*DF*ZMIN                                           
C.....Apply scaling factor.                                            
             DF=DF/ESC                                                 
          END IF                                                      
          IF((IDAMPTYPE.EQ.1).AND.(I.GE.NDAMP)) THEN                    
             HCAV(I) = DCMPLX(VVN(I),DF)*C00(I)                      
          ELSE                                                          
             HCAV(I) = (VVN(I)*C00(I))                                
          END IF                                                        
032    CONTINUE                                                      
C                                                                       
C.....|HCAV> = |HCAV> + |FF1>                                          
       DO 031 I=1,NX                                                 
          HCAV(I) = HCAV(I) + FF1(I)                              
031    CONTINUE                                                      
C                                                                       
C.....Now multiply HCAV by i (=EYE).                                    
       DO 033 I=1,NX                                                 
          HCAV(I)=HCAV(I)*EYE                                        
033    CONTINUE                                                     
C.....|HCAV> = i * H * |PSI>                                        
C                                                                       
C.....Zero out first and last two grid points.                          
      HCAV(1)=CZERO                                                    
      HCAV(2)=CZERO                                                    
      HCAV(NX-1)=CZERO                                                 
      HCAV(NX)=CZERO                                                   
C                                                                       
      RETURN                                                            
      END SUBROUTINE HCAFI
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE HEX2FI(EYE,C00,HCAV,CHEX,NX)                           
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMPLEX*16 EYE,C00(NX),HCAV(NX),CHEX                              
      CHEX=DCMPLX(0.D0,0.D0)                                            
C                                                                       
       DO 034 I=1,NX                                                    
          CHEX=CHEX + DCONJG(C00(I))*HCAV(I)                            
034    CONTINUE                                                         
      CHEX = -EYE*CHEX                                                  
      RETURN                                                            
      END SUBROUTINE HEX2FI                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE WAVENORM(CC,NX,RNORM)                                  
C.....Calculates the norm of the wavefunction.                          
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMPLEX*16 CC(NX)                                                 
C                                                                       
      RNORM=0.0D0                                                    
       DO 036 I=1,NX                                                
          TEMP = DREAL(DCONJG(CC(I))*CC(I))                          
          RNORM = RNORM + TEMP                                      
036    CONTINUE                                                    
C                                                                       
      RETURN                                                            
      END SUBROUTINE WAVENORM                                           
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE ACEV(C0,C,CD,NX)                                       
C.....This subroutine calculates the autocorrelation function.          
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMPLEX*16 C0(NX),C(NX),CD                                        
      CD=DCMPLX(0.D0,0.D0)                                              
       DO 038 I=1,NX                                                  
          CD=CD+DCONJG(C0(I))*C(I)                                     
038    CONTINUE                                                       
      RETURN                                                            
      END SUBROUTINE ACEV                                            
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE DIPMOM(R,DIP)                                          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)   
      DATA E,F,G/0.279D0,-0.905D0,1.029D0/                              
      DATA beta,RE/0.687D0,2.555D0/                                     
      X=beta*(R-RE)                                                     
      Y=TANH(X)                                                         
      Y=1.0D0-Y                                                         
      Y2=Y*Y                                                            
      Y3=Y*Y2                                                           
      DIP=E*(Y+F*Y2+G*Y3)                                               
C.....Set min value for DIP.                                            
       IF(DIP.LT.0.0D0) DIP=0.0D0                                    
      RETURN                                                            
      END SUBROUTINE DIPMOM                                        
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE AUTOCRN(CTEMP,ATC)                                     
C.....Calculates absolute value of the autocorrelation function.        
C                                                                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      COMPLEX*16 CTEMP,VAR1                                             
C                                                                       
      VAR1=DCONJG(CTEMP)*CTEMP                                          
      ATC=CDABS(VAR1)                                                   
C                                                                       
      RETURN                                                            
      END SUBROUTINE AUTOCRN                                       
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE ABSSPEC(DT,NR,FAC,EIG,IUNITA,IWRIT) 
C.....This subroutine calculates the absorption spectrum from the real  
C.....part of the half Fourier transform of the autocorrelation function
C     
C************************************************************************
C                                                                  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      PARAMETER (NSIZE=4096)                                            
      COMPLEX*16 FACU(NSIZE),FAC(NSIZE)                                 
      DIMENSION FF(NSIZE),YY(NSIZE)                                     
C                                                                       
      DATA PI/3.141592653589793D0/                                      
      DATA CZERO/(0.0D0,0.0D0)/                                         
      CLIGHT=29.97925D0*2.418884D0/0.529177D0                           
      TOTZEIT=DFLOAT(NSIZE)*DT                                          
      DF=2.0D0*PI/TOTZEIT                                               
      FACTOR=2.0*PI*TOTZEIT/(3*CLIGHT*DSQRT(DFLOAT(NSIZE)))             
C                                                                       
C.....Conversion factors : Hartrees --> cm-1, a.u. --> Angstroms 2      
      HTCM=2.194753236D5                                                
      ATAS=2.799647996D-1                                               
C                                                                       
      WRITE(6,*)' '                                                   
      WRITE(6,*)'* Absorption cross section (FT analysis) *'          
      WRITE(6,*)' '                                                   
      IF(IWRIT.GT.0) THEN                                               
         WRITE(6,'('' Total number of time steps computed = '',
     & I3)') NR                                                     
         WRITE(6,'('' Length of each time step = '',F19.16)') DT     
         WRITE(6,'('' Total time represented on FT grid = '',
     & F19.12)') TOTZEIT                                               
         WRITE(6,'('' Maximum energy (cm-1) represented on FT grid = '',
     & F19.12)') DF*DFLOAT(NSIZE)*HTCM                                
       END IF                                                         
C                                                                       
      N=2048                                                            
      NH=N/2                                                            
      N=2*N                                                             
      NRP1=NR+1                                                         
      NRP2=NR+2                                                         
C                                                                       
C.....EIG is the energy level (Eigenvalue) for the vibrational state    
C.....selected (obtained from the FGH program).                         
C                                                                       
C.....Put FAC in FACU.                                                  
       DO 041 I=1,NRP1                                                 
          FACU(I)=FAC(I)                                                
041    CONTINUE                                                        
C                                                                       
C.....Zero out empty elements of array, adding large numbers of zeros.  
C.....This has the effect of making the absorption spectrum smoother.   
C                                                                       
       DO 042 I=NRP2,N                                                 
          FACU(I)=CZERO                                                 
042    CONTINUE                                                        
C                                                                       
C.....Fourier transform coefficients.                                   
C                                                                       
      FACU(1)=0.5D0*FACU(1)                                             
      CALL FOUR1(FACU,N,1)                                              
C                                                                       
C.....Calculate energies.                                               
      ABSMAX=0.0D0                                                      
      FMIN=0.0D0                                                        
       DO 043 I=1,N                                                    
          ETOT=FMIN+DFLOAT(I-1)*DF                                      
          EFREQ=ETOT-EIG                                                
          FF(I)=EFREQ*HTCM                                              
          ZMFAC=DREAL(FACU(I))                                          
          ZMFAC=FACTOR*ATAS*EFREQ*ZMFAC                                 
C.....Multiply by degeneracy factor.                                    
          YY(I)=2.0D0*ZMFAC                                             
          IF(YY(I).GT.ABSMAX) THEN                                     
             ISAVE=I                                                  
             ABSMAX=YY(I)                                           
          END IF                                                  
043    CONTINUE                                                        
C                                                                       
      WRITE(6,'('' Maximum absorption cross section (Ang.sqrd) = '',   
     & F11.8)') ABSMAX                                               
      WRITE(6,'('' Corresponding photon energy (cm-1) = '',            
     & F11.4)') FF(ISAVE)                                           
C                                                                       
C..... Write out cross section values if they are greater than 0.1%  of 
C..... maximum value.                                                   
      ABSMAX=0.001D0*ABSMAX                                             
      WRITE(IUNITA,*)'Total absorption cross section data'              
      WRITE(IUNITA,*)'Photon energy (cm-1) vs cross section (Ang.sqrd.)'
      WRITE(IUNITA,*)' '                                                
      WRITE(6,*)'Total absorption cross section data'                   
      WRITE(6,*)'Photon energy (cm-1) vs cross section (Ang.sqrd.)'     
      WRITE(6,*)' '                                                     
       DO 045 I=1,N                                                     
          IF(DABS(YY(I)).GT.ABSMAX) THEN                           
             WRITE(IUNITA,*) FF(I),YY(I)                             
             WRITE(6,*) FF(I),YY(I)                                
          END IF                                                    
045    CONTINUE                                                       
C                                                                       
      RETURN                                                            
      END SUBROUTINE ABSSPEC                                            
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE FGH1D(NDIM,NX,ZMU,ZL,VV,ZR,WCH,FV1,FV2,AR,IWRIT)       
C.....PROGRAM TO PERFORM FGH CALCULATION TO SET UP A 1 DIMENSIONAL      
C.....WAVEFUNCTION ON A GRID.                                           
C.....BASED ON C.C. MARSTON AND G.G. BALINT-KURTI, J. CHEM. PHYS.,      
C.....91, 3571 (1989).                                            
C.....THIS PROGRAM USES AN EVEN NUMBER OF GRID POINTS, AND THE          
C.....THEORY OF THE ABOVE PAPER HAS BEEN EXTENDED TO ACCOMMODATE        
C.....THIS.THE PROGRAM IS TO BE USED EITHER TO                          
C.....SET UP AN INITIAL VIBRATIONAL WAVEFUNCTION FOR A TIME DEPENEDENT  
C.....QUANTUM DYNAMICAL CALCULATION, OR IT MAY BE USED TO SET UP THE    
C.....COEFFICIENT ARRAY FOR THE FINAL STATE ANALYSIS.                   
C.....THE FOLLOWING ARRAYS REPRESENT :-                               
C.....NDIM -- DIMENSION OF ARRAYS.                                      
C.....AR   -- THE HAMILTONIAN MATRIX                                      
C.....NX   -- NUMBER OF POINTS IN GRID, THIS IS 1 GREATER THAN THAT       
C.....        USED IN THE REST OF THE PROGRAM TO MAKE IT ODD.           
C.....ZMU  -- REDUCED MASS.                                              
C.....ZL   -- LENGTH OF INTERVAL (ONE INCREMENT LONGER THAN IN MAIN       
C.....        PART OF CODE.                                               
C.....VV   -- POTENTIAL VECTOR.                                           
C.....ZR   -- THE EIGENVECTOR (X,Y) : X = WAVEFUNCTION                    
C.....                                Y = THE ENERGg LEVEL                
C.....WCH  -- THE EIGENVALUES FOR RESPECTIVE ENERGY LEVELS                
C.....FV1  -- (2/NX)*T(L)  WHERE L = AN INTEGER                      
C.....                          NX = NUMBER OF GRID POINTS          
C.....FV2  -- COS((2*PI*L)/NX)                                       
C
C********************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION AR(NDIM,NX),WCH(NX),ZR(NDIM,NX),VV(NX)                  
      DIMENSION FV1(NX),FV2(NX)                                         
C.....THE FOLLOWING DATA CONSTANTS REPRESENT :-                       
      DATA PI/3.141592653589793D0/                                      
C.....NPRIN - 0 OR 1 FOR EIGENVALUE OR EIGENVECTOR RESPECTIVELY 
      DATA NPRIN/1/                                                     
C.....TITLE FOR DATA TO BE PRINTED                                    
       IF(IWRIT.GE.10) THEN                                         
          WRITE(6,*)'LOWEST 10 ENERGY LEVELS FOR THE DIATOMIC FRAGMENT'
       END IF                                                       
C.....TEST THAT NX IS EVEN.                                           
      ITEST=MOD(NX,2)                                                   
       IF(ITEST.NE.0) THEN                                           
          WRITE(6,*) '**** NX MUST BE EVEN __ FATAL ERROR  ****'     
         STOP                                                           
       END IF                                                        
      NHALF=NX/2                                                        
      NHAM1=NHALF-1                                                     
C.....NOW SET UP HAMILTONIAN MATRIX                                   
      DARG=2.0D0*PI/DFLOAT(NX)                                          
      TARG=4.0D0*((PI/ZL)**2.0D0)/(ZMU*DFLOAT(NX))                      
C.....PRECALCULATE KINETIC ENERGY FACTOR.                             
       DO 001 L=1,NX                                                
          FV1(L)=TARG*DFLOAT(L**2)                                 
          FV2(L)=DCOS(DARG*DFLOAT(L))                            
001    CONTINUE                                                    
C.....INITIALISE VARIABLES                                            
C.....CONST - ((-1)**(I-J))*T(NX/2)                                 
C.....EQUIJ - WHEN I=J THEN SUMMATION ONLY THAT OF T(L)             
C.....INIJ  - (I-J)                                                 
C.....SUM   - TOTAL SUMMATION OF EQUATION WITHIN SUM                
      CONST=0.0D0                                                       
      EQUIJ=0.0D0                                                       
C.....PRECALCULATE THE SUM OF THE T(L) WHEN I=J                       
       DO 002 L=1,NHAM1                                             
          EQUIJ=EQUIJ+FV1(L)                                      
002    CONTINUE                                                     
       DO 003 I=1,NX                                              
          DO 004 J=1,I                                       
             IJ=(I-J)                                             
             INIJ=IJ                                               
             SUM=0.0D0                                       
             IF(IJ.EQ.0) THEN                                   
                SUM=EQUIJ                                          
             ELSE                                                     
                DO 005 L=1,NHAM1                                  
                   SUM=SUM+FV1(L)*FV2(IJ)                        
                   IJ=IJ+INIJ                                      
C.....COSINE IS PERIODIC, SO ONLY NEED VALUE WITHIN ONE CYCLE/PERIOD  
                   IF(IJ.GT.NX) IJ=MOD(IJ,NX)                          
005             CONTINUE                                             
             END IF                                                 
C.....ADD IN T(NX/2)                                                  
             CONST=((-1)**INIJ)*FV1(NHALF)/2                         
             AR(I,J)=SUM+CONST                                      
004       CONTINUE                                                  
C.....ADD THE POTENTIAL VALUE WHEN KNONICLEAR DELTA FUNCTION          
C.....EQUALS 1, I.E. WHEN I AND J ARE EQUAL                        
          AR(I,I)=AR(I,I)+VV(I)                                     
003    CONTINUE                                                    
C.....NOW FILL OUT HAMILTONIAN MATRIX.                                
       DO 006 I=1,NX                                                
          DO 007 J=1,I                                               
             AR(J,I)=AR(I,J)                                         
007       CONTINUE                                                    
006    CONTINUE                                                     
C.....NOW CALL EIGENVALUE SOLVER                                      
      CALL RS(NDIM,NX,AR,WCH,NPRIN,ZR,FV1,FV2,IERR)                     
      I=30                                                              
C.....NOW WRITE OUT LOWEST 10 EIGENVALUES AND EIGENVECTORS.           
       DO 009 I=1,10                                                  
          IF(IWRIT.GE.10) THEN                                         
             WRITE(6,*)' ENERGY LEVEL NO ',I,' EIGENVALUE= ',WCH(I)   
             IF(NPRIN.EQ.1) THEN                                     
                DO 010 J=1,NX                                        
                   WRITE(6,*) J,' WAVEFUNCTION=',ZR(J,I)             
010             CONTINUE                                              
             END IF                                                    
          END IF                                                      
009    CONTINUE                                                       
      RETURN                                                            
      END SUBROUTINE FGH1D                                            
C                                                                       
C-----------------------------------------------------------------------                
C                                              
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)                       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)                      
C                                                                       
C.....THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF                 
C.....SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)     
C.....TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)             
C.....OF A REAL SYMMETRIC MATRIX.                                       
C                                                                       
C.....ON INPUT-                                                         
C                                                                       
C.....NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL    
C.....ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM            
C.....DIMENSION STATEMENT,                                           
C                                                                       
C.....N  IS THE ORDER OF THE MATRIX  A,                              
C                                                                        
C.....A  CONTAINS THE REAL SYMMETRIC MATRIX,                         
C                                                                       
C.....MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF              
C.....ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO          
C.....ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.    
C                                                                       
C.....ON OUTPUT-                                                        
C                                                                       
C.....W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,                
C                                                                       
C.....Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,              
C                                                                       
C.....IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN            
C.....ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE           
C.....DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,            
C                                                                       
C.....FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.                   
C                                                                       
C.....QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C.....APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C********************************************************************
C                                                                       
       IF(N.LE.NM) GO TO 10                                           
          IERR = 10 * N                                               
          GO TO 50                                                    
C                                                                       
   10  IF(MATZ.NE.0) GO TO 20                                      
C*************** FIND EIGENVALUES ONLY **********                       
      CALL TRED1(NM,N,A,W,FV1,FV2)                                     
      CALL TQLRAT(N,W,FV2,IERR)                                        
       GO TO 50                                                     
C*************** FIND BOTH EIGENVALUES AND EIGENVECTORS **********      
   20 CALL TRED2(NM,N,A,W,FV1,Z)                                       
      CALL TQL2(NM,N,W,FV1,Z,IERR)                                     
   50 RETURN                                                            
C*************** LAST CARD OF RS **********                             
      END SUBROUTINE RS
C                                                
C--------------------------------------------------------------------
C                                                                       
      SUBROUTINE TRED1(NM,N,A,D,E,E2)                                   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION A(NM,N),D(N),E(N),E2(N)                                 
C                                                                       
C.....THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,    
C.....NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   
C.....HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   
C                                                                       
C.....THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX                   
C.....TO A SYMMETRIC TRIDIAGONAL MATRIX USING                           
C.....ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            
C                                                                       
C.....ON INPUT-                                                         
C                                                                       
C.....NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C.....ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C.....DIMENSION STATEMENT,                                         
C                                                                       
C.....N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C.....A CONTAINS THE REAL SYMMETRIC INPUT MATRIX. ONLY THE          
C.....LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               
C                                                                       
C.....ON OUTPUT-                                                        
C                                                                       
C.....A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-             
C.....FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER         
C.....TRIANGLE. THE FULL UPPER TRIANGLE OF A IS UNALTERED,        
C                                                                       
C.....D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,    
C                                                                       
C.....E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         
C.....MATRIX IN ITS LAST N-1 POSITIONS. E(1) IS SET TO ZERO,      
C                                                                       
C.....E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    
C.....E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.        
C                                                                       
C.....QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C.....APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C********************************************************************
C                                                                       
       DO 100 I = 1, N                                               
  100     D(I) = A(I,I)                                               
C*************** FOR I=N STEP -1 UNTIL 1 DO -- **********               
       DO 300 II = 1, N                                             
          I = N + 1 - II                                        
          L = I - 1                                              
          H = 0.0D0                                                
          SCALE = 0.0D0                                            
          IF(L.LT.1) GO TO 130                                        
C*************** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********       
          DO 120 K = 1, L                                           
  120     SCALE = SCALE + DABS(A(I,K))                              
C                                                                       
          IF(SCALE.NE.0.0D0) GO TO 140                                
  130     E(I) = 0.0D0                                               
          E2(I) = 0.0D0                                           
          GO TO 290                                                
C                                                                       
  140     DO 150 K = 1, L                                         
             A(I,K) = A(I,K) / SCALE                               
             H = H + A(I,K) * A(I,K)                                 
  150     CONTINUE                                                   
C                                                                       
          E2(I) = SCALE * SCALE * H                                 
          F = A(I,L)                                                
          G = -DSIGN(DSQRT(H),F)                                    
          E(I) = SCALE * G                                          
          H = H - F * G                                             
          A(I,L) = F - G                                           
          IF(L.EQ.1) GO TO 270                                        
             F = 0.0D0                                              
C                                                                       
          DO 240 J = 1, L                                         
             G = 0.0D0                                            
C*************** FORM ELEMENT OF A*U **********                         
             DO 180 K = 1, J                                       
  180        G = G + A(J,K) * A(I,K)                              
C                                                                       
             JP1 = J + 1                                          
             IF(L.LT.JP1) GO TO 220                            
C                                                                      
             DO 200 K = JP1, L                                       
  200        G = G + A(K,J) * A(I,K)                               
C*************** FORM ELEMENT OF P **********                           
  220        E(J) = G / H                                           
             F = F + E(J) * A(I,J)                                
  240     CONTINUE                                                  
C                                                                       
          H = F / (H + H)                                         
C*************** FORM REDUCED A **********                              
          DO 260 J = 1, L                                           
             F = A(I,J)                                             
             G = E(J) - H * F                                     
             E(J) = G                                             
C                                                                       
             DO 260 K = 1, J                                        
                A(J,K) = A(J,K) - F * E(K) - G * A(I,K)             
  260     CONTINUE                                                   
C                                                                       
  270     DO 280 K = 1, L                                           
  280     A(I,K) = SCALE * A(I,K)                                  
C                                                                       
  290     H = D(I)                                                  
          D(I) = A(I,I)                                           
          A(I,I) = H                                               
  300  CONTINUE                                                  
C                                                                       
      RETURN                                                            
C************** LAST CARD OF TRED1 **********                          
      END SUBROUTINE TRED1                                     
C                                                                       
C------------------------------------------------------------------
C                                                                       
      SUBROUTINE TRED2(NM,N,A,D,E,Z)                                    
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION A(NM,N),D(N),E(N),Z(NM,N)                               
C                                                                       
C.....THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,    
C.....NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   
C.....HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   
C                                                                       
C.....THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A              
C.....SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING               
C.....ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            
C                                                                       
C.....ON INPUT-                                                         
C                                                                       
C.....NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C.....ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C.....DIMENSION STATEMENT,                                         
C                                                                       
C.....N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C.....A CONTAINS THE REAL SYMMETRIC INPUT MATRIX. ONLY THE          
C.....LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               
C                                                                       
C.....ON OUTPUT-                                                        
C                                                                       
C.....D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,    
C                                                                       
C.....E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         
C.....MATRIX IN ITS LAST N-1 POSITIONS. E(1) IS SET TO ZERO,      
C                                                                       
C.....Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX                
C.....PRODUCED IN THE REDUCTION,                                   
C                                                                       
C.....A AND Z MAY COINCIDE. IF DISTINCT, A IS UNALTERED.            
C                                                                       
C.....QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C.....APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C********************************************************************
C                                                                       
       DO 100 I = 1, N                                               
          DO 100 J = 1, I                                            
             Z(I,J) = A(I,J)                                      
  100  CONTINUE                                                    
C                                                                       
       IF(N.EQ.1) GO TO 320                                           
C*************** FOR I=N STEP -1 UNTIL 2 DO -- **********               
       DO 300 II = 2, N                                             
          I = N + 2 - II                                             
          L = I - 1                                                  
          H = 0.0D0                                                 
          SCALE = 0.0D0                                             
          IF(L.LT.2) GO TO 130                                    
C*************** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********       
          DO 120 K = 1, L                                            
  120     SCALE = SCALE + DABS(Z(I,K))                               
C                                                                       
          IF(SCALE.NE.0.0D0) GO TO 140                          
  130     E(I) = Z(I,L)                                             
          GO TO 290                                                   
C                                                                       
  140     DO 150 K = 1, L                                           
             Z(I,K) = Z(I,K) / SCALE                                
             H = H + Z(I,K) * Z(I,K)                                
  150     CONTINUE                                                  
C                                                                       
          F = Z(I,L)                                                
          G = -DSIGN(DSQRT(H),F)                                    
          E(I) = SCALE * G                                           
          H = H - F * G                                             
          Z(I,L) = F - G                                             
          F = 0.0D0                                                
C                                                                       
          DO 240 J = 1, L                                           
             Z(J,I) = Z(I,J) / H                                    
             G = 0.0D0                                             
C*************** FORM ELEMENT OF A*U **********                         
             DO 180 K = 1, J                                        
  180        G = G + Z(J,K) * Z(I,K)                                
C                                                                       
             JP1 = J + 1                                           
             IF(L.LT.JP1) GO TO 220                             
C                                                                       
             DO 200 K = JP1, L                                    
  200        G = G + Z(K,J) * Z(I,K)                               
C*************** FORM ELEMENT OF P **********                           
  220        E(J) = G / H                                            
             F = F + E(J) * Z(I,J)                                  
  240     CONTINUE                                                 
C                                                                      
          HH = F / (H + H)                                          
C*************** FORM REDUCED A **********                              
          DO 260 J = 1, L                                           
             F = Z(I,J)                                              
             G = E(J) - HH * F                                     
             E(J) = G                                              
C                                                                       
             DO 260 K = 1, J                                         
                Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)              
  260     CONTINUE                                                  
C                                                                       
  290     D(I) = H                                                  
  300  CONTINUE                                                     
C                                                                       
  320  D(1) = 0.0D0                                                  
       E(1) = 0.0D0                                                 
C*************** ACCUMULATION OF TRANSFORMATION MATRICES **********     
       DO 500 I = 1, N                                               
          L = I - 1                                                 
          IF(D(I).EQ.0.0D0) GO TO 380                            
C                                                                       
          DO 360 J = 1, L                                            
             G = 0.0D0                                               
C                                                                       
             DO 340 K = 1, L                                          
  340        G = G + Z(I,K) * Z(K,J)                                  
C                                                                       
             DO 360 K = 1, L                                          
                Z(K,J) = Z(K,J) - G * Z(K,I)                          
  360     CONTINUE                                                    
C                                                                       
  380     D(I) = Z(I,I)                                              
          Z(I,I) = 1.0D0                                              
          IF(L.LT.1) GO TO 500                                     
C                                                                       
          DO 400 J = 1, L                                             
             Z(I,J) = 0.0D0                                           
             Z(J,I) = 0.0D0                                           
  400     CONTINUE                                                    
C                                                                       
  500  CONTINUE                                                      
C                                                                       
      RETURN                                                            
C*************** LAST CARD OF TRED2 **********                          
      END SUBROUTINE TRED2   
C                         
C--------------------------------------------------------------------
C                                              
      SUBROUTINE TQLRAT(N,D,E2,IERR)                                    
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION D(N),E2(N)                                              
      REAL*8 MACHEP                                                     
C                                                                       
C.....THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,   
C.....ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.                
C                                                                       
C.....THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC              
C.....TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.                     
C                                                                       
C.....ON INPUT-                                                         
C                                                                       
C.....N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C.....D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          
C                                                                       
C.....E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE     
C.....INPUT MATRIX IN ITS LAST N-1 POSITIONS. E2(1) IS ARBITRARY. 
C                                                                       
C.....ON OUTPUT-                                                       
C                                                                       
C.....D CONTAINS THE EIGENVALUES IN ASCENDING ORDER. IF AN          
C.....ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND          
C.....ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE            
C.....THE SMALLEST EIGENVALUES,                                    
C                                                                       
C.....E2 HAS BEEN DESTROYED,                                         
C                                                                       
C.....IERR IS SET TO                                                 
C.....     ZERO       FOR NORMAL RETURN,                                
C.....     J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C.....                DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C.....QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C.....APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C********************************************************************
C                                                                       
C********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C********** THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                                                                       
       MACHEP = 2.D0**(-26)                                          
C                                                                       
       IERR = 0                                                    
       IF(N.EQ.1) GO TO 1001                                     
C                                                                       
       DO 100 I = 2, N                                              
  100  E2(I-1) = E2(I)                                               
C                                                                    
       F = 0.0D0                                                       
       B = 0.0D0                                                    
       C = 0.0D0                                                   
       E2(N) = 0.0D0                                                 
C                                                                       
       DO 290 L = 1, N                                               
          J = 0                                                       
          H = MACHEP * (DABS(D(L)) + DSQRT(E2(L)))                  
          IF(B.GT.H) GO TO 105                                   
          B = H                                                         
          C = B * B                                                   
C*************** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ********** 
  105     DO 110 M = L, N                                           
             IF(E2(M).LE.C) GO TO 120                            
C*************** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT              
C                THROUGH THE BOTTOM OF THE LOOP **********              
  110     CONTINUE                                                  
          WRITE(6,*) '**** FATAL ERROR IN TQLRAT ****'                  
          WRITE(6,*) '**** FALLEN THROUGH BOTTOM OF LOOP 110 ***'       
          STOP                                                          
C                                                                       
  120     IF(M.EQ.L) GO TO 210                                     
  130     IF(J.EQ.30) GO TO 1000                                  
          J = J + 1                                                  
C*************** FORM SHIFT **********                                  
          L1 = L + 1                                                 
          S = DSQRT(E2(L))                                           
          G = D(L)                                                    
          P = (D(L1) - G) / (2.0D0 * S)                              
          R = DSQRT(P*P+1.0D0)                                      
          D(L) = S / (P + DSIGN(R,P))                                
          H = G - D(L)                                               
C                                                                   
          DO 140 I = L1, N                                          
  140     D(I) = D(I) - H                                           
C                                                                       
          F = F + H                                                  
C*************** RATIONAL QL TRANSFORMATION **********                  
          G = D(M)                                                 
          IF(G.EQ.0.0D0) G = B                                 
          H = G                                                    
          S = 0.0D0                                                 
          MML = M - L                                                
C*************** FOR I=M-1 STEP -1 UNTIL L DO -- **********             
          DO 200 II = 1, MML                                         
             I = M - II                                              
             P = G * H                                              
             R = P + E2(I)                                         
             E2(I+1) = S * R                                       
             S = E2(I) / R                                         
             D(I+1) = H + S * (H + D(I))                            
             G = D(I) - E2(I) / G                                 
             IF(G.EQ.0.0D0) G = B                                 
             H = G * P / R                                           
  200     CONTINUE                                                    
C                                                                       
          E2(L) = S * G                                             
          D(L) = H                                                    
C*************** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ********** 
          IF(H.EQ.0.0D0) GO TO 210                                 
          IF(DABS(E2(L)).LE.DABS(C/H)) GO TO 210                 
          E2(L) = H * E2(L)                                         
          IF(E2(L).NE.0.0D0) GO TO 130                            
  210     P = D(L) + F                                               
C*************** ORDER EIGENVALUES **********                           
          IF(L.EQ.1) GO TO 250                                     
C*************** FOR I=L STEP -1 UNTIL 2 DO -- **********               
          DO 230 II = 2, L                                            
             I = L + 2 - II                                          
             IF(P.GE.D(I-1)) GO TO 270                            
             D(I) = D(I-1)                                           
  230     CONTINUE                                                   
C                                                                       
  250     I = 1                                                      
  270     D(I) = P                                                   
  290  CONTINUE                                                      
C                                                                       
       GO TO 1001                                                    
C*************** SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS **********              
 1000  IERR = L                                                       
 1001  RETURN                                                         
C*************** LAST CARD OF TQLRAT **********                         
      END SUBROUTINE TQLRAT             
C                                                                       
C--------------------------------------------------------------------
C                                             
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION D(N),E(N),Z(NM,N)                                       
      REAL*8 MACHEP                                                     
C                                                                       
C.....THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     
C.....NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     
C.....WILKINSON.                                                        
C.....HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   
C                                                                      
C.....THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            
C.....OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               
C.....THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              
C.....BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  
C.....FULL MATRIX TO TRIDIAGONAL FORM.                                  
C                                                                       
C.....ON INPUT-                                                         
C                                                                       
C.....NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C.....ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C.....DIMENSION STATEMENT,                                         
C                                                                       
C.....N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C.....D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          
C                                                                       
C.....E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        
C.....IN ITS LAST N-1 POSITIONS. E(1) IS ARBITRARY,               
C                                                                       
C.....Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           
C.....REDUCTION BY  TRED2, IF PERFORMED. IF THE EIGENVECTORS      
C.....OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        
C.....THE IDENTITY MATRIX.                                         
C                                                                       
C.....ON OUTPUT-                                                       
C                                                                       
C.....D CONTAINS THE EIGENVALUES IN ASCENDING ORDER. IF AN          
C.....ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          
C.....UNORDERED FOR INDICES 1,2,...,IERR-1,                        
C                                                                       
C.....E HAS BEEN DESTROYED,                                          
C                                                                       
C.....Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           
C.....TRIDIAGONAL (OR FULL) MATRIX. IF AN ERROR EXIT IS MADE,     
C.....Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       
C.....EIGENVALUES,                                                 
C                                                                       
C.....IERR IS SET TO                                                 
C.....     ZERO       FOR NORMAL RETURN,                                
C.....     J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C.....                DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C.....QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C.....APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C********************************************************************
C                                                                       
C********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C           THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                                                                       
      MACHEP = 2.D0**(-26)                                              
C                                                                       
      IERR = 0                                                          
       IF(N.EQ.1) GO TO 1001                                      
C                                                                       
       DO 100 I = 2, N                                              
  100  E(I-1) = E(I)                                                 
C                                                                       
      F = 0.0D0                                                         
      B = 0.0D0                                                         
      E(N) = 0.0D0                                                      
C                                                                       
       DO 240 L = 1, N                                              
          J = 0                                                    
          H = MACHEP * (DABS(D(L)) + DABS(E(L)))                    
          IF(B.LT.H) B = H                                        
C********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********         
          DO 110 M = L, N                                            
             IF(DABS(E(M)).LE.B) GO TO 120                        
C********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               
C           THROUGH THE BOTTOM OF THE LOOP **********              
  110     CONTINUE                                                   
C                                                                       
  120     IF(M.EQ.L) GO TO 220                                    
  130     IF(J.EQ.30) GO TO 1000                                   
          J = J + 1                                                  
C********** FORM SHIFT **********                                  
          L1 = L + 1                                                 
          G = D(L)                                                   
          P = (D(L1) - G) / (2.0D0 * E(L))                            
          R = DSQRT(P*P+1.0D0)                                        
          D(L) = E(L) / (P + DSIGN(R,P))                              
          H = G - D(L)                                                
C                                                                       
          DO 140 I = L1, N                                            
  140     D(I) = D(I) - H                                              
C                                                                       
          F = F + H                                                   
C********** QL TRANSFORMATION **********                           
          P = D(M)                                                   
          C = 1.0D0                                                  
          S = 0.0D0                                                   
          MML = M - L                                                 
C********** FOR I=M-1 STEP -1 UNTIL L DO -- **********             
          DO 200 II = 1, MML                                         
             I = M - II                                               
             G = C * E(I)                                             
             H = C * P                                               
             IF(DABS(P).LT.DABS(E(I))) GO TO 150                  
             C = E(I) / P                                           
             R = DSQRT(C*C+1.0D0)                                    
             E(I+1) = S * P * R                                       
             S = C / R                                                
             C = 1.0D0 / R                                           
             GO TO 160                                               
  150        C = P / E(I)                                            
             R = DSQRT(C*C+1.0D0)                                    
             E(I+1) = S * E(I) * R                                   
             S = 1.0D0 / R                                           
             C = C * S                                                
  160        P = C * D(I) - S * G                                     
             D(I+1) = H + S * (C * G + S * D(I))                      
C********** FORM VECTOR **********                                 
             DO 180 K = 1, N                                         
                H = Z(K,I+1)                                         
                Z(K,I+1) = S * Z(K,I) + C * H                        
                Z(K,I) = C * Z(K,I) - S * H                          
  180        CONTINUE                                                 
C                                                                       
  200     CONTINUE                                                   
C                                                                       
          E(L) = S * P                                              
          D(L) = C * P                                              
          IF(DABS(E(L)).GT.B) GO TO 130                          
  220     D(L) = D(L) + F                                            
  240  CONTINUE                                                      
C********** ORDER EIGENVALUES AND EIGENVECTORS **********          
       DO 300 II = 2, N                                             
          I = II - 1                                                
          K = I                                                      
          P = D(I)                                                   
C                                                                       
          DO 260 J = II, N                                           
             IF(D(J).GE.P) GO TO 260                            
             K = J                                                   
             P = D(J)                                                 
  260     CONTINUE                                                    
C                                                                       
          IF(K.EQ.I) GO TO 300                                    
          D(K) = D(I)                                                 
          D(I) = P                                                    
C                                                                    
          DO 280 J = 1, N                                            
             P = Z(J,I)                                               
             Z(J,I) = Z(J,K)                                         
             Z(J,K) = P                                               
  280     CONTINUE                                                    
C                                                                     
  300  CONTINUE                                                      
C                                                                       
       GO TO 1001                                                     
C********** SET ERROR -- NO CONVERGENCE TO AN                      
C           EIGENVALUE AFTER 30 ITERATIONS **********              
 1000  IERR = L                                                      
 1001  RETURN                                                        
C********** LAST CARD OF TQL2 **********                           
      END SUBROUTINE TQL2   
C                                                                       
C====================================================================
C                                                                       
      DOUBLE PRECISION FUNCTION BESSJ0(X)                               
C.....G.G. BALINT-KURTI -- CRUDE BUT ACCURRATE VERSION.                 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DATA PI/3.141592653589793D0/                                      
       IF(DABS(X).LT.20.0D0) THEN                                      
          Z=-X**2/4.0D0                                                
          SUM=1.0D0                                                   
          T=1.0D0                                                    
          DO 046 I=1,1000                                             
             T=T*Z/(DFLOAT(I*I))                                      
             SUM=SUM+T                                              
             IF(DABS(T).LT.1.0D-20) go to 10                         
046       CONTINUE                                                    
          WRITE(6,*)' BESSEL FUNCTION J0 HAS NOT CONVERGED '          
          STOP                                                        
10        BESSJ0=SUM                                                 
       ELSE                                                           
          AX=DABS(X)                                                 
          Z=20.0D0/AX                                                
          Y=Z**2                                                   
          XX=AX - PI*0.25D0                                         
          BESSJ0=DSQRT(2.0D0/(AX*PI))*                               
     &         (DCOS(XX)*PFUNC0(AX) - DSIN(XX)*QFUNC0(AX))              
       END IF                                                          
      RETURN                                                            
      END FUNCTION BESSJ0                                      
C                        
C====================================================================
C                                                                             
      DOUBLE PRECISION FUNCTION PFUNC0(X)                               
C.....G.G. BALINT-KURTI .                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      Z = 0.25D0/(X**2)                                                 
      SUM=1.0D0                                                         
      T=1.0D0                                                           
       DO 047 I=1,1000                                                
          TOLD=T                                                      
          T = -T*DFLOAT(((4*I-1)*(4*I-3))**2)*0.0625D0*Z/             
     &        (DFLOAT(2*I*(2*I-1)))                                  
          SUM=SUM+T                                                  
          IF(DABS(T).LT.1.0D-16) go to 10                             
             IF(DABS(T).GT.DABS(TOLD)) THEN                           
             WRITE(6,*) 'ASYMPTOTIC SERIES PFUNC0 NOT CONVERGING'      
             STOP                                                   
          END IF                                                     
047    CONTINUE                                                      
      WRITE(6,*) 'FUNCTION  PFUNC0 HAS NOT CONVERGED'                  
      STOP                                                              
10    PFUNC0=SUM                                                        
      RETURN                                                            
      END FUNCTION PFUNC0                      
C                                                                       
C====================================================================
C                                                                       
      DOUBLE PRECISION FUNCTION QFUNC0(X)                               
C.....G.G. BALINT-KURTI .                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      Z = 0.25D0/(X**2)                                                 
      T= -0.125D0/X                                                     
      SUM=T                                                             
       DO 048 I=1,1000                                               
          TOLD=T                                                     
          T = -T*DFLOAT(((4*I+1)*(4*I-1))**2)*0.0625D0*Z/             
     &        (DFLOAT(2*I*(2*I+1)))                                  
          SUM=SUM+T                                                  
          IF(DABS(T).LT.1.0D-16) go to 10                             
             IF(DABS(T).GT.DABS(TOLD)) THEN                           
             WRITE(6,*) 'ASYMPTOTIC SERIES QFUNC0 NOT CONVERGING'      
             STOP                                                     
          END IF                                                     
048    CONTINUE                                                      
      WRITE(6,*) 'FUNCTION  QFUNC0 HAS NOT CONVERGED'                  
      STOP                                                              
10    QFUNC0=SUM                                                        
      RETURN                                                            
      END FUNCTION QFUNC0
C                                    
C====================================================================
C                                                                       
      DOUBLE PRECISION FUNCTION BESSJ1(X)                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DATA PI/3.141592653589793D0/                                      
       IF(DABS(X).LT.20.0D0) THEN                                     
          Z=-X**2/4.0D0                                              
          BESSJ1=X/2.0D0                                            
          SUM=1.0D0                                                   
          T=1.0D0                                                    
          DO 049 I=1,1000                                             
             T=T*Z/(DFLOAT(I*(I+1)))                                   
             SUM=SUM+T                                                
             IF(DABS(T).LT.1.0D-20) go to 10                           
049       CONTINUE                                                    
          WRITE(6,*) 'BESSEL FUNCTION J1 HAS NOT CONVERGED'           
          STOP                                                         
10        BESSJ1=BESSJ1*SUM                                            
       ELSE                                                           
          AX=DABS(X)                                                  
          Z=20.0D0/AX                                                 
          Y=Z**2                                                      
          XX=AX - PI*0.75D0                                           
          BESSJ1=DSQRT(2.0D0/(AX*PI))*                               
     &           (DCOS(XX)*PFUNC1(AX) - DSIN(XX)*QFUNC1(AX))          
       END IF                                                          
      RETURN                                                            
      END FUNCTION BESSJ1                                
C                                                                       
C====================================================================
C                                               
      DOUBLE PRECISION FUNCTION PFUNC1(X)                               
C.....G.G. BALINT-KURTI .                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      Z = 0.25D0/(X**2)                                                 
      SUM=1.0D0                                                         
      T=1.0D0                                                           
       DO 050 I=1,1000                                                
          TOLD=T                                                      
          T = -T*DFLOAT((4*I+1)*(4*I-1)*(4*I-3)*(4*I-5))*0.0625D0*Z/  
     &        (DFLOAT(2*I*(2*I-1)))                                   
          SUM=SUM+T                                                    
          IF(DABS(T).LT.1.0D-16) go to 10                              
             IF(DABS(T).GT.DABS(TOLD)) THEN                           
                WRITE(6,*) 'ASYMPTOTIC SERIES PFUNC1 NOT CONVERGING'   
                STOP                                                    
             END IF                                                    
050    CONTINUE                                                      
      WRITE(6,*) 'FUNCTION PFUNC1 HAS NOT CONVERGED'                   
      STOP                                                              
10        PFUNC1=SUM                                                  
      RETURN                                                            
      END FUNCTION PFUNC1                   
C                        
C====================================================================
C                                               
      DOUBLE PRECISION FUNCTION QFUNC1(X)                               
C.....G.G. BALINT-KURTI .                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      Z = 0.25D0/(X**2)                                                 
      T= 0.375D0/X                                                      
      SUM=T                                                             
       DO 051 I=1,1000                                                 
          TOLD=T                                                       
          T = -T*DFLOAT((4*I+3)*(4*I+1)*(4*I-1)*(4*I-3))*0.0625D0*Z/   
     &        (DFLOAT(2*I*(2*I+1)))                                     
          SUM=SUM+T                                                    
          IF(DABS(T).LT.1.0D-16) go to 10                              
             IF(DABS(T).GT.DABS(TOLD)) THEN                            
                WRITE(6,*) 'ASYMPTOTIC SERIES QFUNC1 NOT CONVERGING'  
                STOP                                            
             END IF                                                   
051    CONTINUE                                                     
      WRITE(6,*) 'FUNCTION QFUNC1 HAS NOT CONVERGED'                   
      STOP                                                          
10        QFUNC1=SUM                                                  
      RETURN                                                            
      END FUNCTION QFUNC1                       
C                                                                       
C====================================================================
C                                               
      DOUBLE PRECISION FUNCTION BESSJ(BJ0,BJ1,N,X)                      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      PARAMETER (IACC=200,BIGNO=1.D20,BIGNI=1.D-20)                     
       IF(N.LT.2) PAUSE 'BAD ARGUMENT N IN BESSJ'                     
      TOX=2.0D0/X                                                       
       IF(X.GT.DFLOAT(N)) THEN                                        
          BJM=BJ0                                                     
          BJ=BJ1                                                      
CC          BJM=BJ0                                                         
CC          BJ=BJ1                                                          
          DO 11 J=1,N-1                                                
             BJP=J*TOX*BJ-BJM                                         
             BJM=BJ                                                    
             BJ=BJP                                                  
11        CONTINUE                                                    
          BESSJ=BJ                                                   
       ELSE                                                          
          M=2*((N+INT(SQRT(FLOAT(IACC*N))))/2)                      
          BESSJ=0.0D0                                                
          JSUM=0                                                     
          SUM=0.0D0                                                   
          BJP=0.0D0                                                   
          BJ=1.0D0                                                   
          DO 12 J=M,1,-1                                              
             BJM=J*TOX*BJ-BJP                                         
             BJP=BJ                                                   
             BJ=BJM                                                  
             IF(DABS(BJ).GT.BIGNO) THEN                               
                BJ=BJ*BIGNI                                           
                BJP=BJP*BIGNI                                        
                BESSJ=BESSJ*BIGNI                                    
                SUM=SUM*BIGNI                                          
             END IF                                                    
             IF(JSUM.NE.0) SUM=SUM+BJ                                
             JSUM=1-JSUM                                              
             IF(J.EQ.N) BESSJ=BJP                                     
12        CONTINUE                                                   
          SUM=2.0D0*SUM-BJ                                          
          BESSJ=BESSJ/SUM                                            
       END IF                                                          
      RETURN                                                            
      END FUNCTION BESSJ                                    
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      SUBROUTINE FOUR1(DATA,NN,ISGN)                                    
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                
      DIMENSION DATA(2*NN)                                              
      DATA PI/3.141592653589793D0/                                      
      TPI=2.0D0*PI                                                      
      N=2*NN                                                            
      C=1.0D0/DSQRT(DFLOAT(NN))                                         
      J=1                                                               
       DO 011 I=1,N,2                                                
          IF(J.GT.I) THEN                                             
             TEMPR=DATA(J)                                            
             TEMPI=DATA(J+1)                                        
             DATA(J)=DATA(I)                                          
             DATA(J+1)=DATA(I+1)                                      
             DATA(I)=TEMPR                                           
             DATA(I+1)=TEMPI                                         
          END IF                                                       
          M=N/2                                                        
001       IF((M.GE.2).AND.(J.GT.M)) THEN                              
             J=J-M                                                     
             M=M/2                                                  
             GO TO 001                                               
          END IF                                                      
          J=J+M                                                      
011    CONTINUE                                                       
      MMAX=2                                                            
002    IF(N.GT.MMAX) THEN                                            
          ISTEP=2*MMAX                                               
          THETA=TPI/(ISGN*MMAX)                                       
          WPR=-2.0D0*SIN(0.5D0*THETA)**2                             
          WPI=SIN(THETA)                                             
          WR=1.                                                       
          WI=0.                                                      
          DO 013 M=1,MMAX,2                                           
             DO 012 I=M,N,ISTEP                                        
                J=I+MMAX                                               
                TEMPR= WR*DATA(J)- WI*DATA(J+1)                      
                TEMPI= WR*DATA(J+1)+ WI*DATA(J)                    
                DATA(J)=DATA(I)-TEMPR                              
                DATA(J+1)=DATA(I+1)-TEMPI                           
                DATA(I)=DATA(I)+TEMPR                               
                DATA(I+1)=DATA(I+1)+TEMPI                            
012          CONTINUE                                                
             WTEMP=WR                                             
             WR=WR*WPR-WI*WPI+WR                              
             WI=WI*WPR+WTEMP*WPI+WI                                
013       CONTINUE                                               
          MMAX=ISTEP                                             
          GO TO 002                                                  
       END IF                                                        
       DO 014 I=1,N                                                   
          IR=2*I-1                                                    
          II=2*1                                                       
          DATA(I)=DATA(I)*C                                            
014    CONTINUE                                                     
      RETURN                                                            
      END SUBROUTINE FOUR1                                         


