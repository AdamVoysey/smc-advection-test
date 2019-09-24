!!
!! Adapted for multiple cell 2D advection tests using UNO schemes.
!!                    J G Li   26 Jul 2007
!! Adapted for global multiple cell advection tests using UNO2 scheme.
!!                    J G Li   22 Nov 2007
!! Modified for global SMC grid extended to cover the N Pole.
!!                    J G Li   26 Nov 2008
!! Adapted for 2-part, 3-step, 3-resolution G6kSMC grid spectral transport.
!!                    J G Li    5 Mar 2010
!! Changed for global only SMC625 grid spectral transport.
!!                    J G Li   12 Dec 2011
!! Automatic setting of multi-resolution loops  with MRL and MRFct.
!!                    J G Li   28 Feb 2014
!! Add hybrid MPI-OpenMP lines for hybrid test on Cray XC40.
!!                    J G Li   12 Dec 2017
!! Separate W3PSMCMD for easy comparison with WW3 code. 
!!                    J G Li   10 Jul 2019
!! Fix bugs, update with OpenMP Atomic lines and tidy up parameters. 
!!                    J G Li   12 Jul 2019
!!


      PROGRAM SMCPrDRG 
       USE Constants
       USE W3PSMCMD
       IMPLICIT NONE

       REAL:: CNST, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6, CNST8
       REAL (Kind=8) :: TM1, TM2, GTM1, GTM2

       REAL :: t1,t2
!  !$ACC Routine(w3krtn) SEQ
! Initialize MPI with threading
       call MPI_INIT_THREAD(required, provided, ierr)
       MPI_COM = MPI_COMM_WORLD

       call MPI_COMM_SIZE(MPI_COM, nprocs, ierr)
       call mpi_comm_rank(MPI_COM, myrank, ierr)

       CALL get_environment_variable("SCRDIR", CelPath)

!  Setup allocatable arrays for broadcast among ranks
       ALLOCATE( INTALLOC(15), INTALLO2(0:MRL, 3), STAT=malloc )

!  Read Global and Arctic part Multiple-Cell info on rank 0 
       if (myrank .eq. 0) then
           CALL READCELL

           INTALLOC( 1)=NC
           INTALLOC( 2)=NU
           INTALLOC( 3)=NV
           INTALLOC( 4)=NGLo
           INTALLOC( 5)=NGLA
           INTALLOC( 6)=NGLB
           INTALLOC( 7)=NArc
           INTALLOC( 8)=NArA
           INTALLOC( 9)=NArB
           INTALLOC(10)=NUGL
           INTALLOC(11)=NUAr
           INTALLOC(12)=NVGL
           INTALLOC(13)=NVAr
           INTALLOC(14)=NNor
           INTALLOC(15)=NSou

           INTALLO2(0:MRL,1)=NRLCel(0:MRL)
           INTALLO2(0:MRL,2)=NRLUFc(0:MRL)
           INTALLO2(0:MRL,3)=NRLVFc(0:MRL)
       endif

!  Broadcast cell and face parameters to all ranks
       CALL MPI_Bcast(INTALLOC, 15, MPI_INTEGER, 0, MPI_COM, ierr)
       m = (MRL+1)*3
       CALL MPI_Bcast(INTALLO2,  m, MPI_INTEGER, 0, MPI_COM, ierr)

       CALL MPI_BARRIER(MPI_COM, ierr)

!  Assign parameters back to original names
       IF( myrank .NE. 0 ) THEN
           NC=  INTALLOC( 1)
           NU=  INTALLOC( 2)
           NV=  INTALLOC( 3)
           NGLo=INTALLOC( 4)
           NGLA=INTALLOC( 5)
           NGLB=INTALLOC( 6)
           NArc=INTALLOC( 7)
           NArA=INTALLOC( 8)
           NArB=INTALLOC( 9)
           NUGL=INTALLOC(10)
           NUAr=INTALLOC(11)
           NVGL=INTALLOC(12)
           NVAr=INTALLOC(13)
           NNor=INTALLOC(14)
           NSou=INTALLOC(15)

           NRLCel(0:MRL)=INTALLO2(0:MRL,1)
           NRLUFc(0:MRL)=INTALLO2(0:MRL,2)
           NRLVFc(0:MRL)=INTALLO2(0:MRL,3)
       ENDIF

! Broadcast cell and face arrays among all ranks.
       n=5*(NCL+10)
       CALL MPI_Bcast(ICE(1,-9), n, MPI_INTEGER, 0, MPI_COM, ierr)
       m=7*NFC
       CALL MPI_Bcast(ISD(1, 1), m, MPI_INTEGER, 0, MPI_COM, ierr)
       k=8*NFC
       CALL MPI_Bcast(JSD(1, 1), k, MPI_INTEGER, 0, MPI_COM, ierr)

       CALL MPI_BARRIER(MPI_COM, ierr)

       DEALLOCATE( INTALLOC, INTALLO2, STAT=malloc )

!!  Evaluate cell centre and V-face latitude cosine
          CNST=0.5*DLAT
!$OMP Parallel DO Private(CNST2)
       DO n=1, NC-NPol
          CNST2=Real( 2*ICE(2,n)+ICE(4,n) )*CNST + ZrLat
          CLats(n)=COS(  CNST2*D2RAD )
          CTHG0S(n)= - TAN( CNST2*D2RAD ) / REARTH
       ENDDO
!$OMP END Parallel DO

!$OMP Parallel DO Private(CNST3)
       DO j=1, NV
          CNST3=Real( JSD(2,j) )*DLat + ZrLat
          CLatF(j)=COS( CNST3*D2RAD )
       ENDDO
!$OMP END Parallel DO
     
!!  Define cell area for cell update
       DX0=DLON*D2RAD
       DY=DLAT*D2RAD
       DYR=1.0/DY

       IF( NPol .GT. 0 ) THEN 
!!  Polar cells are round cell of radius R=PCRDY*DY*ICE(4,NNor)
!!  The RCELA(NC) represent the net V-flux factor for the polar cell.
!!  Net V-flux are all divided by the CLats factor for cell update
!!  but it is not needed for polar cell.  To avoid zero-dividing 
!!  polar cell CLats values are evaluated at cell edge rather than centre 
!!  and are set equal for south polar cell if any. 
!!  Add polar cell GCT factor CTHG0S though it will be overwritten in sub ArctAngd.
          CNST2=Real( ICE(2,NNor) )*DLat + ZrLat
          CLats(NNor)=COS(  CNST2*D2RAD )
          CNST1=Real(ICE(4,NNor))*DY*PCRDY
          CNST2=Pie*CNST1*CNST1
          RCELA(NNor)=DX0*DY*CLats(NNor)/CNST2 
          CTHG0S(NNor)= 0.0 
       ENDIF
!!  Duplicate north polar cell area for south polar cell if NPol = 2
       IF( NPol .GT. 1 ) THEN
          CLats(NSou)= CLats(NNor)
          RCELA(NSou)= RCELA(NNor)
          CTHG0S(NSou)= 0.0 
       ENDIF

!!   Evaluate other cell dx length in rad and cell area in rad**2
!!   except for the polar cell(s).
!$OMP Parallel DO
       DO L=-9, NC-NPol
          RCELA(L)=1.0/Real( ICE(3,L)*ICE(4,L) )
       ENDDO
!$OMP END Parallel DO

!!  Directional bins in rad, shifted by half-width from axises
        DThta=2.0*Pie/FLOAT(NDIR)
      DO K=1, NDIR
         Theta(K)=(FLOAT(K) - 0.5)*DThta
         ECOS(K)=COS(Theta(K))
         ESIN(K)=SIN(Theta(K))
         EC2(K) = ECOS(K)**2
         ES2(K) = ESIN(K)**2
         ESC(K) = ECOS(K)*ESIN(K)
      ENDDO

!!  Convert diffusivity into radian^2 s-1 and multiply by 2*DT
!!  The 2.0 factor is added to cancel the 0.5 factor in the gradient.
        AKHDT2 = 2.0*AKH*DTG/(REARTH*REARTH)

!!  Set up frequency bands and work out wave number and group speed
        CALL SgmCgWnK

!!  Assign spectral components to MPI ranks for balanced load.
        CALL Spctinit

!   Whole array assignment for i=-8,NCL
        C=0.0

!!  Multiple factor for multi-resolution loops
        MRFct=2**(MRL - 1)
!       WRITE(6,*) ' Multi-Resolution level and factor=', MRL, MRFct

!!  Some duplicated variables for WW3 model
        SX = DLon*MRFct
        SY = DLat*MRFct
        NSEA = NC
      
!!  Workout number of spectra or sea points per MPI rank
        npsea = MAX(1, NSEA/nprocs)
        npspc = MAX(1, NSpc/nprocs)

!  Initialise AngU/V/CD to be zero.
        AngU=0.0
        AngV=0.0
        AngCD=0.0
        ELaCD=0.0

!  Generate flux face and cell rotation Angle for Arctic cells
        IF( Arctic )  CALL ArctAngd

!  Initialise wave spectra and bin velocities,
!  including GCT and refraction Courant numbers.
        CALL SPECUUVV

!  Read restart time step or set default = 0
        NS = 0
!       WRITE(6,*) " Restart time step NS set to be ", NS

       if (myrank .eq. 0) then
!  Save a copy of initial values at writeup time step.
        NT=NS
        write(FL9NM(3:7), FMT='(i5)' )  10000+NT
        OPEN(UNIT=26, FILE=FL9NM, STATUS='NEW',IOSTAT=nn)
        IF(nn /= 0) PRINT*,' File FL9NM was not opened! '
        WRITE(UNIT=6,FMT='(2x,"NT= ",i6,3x,A9)') NT,FL9NM

!!   Integration of spectral energy for initial field.
           ii=0
        DO i=1, NC
           CTT=0.0
           DO j=1, NFrq
           DO k=1, NDIR
              CTT = CTT + WSpc(k,j,i)*DSIP(j)
           ENDDO
           ENDDO
           C(i)=CTT*DThta
           D(i) = SIGN( SQRT( ABS(C(i)) ), C(i) )

!!   Filter very small C(n) value so it is greater than E-36
           IF( Abs(D(i)) .LT. 1.0E-36 ) THEN
               D(i)=SIGN(1.0E-36, C(i))
!Li  Resetting NaNQ VQ to zero if any.   JGLi14Nov2017
           ELSE IF ( .NOT. (C(i) .EQ. C(i)) ) THEN
               D(i) = 0.0
!              WRITE(6,*) "NaN found at i=", i, C(i)
           ENDIF

        ENDDO

!    Central basic cells only, all at end of cell array
        WRITE(UNIT=26, FMT='(2x,2i8)' )  NT, NC
        WRITE(UNIT=26, FMT=7113)  (D(n), n=1,NC)

        CLOSE(26)

!!  End of rank 0 output of initial wave field
       endif

!!  Calculate bathymetry gradient DHDX/DY and refraction coefficient.
        CALL SMCDHXY 

!!  Read current velocity and calculate gradients if available.
!   Moved inside time loop to read for each hour.  JGLi23Nov2017
!       IF ( FLCUR ) THEN
!          CALL READCURNT(NS)
!          CALL SMCDCXY
!       ENDIF

!    Define cell and face sub-step counts.  JGLi18Jan2012
         NRLCel(0)=0
         NRLUFc(0)=0
         NRLVFc(0)=0
         DO i=1,MRL-1
         NRLCel(i)=NRLCel(i)+NRLCel(i-1)
         NRLUFc(i)=NRLUFc(i)+NRLUFc(i-1)
         NRLVFc(i)=NRLVFc(i)+NRLVFc(i-1)
         ENDDO
         NRLCel(MRL)=NC
         NRLUFc(MRL)=NU
         NRLVFc(MRL)=NV
         
!$     TM1= OMP_GET_WTIME()
!$OMP Parallel  Private(i, j)
!$    i=omp_get_thread_num()
!$    j=omp_get_num_procs()
!$    m=omp_get_max_threads()
!$    n=omp_get_num_threads()
!$    IF(i .EQ. 0 .AND. myrank .EQ. 0) THEN 
!$       WRITE( 6,*) "Num_Threads to be used =", n, m, j 
!$    ENDIF
!$OMP END Parallel
         nthreads = n

       if (myrank .eq. 0) then
!     Open files to store writups
       OPEN(UNIT=16,FILE='CMesgs.txt',STATUS='UNKNOWN',IOSTAT=nn, &
        &           ACTION='WRITE')
       IF(nn /= 0) PRINT*,' File CMesgs.txt was not opened! '

!     Header messages and configuration information 
       WRITE(UNIT=16,FMT='(1x/   &
        &  "  Global Multiple-Cell 2-D Spherical Advection Model" /    &
        &  "         SMC Version 2.0   J G  Li  Mar 2010  " /)' )

       CALL DATE_AND_TIME(CDate, CTime)
       WRITE(UNIT=16,FMT="(1x,' Run time date ',A10,2x,A10)") CTime, CDate

       WRITE(UNIT=16,FMT='(1x," Size-1 Units DLON DLAT = ",2f14.10)')  DLON, DLAT
       WRITE(UNIT=16,FMT='(1x," Equatorial PoLat PoLon = ",2f8.2)')  PoLat, PoLon
       WRITE(UNIT=16,FMT='(1x," Standard grid ZrLatLon = ",2f9.5)')  ZrLat, ZrLon
       WRITE(UNIT=16,FMT='(1x," Horizontal diffusivity = ",f8.1)' )  AKH
       WRITE(UNIT=16,FMT='(1x," Global time step DTG s = ",f8.1)' )  DTG
       WRITE(UNIT=16,FMT='(1x," CFL limited step DTCFL = ",f8.1)' )  DTCFL 
       WRITE(UNIT=16,FMT='(1x," Maximum grid speed s-1 = ",ES12.3)') UMX
       WRITE(UNIT=16,FMT='(1x," Maximum Courant number = ",f8.3)' )  CMX
       WRITE(UNIT=16,FMT='(1x," Max GCT Courant number = ",f8.3)' )  CGCMX
       WRITE(UNIT=16,FMT='(1x," Max Rfr Courant number = ",f8.3)' )  CRFMX
       WRITE(UNIT=16,FMT='(1x," Initial integrated SWH = ",f8.3)' )  SWH0
       WRITE(UNIT=16,FMT='(1x," First Freqncy Frqc0 Hz = ",f8.4)' )  Frqc0
       WRITE(UNIT=16,FMT='(1x," Init Wave Peak Freqncy = ",f8.4)' )  PkFrq
       WRITE(UNIT=16,FMT='(1x," Nos. of Dirs and Frqcy = ",2i8)' ) NDir, NFrq
       WRITE(UNIT=16,FMT='(1x," Start & Total timestes = ",2i8)' ) NS, NTS
       WRITE(UNIT=16,FMT='(1x," Writeup timestep every = ",i8)' )  NWP
       WRITE(UNIT=16,FMT='(1x," Multi-reso levl factor = ",2i8)')  MRL, MRFct
       WRITE(UNIT=16,FMT='(1x," Horizontal cell number = ",2i8)')  NC, NCL
       WRITE(UNIT=16,FMT='(1x," Globl/Arctic cell No.s = ",2i8)')  NGLo, NArc
       WRITE(UNIT=16,FMT='(1x," Globl/Arctic bndy No.s = ",2i8)')  NGLB, NArB
       WRITE(UNIT=16,FMT='(1x," Total number of U-face = ",2i8)')  NU, NFC
       WRITE(UNIT=16,FMT='(1x," Globl/Arctic U-face No = ",2i8)')  NUGL, NUAr
       WRITE(UNIT=16,FMT='(1x," Total number of V-face = ",i8)' )  NV
       WRITE(UNIT=16,FMT='(1x," Globl/Arctic V-face No = ",2i8)')  NVGL, NVAr
       WRITE(UNIT=16,FMT='("Sub-step cell count NRLCel(0:MRL)=",6i8)') NRLCel
       WRITE(UNIT=16,FMT='("Sub-step Ufce count NRLUFc(0:MRL)=",6i8)') NRLUFc
       WRITE(UNIT=16,FMT='("Sub-step Vfce count NRLVFc(0:MRL)=",6i8)') NRLVFc
       WRITE(UNIT=16,FMT='(1x," Total number MPI ranks = ",i8)' )  nprocs
       WRITE(UNIT=16,FMT='(1x," No. PRk OpenMP threads = ",i8)' )  nthreads
       WRITE(UNIT=16,FMT='(1x," Spctr/Sea pts per rank = ",2i8)')  npspc, npsea

 3912 FORMAT(1x,i4,3F9.1,ES12.3)

       WRITE(16,FMT='(/1x," Frequency 0:NFrq+1 in Hz")' ) 
       WRITE(16,FMT='((10F8.5))')  (0.5*Sig(k)/Pie,k=0,NFrq+1)
       WRITE(16,FMT='(/1x," Initial PM frequency factor")' ) 
       WRITE(16,FMT='((10F8.5))')  (SpecPM(k),k=0,NFrq+1)
       WRITE(16,FMT='(/1x," First and last ICE values ")' ) 
       WRITE(16,FMT='(1x,6i8)')  1,(ICE(i, 1),i=1,5)
       WRITE(16,FMT='(1x,6i8)') NC,(ICE(i,NC),i=1,5)
       WRITE(16,FMT='(/1x," First and last ISD values ")' ) 
       WRITE(16,FMT='(1x,8i8)')  1,(ISD(i, 1),i=1,7)
       WRITE(16,FMT='(1x,8i8)') NU,(ISD(i,NU),i=1,7)
       WRITE(16,FMT='(/1x," First and last JSD values ")' ) 
       WRITE(16,FMT='(1x,9i8)')  1,(JSD(i, 1),i=1,8)
       WRITE(16,FMT='(1x,9i8)') NV,(JSD(i,NV),i=1,8)
       WRITE(16,FMT='(1x,8i8)') 

       CALL DATE_AND_TIME(CDate, CTime)
       WRITE(UNIT=16,FMT="(1x,' Loop start time',A10)") CTime 
       WRITE(UNIT= 6,FMT="(1x,' Loop start time',A10)") CTime 
!      GTM1 = REAL(CTime)

!!  End of rank 0 write out to message file.
       endif

!!  Work out number of procs on each rank for spectral loops
!!  Assume that nprocs <= NSpc, and npspc=MAX(1, NSpc/nprocs)
!     npstar =  myrank*npspc  
!     npsend = (myrank+1)*npspc - 1
!     IF( npsend .GT. NSpc - 1 ) npsend = NSpc - 1
!!  Replaced by sub Spctinit as in WW3 model.

!!  Work out number of procs on each rank for spatial loop
!!  Assume that nprocs <= NC-NPol, npsea=MAX(1, NSea/nprocs)
      npseatr =  myrank*npsea + 1  
      npseand = (myrank+1)*npsea 
      IF( npseand .GT. NSea ) npseand = NSea

!!  Allocate some temporary space for passing propagated spectral components.

      ALLOCATE( INTALLOC(3), INTALLO2(3, nprocs), STAT=malloc )
      ALLOCATE( REALandN(NCL, nprocs), REALLOC2(NSpc, nprocs), STAT=malloc )

!     Start of major time step loop
 TSLoop:  DO  NT=NS, NS+NTS

!!  Read current velocity and calculate gradients if available.
        IF ( FLCUR .AND. (MOD(NT, NHr) .EQ. 0) .AND. (myrank .EQ. 0) ) THEN
           CALL READCURNT(NT)
        ENDIF

!!  Broadcast current velocity to all ranks.
        IF ( FLCUR .AND. (MOD(NT, NHr) .EQ. 0) ) THEN
       CALL MPI_Bcast(CX(1), NSEA, MPI_REAL, 0, MPI_COM, ierr)
       CALL MPI_Bcast(CY(1), NSEA, MPI_REAL, 0, MPI_COM, ierr)
       CALL MPI_BARRIER(MPI_COM, ierr)

           CALL SMCDCXY
        ENDIF

!!    Space propagation for every spectral component.
!$ACC data copy(clats,WSpc,CGrp), &
!$ACC  copy(isd), &
!$ACC  copy(ice), &
!$ACC  copy(clatf), &
!$ACC  copy(jsd)

       IF ( FLCXY ) THEN

!!  Parallelised spectral loop for each rank
! SpcLop:  DO  NP=npstar, npsend
!!  Select assigned spectral components for own rank.

t1 = MPI_WTIME()

  SpcLop:  DO NP=1, NSpc
              IF( IAPPRO(NP) .EQ. myrank+1 ) THEN 
                  NF=(NP+NDIR-1)/NDir
                  ND=MOD(NP-1, NDir) + 1
!!    Propagation for the given spectral component
                  CALL W3PSMC( ND, NF, NT )           

              END IF
!!   End of IF( IAPPRO(NP) .EQ. myrank+1 ) block.

!!    End of spectral loops
           ENDDO  SpcLop


t2 = MPI_WTIME()
if (t2-t1>0) THEN
write (6,*) "Outer Time = ",(t2-t1)
end if

!!    Wait all ranks finish their spatial propagation for 
!!    all their assigned spectral components.
           CALL MPI_BARRIER(MPI_COM, ierr)

!!   Distribute propagation results to other ranks.
  Dislop:  DO NP=1, NSpc
              NN = IAPPRO(NP)
              NF=(NP+NDIR-1)/NDir
              ND=MOD(NP-1, NDir) + 1
!!    Save propagated spectral component in one array and broadcast to others
!$ACC data copy(REALandN(:,NN))
              write(*,*) "NN,ND, NF",nn,nd,nf
              IF( NN .EQ. myrank+1 ) THEN
!$ACC Kernels
                  REALandN(:,NN)=WSpc(ND, NF, :)
!$ACC End kernels
              ENDIF
!$ACC End data

!!    Wait the specific rank to finish its propagation results assignment.
              CALL MPI_BARRIER(MPI_COM, ierr)

!!   Broadcast to all other ranks from each rank.
              CALL MPI_Bcast(REALandN(1, NN), NCL , MPI_REAL, NN-1, MPI_COM, ierr)
!$ACC data copyin(REALandN(:,NN))
!$ACC Kernels
              IF( myrank + 1 .NE. NN ) THEN
                WSpc( ND, NF, : ) = REALandN(:, NN) 
              ENDIF
!$ACC End kernels
!$ACC End data

!!    Wait all ranks finish their propagation results storage.
              CALL MPI_BARRIER(MPI_COM, ierr)

!!    End of distribution loops
           ENDDO  DisLop
!!    End of spatial propagation FLCXY block.
       ENDIF 

!!    Wait all ranks finish their spatial propagation.
!      CALL MPI_BARRIER(MPI_COM, ierr)

!!    Refraction part for every cell except for polar cell(s).
       IF ( FLCTH .OR. FLCK ) THEN

! CelLop:  DO  NE=1, NC-NPol
!!  New distributed local loop for individual rank.  Note the subroutines 
!!  SMCkUNO2 and SMCGtCrfr called in W3KRTN are OpenMP parallelised.  JGLi10Jul2019
! !$ACC Kernels
! CelLop:  DO  NE=npseatr, npseand

!!    Great circle turning (GCT) for lower part cells at 4 times of substeps,
!!    Extended to include Arctic part if any
!      CALL GMCGtCrTn(1, NC, MRFct)
!      CALL W3KRTN( NE, NT )
       CALL W3KRTN( NT )

!!    End of refraction cell loop.
!          ENDDO  CelLop
! !$ACC End kernels
!!    Wait all ranks finish their spectral turning.
         CALL MPI_BARRIER(MPI_COM, ierr)

!!   Broadcast to all other ranks from each rank.
         DO MM = 0, nprocs-1
            IF( myrank .EQ. MM ) THEN
!!    Save propagated spectral component in one array and broadcast to others
                INTALLOC(1)=npseatr
                INTALLOC(2)=npseand
                INTALLOC(3)=(npseand - npseatr + 1)*NSpc
            ENDIF
            CALL MPI_Bcast(INTALLOC(1), 3, MPI_INTEGER, MM, MPI_COM, ierr)
!$ACC Update Host(WSpc)
            CALL MPI_Bcast(WSpc(1,1, INTALLOC(1)), INTALLOC(3),   &
       &                   MPI_REAL, MM, MPI_COM, ierr)
            CALL MPI_BARRIER(MPI_COM, ierr)
         ENDDO

!!    End of refraction FLCTH or FLCK block.
       ENDIF 
!$ACC End data

!!    Update boundary cells after proper rotation if Arctic part is
!!    included. 
       IF( Arctic ) THEN

!!    Arctic cells for global boundary cells
       DO i=1,NGLB
          ii=i+NGLA
          kk=MBGLo(i)
          DO j=1,NFrq

!!   Rotate the Arctic spectra by -AnglD before assigning to the lower part
!!   Note that it is equivalent to rotated the directional bins by AnglD.
          Spectr=WSpc(1:NDIR,j,kk)
          Alpha=  AngCD(kk)

          CALL Specturn( NDir, 1, Alpha, Spectr )

          WSpc(1:NDIR,j,ii)=Spectr

          ENDDO
       ENDDO

!!    Global cells for Arctic boundary cells
       DO i=1,NArB
          ii=i+NArA
          kk=MBArc(i)
          DO j=1,NFrq

!!   Rotate the lower part spectra by AnglD before assigning to the Arctic part
!!   Or turn the directional bins by -AnglD.   21 Jul 2009
          Spectr=WSpc(1:NDIR,j,kk)
!!   Note only the Arctic cells are assigned the AngCD value
          Alpha= -AngCD(ii)

          CALL Specturn( NDir, 1, Alpha, Spectr )

          WSpc(1:NDIR,j,ii)=Spectr

          ENDDO
       ENDDO

!!   End of updating boundary cells IF( Arctic ). 
       ENDIF


       if (myrank .eq. 0) then
!  Output tracer concentration if at selected writeup time steps
      IF( (NT+1 .LT. 10*NWP .AND. MOD(NT+1,NWP/2) .eq. 0) .OR.    &
     &    (MOD(NT+1,NWP) .eq. 0) )  THEN
        write(FL9NM(3:7), FMT='(i5)' )  10000+NT+1
        OPEN(UNIT=26, FILE=FL9NM, STATUS='NEW',IOSTAT=nn)
        IF(nn /= 0) PRINT*,' File FL9NM was not opened! '
        WRITE(UNIT=6,FMT='(2x,"NT= ",i6,3x,A9)') NT+1,FL9NM

           ii=0
        DO i=1, NC
           CTT=0.0
           DO j=1, NFrq
           DO k=1, NDIR
              CTT = CTT + WSpc(k,j,i)*DSIP(j)
           ENDDO
           ENDDO
           C(i)=CTT*DThta
           D(i) = SIGN( SQRT( ABS(C(i)) ), C(i) )

!!   Filter very small C(n) value so it is greater than E-36
           IF( Abs(D(i)) .LT. 1.0E-36 ) THEN
               D(i)=SIGN(1.0E-36, C(i))
           ENDIF

        ENDDO

!    All cells are saved 
        WRITE(UNIT=26, FMT='(2x,2i8)' )  NT+1, NC
        WRITE(UNIT=26, FMT=7113)  (D(n),  n=1, NC)

        CLOSE(26)
      ENDIF
 7113 FORMAT( 1x, 7ES11.3 )

!!  End of rank 0 write out at NT time step.
       endif

!!  End of time step loop
      ENDDO  TSLoop

       if (myrank .eq. 0) then
       CALL DATE_AND_TIME(CDate, CTime)
       WRITE(UNIT=16,FMT="(1x,' Loop ended time',A10)") CTime 
       WRITE(UNIT= 6,FMT="(1x,' Loop ended time',A10)") CTime 

!$     TM2= OMP_GET_WTIME()
!$     WRITE(UNIT= 6,FMT='(1x," Total loop time (s) = ",ES16.5)' ) TM2-TM1
!$     WRITE(UNIT=16,FMT='(1x," Total loop time (s) = ",ES16.5)' ) TM2-TM1

!!  End of rank 0 write out CTIME.
       endif

 9999  PRINT*, ' SMCPrDRG completed from rank ', myrank

       if (myrank .eq. 0) then
       CALL DATE_AND_TIME(CDate, CTime)
!      GTM2 = REAL(CTime)
       WRITE(UNIT= 6,FMT="(1x,' End time date ',A10,2x,A10)") CTime, CDate
       WRITE(UNIT=16,FMT="(1x,' End time date ',A10,2x,A10)") CTime, CDate
!      WRITE(UNIT= 6,FMT="(1x,' Net time date ',F10.1,A10)") GTM2-GTM1, CDate
!      WRITE(UNIT=16,FMT="(1x,' Net time date ',F10.1,A10)") GTM2-GTM1, CDate
       WRITE(UNIT=16,FMT="(1x)") 
!!  End of rank 0 write out CTIME.
       endif

!! Teminate MPI session
       CALL MPI_Finalize(ierr)

       END PROGRAM SMCPrDRG 
!  End of main program


! Subroutine that initialise frequecy and direction bins, depth related group speed, 
!    wave number, and refraction coefficient.       JGLi20Jun2017 
       SUBROUTINE SgmCgWnK
         USE Constants
         USE W3DISPMD
         IMPLICIT NONE
         REAL:: CNST, CNST0, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6
         LOGICAL:: WW3DIS = .true. 

!!   Setup frequency bands.
         CNST5=0.5*(XFR - 1.0/XFR)
          SIG(0)=Frqc0*2.0*Pie/XFR
         DSIP(0)=SIG(0)*CNST5
         DO n=1, NFrq+1
             SIG(n) = SIG(n-1)*XFR 
            DSIP(n) = SIG(n)*CNST5 
         END DO

!!   Assign water depth to HCel from ICE(5,:) integer values but DMIN.
!!   Set all depth half minimum depth for negative cells.
         HCel(-9:0) = 0.5*DMIN
         HCel(1:NC)=FLOAT( ICE(5,1:NC) )
         DO n=1, NC
            IF(HCel(n) .LT. DMIN) HCel(n)=DMIN 
         ENDDO

!!   Assign group speed to be zero for negative cells.  
         CGRP(0:NFrq+1,-9:0)=0.0

!!   Calculate group speed and refaction factor for all cells
         CNST0=1.0E-6

!!   Use WW3 interpolation scheme if selected
         IF( WW3DIS ) THEN
!!   Setup interpolation table for dispersion relationship. 
            CALL DISTAB
         ENDIF 

!!   Frequency loop
         DO k=0, NFrq+1

!      CNST=2.0*GeoPie*Frqcy
         CNST=SIG(k)
         CNST4=CNST*CNST/GRAV

!!   Cell loop
         DO n=1, NC
            CNST3=HCel(n)

!!  Use WW3 interpolation scheme 
          IF( WW3DIS ) THEN

!         Calculate wavenumbers and group velocities.
            CALL WAVNU1(CNST,CNST3, WNmk(k,n),CGrp(k,n))

!!  Othewise, use iteration scheme.
          ELSE

!!  Iteration to calculate kH value for this cell
            CNST1=CNST4/TANH(CNST4*CNST3)
            CNST2=CNST4/TANH(CNST1*CNST3)
            DO WHILE ( ABS(CNST2 - CNST1) .GT. CNST0 )
               CNST1=CNST2
               CNST2=CNST4/TANH(CNST1*CNST3)
            ENDDO
            
!!  Save wave number 
            Wnmk(k,n) = CNST2

            CNST1=CNST2*CNST3 
            CNST2=1.0/COSH(CNST1)
!!  Group speed
            CGrp(k,n)=(GRAV*0.5/CNST)*(TANH(CNST1)+CNST2*CNST2*CNST1)

!!  End of selected dispersion scheme
          ENDIF

!!  Refraction rate factor, without the dt/d(theta) factor
            REFR(k,n)=CNST/SINH(2.0*CNST3*Wnmk(k,n))

!!  End of both cell and frequency loops
         ENDDO
         ENDDO

! 999  PRINT*, ' Sub SgmCgWnK ended.'

       RETURN
       END SUBROUTINE SgmCgWnK


! Subroutine that assigns spectral components to MPI ranks for spatial propagation.
!    Adapted from WW3 model subroutine w3init in module w3initmd.ftn.
!    Wave speed is used to balance the computing load.    JGLi26Jun2019 
       SUBROUTINE Spctinit
         USE Constants
         IMPLICIT NONE
         REAL:: CNST, CNST0, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6
         INTEGER:: NTTOT, NTLOC, NTTARG, NTTMAX, NTTSpc(NSpc), ISTEP, ISP 
!
! 2.c.3 Calculated expected number of prop. calls per processor
!
         NTTOT = 0
         DO K=1, NFrq
            NTLOC = 1 + INT(DTG/(DTCFL*SIG(K)/SIG(1))-0.001)
            NTTOT = NTTOT + NTLOC*NDir
         END DO
         NTTARG = 1 + (NTTOT-1)/nprocs
         NTTARG = NTTARG + INT(DTG/(DTCFL*SIG(NFrq)/SIG(1))-0.001)
         NTTMAX = NTTARG + 5
!
! 2.c.4 Initialize IAPPRO and NTTSpc
!
         IAPPRO = 1
         NTTSpc = NTTOT
!
! 2.c.5 First sweep filling IAPPRO
!
        DO i = 1, nprocs
           ISTEP  = i
           ISP    = 0
           NTTSpc(i) = 0
           DO j=1, 1+NSpc/nprocs
             ISP    = ISP + ISTEP
             IF ( MOD(j,2) .EQ. 1 ) THEN
                ISTEP  = 2*(nprocs-i) + 1
             ELSE
                ISTEP  = 2*i - 1
             END IF
             IF ( ISP .LE. NSpc ) THEN
                IJ     = 1 + (ISP-1)/NDir
                NTLOC  = 1 + INT(DTG/(DTCFL*SIG(IJ)/SIG(1))-0.001)
                IF ( NTTSpc(i)+NTLOC .LE. NTTARG ) THEN
                    IAPPRO(ISP) = i
                    NTTSpc(i)   = NTTSpc(i) + NTLOC
                ELSE
                    IAPPRO(ISP) = -1
                END IF
             END IF
           END DO 
        END DO
!
! 2.c.6 Second sweep filling IAPPRO
!
        DO i=1, nprocs
           IF( NTTSpc(i) .LT. NTTARG ) THEN
              DO ISP=1, NSpc
                 IF( IAPPRO(ISP) .EQ. -1 ) THEN
                    IJ     = 1 + (ISP-1)/NDir
                    NTLOC  = 1 + INT(DTG/(DTCFL*SIG(IJ)/SIG(1))-0.001)
                    IF ( NTTSpc(i)+NTLOC .LE. NTTARG ) THEN
                        IAPPRO(ISP) = i
                        NTTSpc(i)   = NTTSpc(i) + NTLOC
                    END IF
                 END IF
              END DO
           END IF
        END DO
!

! 999  PRINT*, ' Sub Spctinit ended.'

       RETURN
       END SUBROUTINE Spctinit


! Subroutine that initialise C U V UC VC for rotating or deform 
!   flow field, using the same grid and time step.

      SUBROUTINE SPECUUVV
        USE Constants
        IMPLICIT NONE
        REAL:: CNST, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6, CNST7
        REAL:: Spec1D(NDIR), Spec2D(NDIR)
        REAL:: Spec1F(NFrq), Spec2F(NFrq)

!!  All velocities are in unit of basic cell length and time step or 
!!  divided by the grid velocity BX/DT or BY/DT

!! Initialise C U V for spherical rotation test in Arctic

!! Whole array assignment for WSpc(NDir,NFrq,NCL)
        WSpc=0.0

!! All spectra are identical to a cosine distribution with main direction
!! at 45 deg NE and total wave energy equal to 25.0 so SQRT(E)=5.0.
!! or at -45 deg SE
!       CNST=50.0/Pie
!! Wave spectral factor is doubled (for SWH0 ~ 7m) as frequency profile is 
!! normalised for full range so its limited range integration is < 1.
        CNST=100.0/Pie
        Spec1D = 0.0
        Spec2D = 0.0
        CNST6 = 0.0
        DO k=1, NDIR
           CNST1=COS(Theta(k) + Pie/4.0)
           CNST2=COS(Theta(k) - Pie/4.0)
           IF(CNST1 .GT. 0.0) THEN
              Spec1D(k)=CNST*CNST1*CNST1
           ENDIF
           IF(CNST2 .GT. 0.0) THEN
              Spec2D(k)=CNST*CNST2*CNST2
           ENDIF
           CNST6 = CNST6 + Spec1D(k)
        ENDDO
        SWH0=SQRT(CNST6*DThta)
        
!! Frequecy variable will be defined using the well-developed 
!! Pierson-Moskowitz (P-M) spectrum and normalised for full range
!! as (4*beta/omega)*(x^4)*EXP(-beta*x^4), where x=Ompeak/omega)
        CNST1=5.0
        CNST2=CNST1/4.0
        CNST3=SIG(NFrq/2)
        PkFrq=CNST3*0.5/Pie

!!   Frequency factors as P-M spectral.
        DO k=0, NFrq+1
           CNST=SIG(k)/CNST3
           CNST4=CNST*CNST*CNST*CNST 
           SpecPM(k)=(CNST1/CNST)*CNST4*EXP(-CNST2*CNST4)
        ENDDO
        Spec1F=SpecPM(1:NFrq)
        Spec2F = Spec1F

!! Initialise a strip below 45S-60S or j = -768 to -1024 in Southern Ocean 
!! and 45N to 60N or j = 768 to 1024 in upper part of the G6-25SMC domain.
!! Also put a non-zero spectral zone in the Arctic above 86.25N or 16 rows.
!! Note ICE(2,*) is on SW cell corner.
!! Equator is on cell face and ICE is equal to 0 at Equator.
!     ij=NINT(45.0/DLat)
!     jk=NINT(60.0/DLat)
!     mn=NINT(86.25/DLat)
!! Two initial wave spectral strips for Med36125 grid, bounded by a latitude 
!! upperlimit for the south and lowerlimit for the north strips. JGLi09Jul2019
      ij=NINT( (SpSouLat-ZrLAT)/DLat)
      kl=NINT( (SpNorLat-ZrLAT)/DLat)

!! Loop over all cells to select strip points. 
      DO i=1,NC
         ll=ICE(2,i)
         mm=ICE(2,i)+ICE(4,i)

!! Spectral frequency factor for all directions
         DO j=1,NFrq

!!  South and west boundaries use Spec2D
         IF( mm .LT. ij ) THEN
             WSpc(1:NDIR,j,i)=Spec2D*Spec2F(j)
             C(i)=SWH0
         ENDIF
!!  Northern boundary use Spec1D
         IF( kl .LT. ll ) THEN
             WSpc(1:NDIR,j,i)=Spec1D*Spec1F(j)
             C(i)=SWH0
         ENDIF

         ENDDO
      ENDDO

!     WRITE(6,*) '  Wind file conversion done!'

! 999  PRINT*, ' Sub SPECUUVV ended.'

      RETURN

      END SUBROUTINE SPECUUVV


!  This subroutine turn the wave spectrum by an fixed angle anti-clockwise
!  so that it may be used in the rotated or stanadard system.
!  First created:   26 Aug 2005   Jian-Guo Li
!  Last modified:   20 Jul 2009   Jian-Guo Li
!
! Subroutine Interface:

      Subroutine Specturn( NDirc, NFreq, Alphad, Spectr )
 
! Description:
!   Rotates wave spectrum anticlockwise by angle alphad
!
! Subroutine arguments
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: NFreq, NDirc         ! No. frequ and direc bins
       REAL,    INTENT(IN) :: Alphad               ! Turning angle in degree
       REAL, INTENT(INOUT) :: Spectr(NDirc,NFreq)  ! Wave spectrum in and out

! Local variables
       INTEGER :: ii, jj, kk, nsft
       REAL    :: Ddirc, frac, CNST
       REAL, Dimension(NFreq)      ::  Wrkfrq, Tmpfrq
       REAL, Dimension(NDirc,NFreq)::  Wrkspc

! Check input bin numbers
       IF( (NFreq .LT. 0) .OR. (NDirc .LT. 0) )  THEN
          PRINT*, " Invalid bin number NF or ND", NFreq, NDirc
          RETURN
       ELSE
          Ddirc=360.0/FLOAT(NDirc)
       ENDIF

! Work out shift bin number and fraction

      CNST=Alphad/Ddirc
      nsft=INT( CNST )
      frac= CNST - FLOAT( nsft )
!     PRINT*, ' nsft and frac =', nsft, frac

! Shift nsft bins if >=1
        IF( ABS(nsft) .GE. 1 )  THEN
      DO ii=1, NDirc

! Wave spectral direction bin number is assumed to increase clockwise from North
! So shift nsft bins anticlockwise results in local bin number increases by nsft
         jj=ii + nsft
 
! As nsft may be either positive or negative depends on alphad, wrapping may
! happen in either ends of the bin number train
         IF( jj > NDirc )  jj=jj - NDirc
         IF( jj < 1     )  jj=jj + NDirc

! Copy the selected bin to the loop bin number
         Wrkspc(ii,:)=Spectr(jj,:)
 
      Enddo

! If nsft=0, no need to shift, simply copy
        ELSE
        Wrkspc = Spectr
        ENDIF

! Pass fraction of wave energy in frac direction
! Positive or anticlock case, larger bin upstream
        IF( frac > 0.0 ) THEN
      Tmpfrq=Wrkspc(1,:)*frac
      DO kk=NDirc, 1, -1
         Wrkfrq=Wrkspc(kk,:)*frac 
         Spectr(kk,:)=Wrkspc(kk,:) - Wrkfrq + Tmpfrq 
         Tmpfrq=Wrkfrq
      ENDDO
        ELSE
! Negative or clockwise case, smaller bin upstream
      Tmpfrq=Wrkspc(NDirc,:)*frac
      DO kk=1, NDirc
         Wrkfrq=Wrkspc(kk,:)*frac
         Spectr(kk,:)=Wrkspc(kk,:) + Wrkfrq - Tmpfrq
         Tmpfrq=Wrkfrq
      ENDDO
        ENDIF

! Specturn completed

       Return 
       End Subroutine Specturn
!

! Subroutine that generates the Arctic reference direction angle
      SUBROUTINE ArctAngd
        USE Constants
        IMPLICIT NONE
        REAL:: CNST, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6
        REAL, ALLOCATABLE, DIMENSION(:)::  XLon, WLat, ELon, ELat, AnglD
        REAL ::  DfPolat, DfPolon, Omega, DfSpeed

!!    Note only the Arctic part needs the rotation angles.
!     Work out u-face central position XLon, WLat in standard grid
      CNST1=DLon*0.5
      CNST2=DLat*0.5
      DfPolat=Polat
      DfPolon=Polon

      ALLOCATE( XLon(NUAr), WLat(NUAr), ELon(NUAr), ELat(NUAr), AnglD(NUAr) )
      WRITE(6,*) " Calculating U component ..."

!!  Initialisation done before call this sub for Arctic part.
!        AngU=0.0
!        AngV=0.0
!        AngCD=0.0

      DO L=1, NUAr
         i=L+NUGL 
!!  U-face latitude with half dlat increase from SW corner
!!  Note j is from -220 to 239 and j=0 corresponds to v-face on the Equator
!!  Longitude is measured half-grid from ZrLon and first cell i=0 coicides with XLon=0.
         XLon(L)= Float( ISD(1,i) )*DLon 
         WLat(L)= Float( ISD(2,i) )*DLat + Float( ISD(3,i) )*CNST2 
      END DO

!!  Convert standard lat/lon into rotated lat/lon for deformation wind
      CALL LLTOEQANGLE( WLat, XLon, ELat, ELon,     &
     &                 AnglD, DfPolat, DfPolon, NUAr)

      DO L=1, NUAr
         i=L+NUGL 
!!  Convert the AnglD into rad and store in AngU(L). True U value for all 
!!  directional bins will be generated from the value later.
         AngU(i)=AnglD(L)*D2RAD 
      END DO

!!  Output AngU for checking
!     WRITE(6,        *)  "(AngU(L+NUGL), L=1, NUAr)"
!     WRITE(6,'(8ES12.3)') (AngU(L+NUGL), L=1, NUAr)

      DEALLOCATE( XLon, WLat, ELon, ELat, AnglD )

        ALLOCATE( XLon(NVAr), WLat(NVAr), ELon(NVAr), ELat(NVAr), AnglD(NVAr) )

!     Work out v-face central position XLon, WLat in standard grid
      DO L=1, NVAr
         j=L+NVGL
!!  V-face latitude with half_dlon*JSD(3) increase from SW corner
         XLon(L)= Float( JSD(1,j) )*DLon + CNST1*Float( JSD(3,j) ) 
         WLat(L)= Float( JSD(2,j) )*DLat
      END DO

!!  Convert standard lat/lon into rotated lat/lon for deformation wind
      CALL LLTOEQANGLE( WLat, XLon, ELat, ELon,     &
     &                 AnglD, DfPolat, DfPolon, NVAr)

      DO L=1, NVAr
         j=L+NVGL
!!  Convert the AnglD into rad and store in AngV(L). True V value for all 
!!  directional bins will be generated from the value later.
         AngV(j)= AnglD(L)*D2RAD
      END DO

!!  Output AngV for checking
!     WRITE(6,        *)  "(AngV(L+NVGL), L=1, NVAr)"
!     WRITE(6,'(8ES12.3)') (AngV(L+NVGL), L=1, NVAr)

      DEALLOCATE( XLon, WLat, ELon, ELat, AnglD )


!! Specific cell centre velocity compenents
      ALLOCATE( XLon(NArc), WLat(NArc), ELon(NArc), ELat(NArc), AnglD(NArc) )
      WRITE(6,*) " Calculating UC, VC component ..."

      CNST1=DLon*0.5
      CNST2=DLat*0.5
!! All cells include the polar cell
!! Note the wlat is not at 90N for the polar cell as direction will be undefined.
!! Here Wlat is half dlat from the polar cell edge and half dlat from the NP.
      DO L=1, NArc-1
         i=L+NGLo

!!  Cell centre latitude equal to west side centre latitude.
!!  Cell centre longitude with half cell width increase from West side centre
!!  Although the polar cell is of angular radius dlat (not dlat/2) the 
!!  transformation location is still used dlat/2 from its SW corner. The error
!!  will be negeligible as only the AnglD is used.
         XLon(L)= Float( ICE(1,i) )*DLon + CNST1*Float( ICE(3,i) ) 
         WLat(L)= Float( ICE(2,i) )*DLat + CNST2*Float( ICE(4,i) )

      END DO

!! North Polar cell centre coincide with NP
         XLon(NArc)=0.0
         WLat(NArc)=90.0
!! AnglD will be undefined with NP location as no local east at NP.

!!  Convert standard lat/lon into rotated lat/lon for deformation wind
      CALL LLTOEQANGLE( WLat, XLon, ELat, ELon,     &
     &                 AnglD, DfPolat, DfPolon, NArc )

      DO L=1, NArc
         i=L+NGLo
!!  Keep the AnglD in Deg and store in AngCD(L).  Spectral rotation for
!!  boundary cell update will use this angle later.
         AngCD(i)=  AnglD(L) 
!!  Save rotated latitude for refraction in Arctic part.  JGLi08Jul2015
         ELaCD(i)=  ELat(L) 
!!Li   Redefine GCT term factor for Arctic part or the netative of 
!!Li   tangient of rotated latitude divided by radius.  JGLi14Sep2015
         CTHG0S(i)= - TAN( ELat(L)*D2RAD ) / REARTH
      END DO

!!  Output AngCD for checking
!     WRITE(6,        *)  "(AngCD(L+NGLo), L=1, NArc)"
!     WRITE(6,'(8ES12.3)') (AngCD(L+NGLo), L=1, NArc)

 999  PRINT*, ' Sub ArctAngd ended.'

      RETURN

      END SUBROUTINE ArctAngd


!Li
!Li  Merged UM source code for rotated grid, consiting the following
!Li  original subroutines in UM 6.1
!Li    LLTOEQ1A  WCOEFF1A  and  LBCROTWINDS1
!Li  The last subroutine is modified to process only one level winds
!Li  cpp directives are removed and required header C_Pi.h inserted.
!Li	    Jian-Guo Li     26 May 2005
!Li
!Li  The WCOEFF1A subroutine is merged into LLTOEQ to reduce repetition
!Li  of the same calculations. Subroutine interface changed to 
!Li  LLTOEQANGLE
!Li	    Jian-GUo Li     23 Aug 2005
!Li
!LL  Subroutine LLTOEQANGLE--------------------------------------------    
!LL                                                                        
!LL  Purpose:  Calculates latitude and longitude on equatorial             
!LL            latitude-longitude (eq) grid used in regional               
!LL            models from input arrays of latitude and                    
!LL            longitude on standard grid. Both input and output           
!LL            latitudes and longitudes are in degrees.                    
!Li	       Also calculate rotation angle in degree to tranform
!Li            standard wind velocity into equatorial wind.
!Li	       Valid for 0<PHI_POLE<90 or new pole in N. hemisphere.
!LL                                                                        
!* Arguments:--------------------------------------------------------    
      SUBROUTINE LLTOEQANGLE( PHI, LAMBDA, PHI_EQ, LAMBDA_EQ,     &  
     &                 ANGLED, PHI_POLE, LAMBDA_POLE, POINTS )       

      IMPLICIT NONE 

      INTEGER:: POINTS    !IN  Number of points to be processed             

      REAL :: PHI_POLE,  & !IN  Latitude of equatorial lat-lon pole
     &        LAMBDA_POLE  !INOUT  Longitude of equatorial lat-lon pole

      REAL, DIMENSION(POINTS) ::         &
     &        PHI,       & !IN  Latitude
     &        LAMBDA,    & !IN  Longitude
     &        ANGLED,    & !OUT turning angle in deg for standard wind
     &        LAMBDA_EQ, & !OUT Longitude in equatorial lat-lon coords
     &        PHI_EQ       !OUT Latitude in equatorial lat-lon coords

! Define local varables:-----------------------------------------------
      REAL :: A_LAMBDA, A_PHI, E_LAMBDA, E_PHI, SIN_PHI_POLE, COS_PHI_POLE,  &
     &        TERM1, TERM2, ARG, LAMBDA_ZERO, LAMBDA_POLE_KEEP
      INTEGER   :: I   
      REAL, PARAMETER :: SMALL=1.0E-6

! Constants from comdecks:---------------------------------------------

      Real, Parameter :: Pi = 3.14159265358979323846  , &
     &                   Pi_Over_180 = Pi/180.0       , &
     &                   Recip_Pi_Over_180 = 180.0/Pi        

!*----------------------------------------------------------------------   

! 1. Initialise local constants
! Scale lambda pole to range -180 to 180 degs
      LAMBDA_POLE_KEEP=LAMBDA_POLE
      IF (LAMBDA_POLE.GT. 180.0) then
          LAMBDA_POLE=LAMBDA_POLE-360.0
      ENDIF

! Latitude of zeroth meridian
      LAMBDA_ZERO=LAMBDA_POLE+180.0
! Sine and cosine of latitude of eq pole
      IF (PHI_POLE >= 0.0) THEN
        SIN_PHI_POLE =  SIN(PI_OVER_180*PHI_POLE)
        COS_PHI_POLE =  COS(PI_OVER_180*PHI_POLE)
      ELSE
        SIN_PHI_POLE = -SIN(PI_OVER_180*PHI_POLE)
        COS_PHI_POLE = -COS(PI_OVER_180*PHI_POLE)
      ENDIF

! 2. Transform from standard to equatorial latitude-longitude

      DO 200 I= 1, POINTS

! Scale longitude to range -180 to +180 degs

      A_LAMBDA=LAMBDA(I)-LAMBDA_ZERO
      IF(A_LAMBDA.GT. 180.0) A_LAMBDA=A_LAMBDA-360.
      IF(A_LAMBDA.LE.-180.0) A_LAMBDA=A_LAMBDA+360.

! Convert latitude & longitude to radians

      A_LAMBDA=PI_OVER_180*A_LAMBDA
      A_PHI=PI_OVER_180*PHI(I)

! Compute eq latitude using equation (4.4)

      ARG=-COS_PHI_POLE*COS(A_PHI)*COS(A_LAMBDA)   &
     &    +SIN_PHI_POLE*SIN(A_PHI)
      ARG=MIN(ARG, 1.0)
      ARG=MAX(ARG,-1.0)
      E_PHI=ASIN(ARG)
      PHI_EQ(I)=RECIP_PI_OVER_180*E_PHI

! Compute eq longitude using equation (4.6)

      TERM1 = SIN_PHI_POLE*COS(A_PHI)*COS(A_LAMBDA)   &
     &       +COS_PHI_POLE*SIN(A_PHI)
      TERM2 = COS(E_PHI)
      IF(TERM2 .LT. SMALL) THEN
        E_LAMBDA=0.0
      ELSE
        ARG=TERM1/TERM2
        ARG=MIN(ARG, 1.0)
        ARG=MAX(ARG,-1.0)
        E_LAMBDA=RECIP_PI_OVER_180*ACOS(ARG)
        E_LAMBDA=SIGN(E_LAMBDA,A_LAMBDA)
      ENDIF

! Scale longitude to range 0 to 360 degs

      IF(E_LAMBDA.GE.360.0) E_LAMBDA=E_LAMBDA-360.0
      IF(E_LAMBDA.LT.  0.0) E_LAMBDA=E_LAMBDA+360.0
      LAMBDA_EQ(I)=E_LAMBDA

!Li  Calculate turning angle for standard wind velocity

      E_LAMBDA=PI_OVER_180*LAMBDA_EQ(I)

! Formulae used are from eqs (4.19) and (4.21)

      TERM2=SIN(E_LAMBDA)
      ARG= SIN(A_LAMBDA)*TERM2*SIN_PHI_POLE      &
     &    +COS(A_LAMBDA)*COS(E_LAMBDA)
      ARG=MIN(ARG, 1.0)
      ARG=MAX(ARG,-1.0)
      TERM1=RECIP_PI_OVER_180*ACOS(ARG)
      ANGLED(I)=SIGN(TERM1,TERM2)
!Li

 200  CONTINUE

! Reset Lambda pole to the setting on entry to subroutine
      LAMBDA_POLE=LAMBDA_POLE_KEEP

      RETURN
      END SUBROUTINE LLTOEQANGLE
!Li


! Subroutine to read cell and face arrays and to define grid variables.
!  First created:    1 Apr 2015   Jian-Guo Li
!  Last modified:    1 Apr 2015   Jian-Guo Li

      SUBROUTINE ReadCell
        USE Constants
        IMPLICIT NONE
        REAL:: CNST, CNST0, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6

!  Read Global and Arctic part Multiple-Cell info

       OPEN(UNIT=8, FILE=TRIM(CelPath)//TRIM(CelFile),  &
                    STATUS='OLD',IOSTAT=nn,ACTION='READ')
       IF(nn /= 0) PRINT*, TRIM(CelPath)//TRIM(CelFile)//' was not opened! '
          READ (8,*) NGLo, NRLCel(1:MRL) 
       DO J=1,NGLo
          READ (8,*) ICE(1,J), ICE(2,J), ICE(3,J), ICE(4,J), ICE(5,J)
       END DO
       CLOSE(8)
       PRINT*, CelFile//' read done ', NGLo, NRLCel(1:MRL) 

!!  Arctic part becomes optional.  JGLi12Dec2011
       IF( Arctic ) THEN

       OPEN(UNIT=8, FILE=TRIM(CelPath)//TRIM(ArcFile),  &
                    STATUS='OLD',IOSTAT=nn,ACTION='READ')
       IF(nn /= 0) PRINT*, ArcFile//' was not opened! '
          READ (8,*) NArc, NArB, NGLB
       DO J=NGLo+1, NGLo+NArc
          READ (8,*) ICE(1,J), ICE(2,J), ICE(3,J), ICE(4,J), ICE(5,J)
       END DO
       CLOSE(8)
       PRINT*, ArcFile//' read done  NArc=', NArc

!!  Total Arctic boundary cells
          NB=NArB+NGLB
          PRINT*, ' With Arctic part', NArc, NArB, NGLB

       ELSE
          NArc = 0
          NArB = 0
          NGLB = 0
          PRINT*, ' No Arctic part', NArc, NArB, NGLB
       ENDIF

!  Total cell number will be sum of two parts
       NC = NGLo + NArc

!!   Set boundary cell counts.  Boundary cells for the global part are at the end
!!   of SMC625Cels.dat and for the Arctic part at the start of SMC625Budy.dat.
!!   Boundary cell will then from NGLo-NGLB+1 to NGLo for lower part and NGLo+1 to NGLo+NArB
!!   NGLA and NArA are the extra numbers to be added for boundary loop 1, NGLB and 1, NArB
       NGLA=NGLo-NGLB
       NArA=NGLo

!!    Work out South and North Pole cell number if NPol = 2
      IF (NPol .EQ. 2) THEN
         IF( ICE(2,NC) .GT. ICE(2,NC-1) ) THEN
             NNor = NC
             NSou = NC - 1
         ELSE
             NSou = NC
             NNor = NC -1
         ENDIF
      ELSE
!!    Assume only North Pole in the Arctic is used (NPol = 1).
         NNor = NC
         NSou = 0
      ENDIF

!  Output a few to check input values
       DO J=1, NC, 10000
          WRITE(6,'(i8,2i6,2i4,i6)') J, ICE(1,J), ICE(2,J), ICE(3,J), ICE(4,J), ICE(5,J)
       END DO

!! Matching boundary cells for links between global and Arctic parts.
!! Assuming global boundary cells are at the end and Arctic boundary
!! cells are at the begginning of their cell list files.
       IF( Arctic ) THEN

!!   Match global boundary cells with Arctic inner cells
       DO i=1, NGLB
          ii=i+NGLA
!!   Search arctic part cells to match global part boundary cells
          mm=0
          DO k=NArA+NArB+1, NC-NPol
             IF(ICE(1,ii) .EQ. ICE(1,k) .AND. ICE(2,ii) .EQ. ICE(2,k)) THEN
                MBGLo(i)=k
                mm=1
             ENDIF
          ENDDO
          IF( mm .EQ. 0 ) PRINT*,' Miss global part boundary cell i=',i
       ENDDO

!!   Match Arctic boundary cells with global inner cells
       DO i=1, NArB
          ii=i+NArA
!!   Search global part to match arctic part boundary cells
          mm=0
          DO k=NGLA-2*NArB, NGLA
             IF(ICE(1,ii) .EQ. ICE(1,k) .AND. ICE(2,ii) .EQ. ICE(2,k)) THEN
                MBArc(i)=k
                mm=1
             ENDIF
          ENDDO
          IF( mm .EQ. 0 ) PRINT*,' Miss Arctic part boundary cell i=',i
       ENDDO
       PRINT*, ' Boundary cells matched for', NGLB, NArB

!!   End of boundary cell matching if (Arctic).
       ENDIF


!    Boundary -9 to 0 cells for cell size 2**n
!    Note the position indice for bounary cell are not used.
       ICE(1,-9:0)=0
       ICE(2,-9:0)=0
       ICE(3,   0)=1
       ICE(4,   0)=1
!!   Restrict boundary cell y-size no more than base cell size
!!          2**(MRL-1).
       mm = 2**(MRL-1)
       DO i=1,9
          ICE(3,-i)=ICE(3,-i+1)*2
          ICE(4,-i)=MIN(mm, ICE(3,-i))
       ENDDO

!!  Read sorted ISD JSD variables for global part.
      OPEN(UNIT=10, FILE=TRIM(CelPath)//TRIM(ISdFile),  &
                    STATUS='OLD',IOSTAT=nn,ACTION='READ')
      IF(nn /= 0) PRINT*, ISdFile//' was not opened! '
      READ(10,*) NUGL, NRLUFc(1:MRL)      
      WRITE(6,*) " Read u face numbers NUGL, NRLUFc(1:MRL)"     
      WRITE(6,*)                       NUGL, NRLUFc(1:MRL)      
      DO I=1,NUGL
         READ(10,*)  (ISD(N,I), N=1,7)
      END DO
      CLOSE(10)

      OPEN(UNIT=11, FILE=TRIM(CelPath)//TRIM(JSdFile),  &
                    STATUS='OLD',IOSTAT=nn,ACTION='READ')
      IF(nn /= 0) PRINT*, JSdFile//' was not opened! '
      READ(11,*) NVGL, NRLVFc(1:MRL)     
      WRITE(6,*) " Read v face numbers NVGL, NRLVFc(1:MRL) "
      WRITE(6,*)                       NVGL, NRLVFc(1:MRL)  
      DO J=1,NVGL
         READ(11,*)  (JSD(N,J), N=1,8)
      END DO
      CLOSE(11)

!!  Read sorted ISD JSD variables for Arctic part.
      IF( Arctic ) THEN

      OPEN(UNIT=10, FILE=TRIM(CelPath)//TRIM(AISFile),  &
                    STATUS='OLD',IOSTAT=nn,ACTION='READ')
      IF(nn /= 0) PRINT*, AISFile//' was not opened! '
      READ(10,*) NUAr
      WRITE(6,*) " Read u face numbers NUAr =", NUAr
      DO I=1,NUAr
         READ(10,*)  (ISD(N,I+NUGL), N=1,7)
      END DO
      CLOSE(10)

      OPEN(UNIT=11, FILE=TRIM(CelPath)//TRIM(AJSFile),  &
                    STATUS='OLD',IOSTAT=nn,ACTION='READ')
      IF(nn /= 0) PRINT*, AJSFile//' was not opened! '
      READ(11,*) NVAr
      WRITE(6,*) " Read v face numbers NVAr =", NVAr
      DO J=1,NVAr
         READ(11,*)  (JSD(N,J+NVGL), N=1,8)
      END DO
      CLOSE(11)

!!  Set total face nubmers
      NU=NUGL+NUAr
      NV=NVGL+NVAr

!!  Reset arctic part cell numbers in I/JSD by adding NGLo for positive cells only.
!!  The 0 and negative cells for boundary useage will be shared by the two parts.
      DO I=NUGL+1, NU
         DO M=4,7
            IF(ISD(M,I) > 0) ISD(M,I)=ISD(M,I)+NGLo
         END DO
      END DO

      DO J=NVGL+1, NV
         DO M=4,7
            IF(JSD(M,J) > 0) JSD(M,J)=JSD(M,J)+NGLo
         END DO
      END DO

      WRITE(6,*) " Arctic u v face cell values have been adjusted."

!!  Without Arctic part, set total NU NV equal to global part.
      ELSE

          NUAr = 0
          NVAr = 0
          NU=NUGL
          NV=NVGL

      ENDIF

 999  PRINT*, ' Sub READCELL ended.'

      RETURN

      END SUBROUTINE READCELL


! Subroutine to read current velocity and calculate their gradients.
!  First created:    20 Apr 2017   Jian-Guo Li
!  Last modified:    23 Nov 2017   Jian-Guo Li

      SUBROUTINE ReadCurnt(NTM)
        USE Constants
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NTM         ! No. of time steps 

        REAL:: CNST, CNST0, CNST1, CNST2, CNST3, CNST4, CNST5, CNST6
        INTEGER :: NwDay, NHour
!       CHARACTER(Len=20) :: CurnFile='Curnt100000.dat'
        CHARACTER(Len=20) :: CurnFile='ww3.14071700.cur'
        CHARACTER(Len=30) :: CurnPath='/data/d02/frjl/SMCUK12cur/'

        NHour=NTM/NHR
        NwDay=NDay+(NHour/24)*100
!  Read current velocity from ww3 output text file.
!       WRITE(CurnFile(6:11), FMT='(i6)' )  100000+NTM
        WRITE(CurnFile(9:12), FMT='(i4.4)' )  NwDay+MOD(NHour, 24) 
!       OPEN(UNIT=8, FILE=TRIM(CurnFile), STATUS='OLD',IOSTAT=nn,ACTION='READ')
        OPEN(UNIT=8, FILE=TRIM(CurnPath)//TRIM(CurnFile),  &
             STATUS='OLD',IOSTAT=nn,ACTION='READ')
        IF(nn /= 0) PRINT*, CurnFile//' was not opened! '
          READ (8,FMT='(30X,I9)') NSEA 
        IF(NC .NE. NSEA) THEN 
          PRINT*, " *** Inconsistent cell number NC, NSEA =", NC, NSEA
          Return
        ENDIF
          READ (8,*) (CX(I), I=1,NSEA), (CY(J), J=1,NSEA)
        CLOSE(8)

 999  PRINT*, ' Sub READCURNT ended for NT =', NTM

      RETURN

      END SUBROUTINE READCURNT

