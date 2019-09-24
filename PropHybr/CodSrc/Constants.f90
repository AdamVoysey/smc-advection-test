!!
!  This module is for SMC Med36125 grid setup and model configuration. 
!  It is for a propagation test and is written in FORTRAN 90 format. 
!                     J G Li   16 Jul 2019
!!
                      
      MODULE Constants
          USE mpi
          USE omp_lib
          IMPLICIT NONE

! Parameters fixed in the program
       INTEGER,PARAMETER::NCL=30000, NFC=36000, NBDY=100, NPol=0,   &
          &               NLat=6144, NLon=8192, MNDPTH=5, MRL=4,    &
          &               NDir=36,   NFrq=30,   NSpc=NDir*NFrq  

       REAL,PARAMETER:: XFR=1.1, Frqc0=0.04118, CTMAX=0.7, DTME=3600.0, &
      &                 CLARMN=85.0, DMIN=10.0, Refran=36.0    

       INTEGER :: NX, NY, NSEA 
       REAL :: SX, SY 

       REAL,PARAMETER:: ZrLAT=0.0, ZrLON=0.0, DLON=0.0439453125, DLAT=0.029296875, &
         &              PoLAT=0.0, PoLON=-180.0, SpSouLat=35.0, SpNorLat=40.0,   &
         &              Pie=3.141592654, RAD2D=180.0/Pie, D2RAD=Pie/180.0, PCRDY=1.0 

!      REAL,PARAMETER:: DTG=300.0, DTR=1.0/DTG, DTCFL=60.0, AKH=1000.0
       REAL,PARAMETER:: DTG=900.0, DTR=1.0/DTG, DTCFL=300.0, AKH=6000.0

!  Writeup interval and model run time in hours.
       INTEGER,PARAMETER:: NHr=INT(3600.0/DTG), NWP=2*NHr, NTS=1*NHr, NDay=1700 

!  Some physical and atmospheric constants
       REAL,PARAMETER:: GRAV=9.806,CPVAP=1004.5,RDRY=287.05, &
      &                 CT0=273.16,CALJO=4.1868,PATM=101325.0,ANGUL=7.2921E-5,  &
      &                 EPSLN=0.6220,CLIGHT=2.99792458E8, GeoPie=3.141592654,   &
      &                 REARTH=6.371E6, Agu36=4.8481E-5 

! Array variables to be used for data storage
       REAL::  AMG, CMX, CTT, UMX, DY, DYR, DX0, DTH, DThta, SWH0, Alpha
       REAL::  AKHDT2, CGCMX, CRFMX, PkFrq
       REAL, DIMENSION(-9:NCL):: A, C, D, F, AU, AV, DX, DXR, RCELA
       REAL, DIMENSION(-9:NCL):: HCel, DHDX, DHDY, AngCD, ELaCD, CLats, CTHG0S 
       REAL, DIMENSION(-9:NCL):: DW, CX, CY, DCXDX, DCXDY, DCYDX, DCYDY
       REAL, DIMENSION( NBDY ):: MBGlo, MBArc
       REAL, DIMENSION( NDir ):: Theta, ESIN, ECOS, EC2, ESC, ES2, Spectr, SpeGCT
       REAL, DIMENSION( 0:NFrq+1 )::  SIG, DSIP, SpecPM 
       REAL, DIMENSION( 0:NFrq+1, -9:NCL)::  REFR, CGrp, Wnmk
       REAL, DIMENSION(NFC)::   AngU, AngV, CLatF
       REAL, DIMENSION(NDir,NCL)::   DHLMT
       REAL, DIMENSION(NDir,NFrq,NCL):: WSpc

       REAL,    Allocatable, Dimension(:)  ::  REALLOC1
       REAL,    Allocatable, Dimension(:,:)::  REALLOC2, REALandN
       INTEGER, Allocatable, Dimension(:)  ::  INTALLOC
       INTEGER, Allocatable, Dimension(:,:)::  INTALLO2

       INTEGER:: icl, jcl, iuf, juf, ivf, jvf, LvR, Lvm 
       INTEGER:: NS, NT, ND, NE, NF, NA, NB, NP, NR, MRFct 
       INTEGER:: NC, NU, NV, NGLo, NGLA, NGLB, NArc, NArA, NArB,    &
      &          NUGL, NUAr, NVGL, NVAr, NNor, NSou
       INTEGER, DIMENSION(5,-9:NCL)::  ICE
       INTEGER, DIMENSION(7,NFC)::  ISD
       INTEGER, DIMENSION(8,NFC)::  JSD
       INTEGER, DIMENSION(0:MRL)::  NRLCel, NRLUFc, NRLVFc
       INTEGER, DIMENSION( NSpc)::  IAPPRO 
       INTEGER:: I,II,IJ,IJK,J,JJ,JK,K,KK,KL,L,LL,LM,LMN,M,MM,MN,N,NN

!  MPI related varialbes
       INTEGER:: MPI_COM, ierr, malloc, myrank, nprocs, nthreads,   &
      &       npspc, npstar, npsend, npsea, npseatr, npseand, provided
       INTEGER:: mpistat(MPI_STATUS_SIZE)         ! MPI status for Recv
       INTEGER,PARAMETER:: required=MPI_THREAD_FUNNELED 
       CHARACTER(len=MPI_MAX_PROCESSOR_NAME):: pname   ! Processor name

       LOGICAL:: Arctic = .false., FVERG = .false., FLCUR = .false., &
                  FLCXY = .true.,  FLCTH = .true.,  FLCK = .true.

       CHARACTER(LEN=9):: FL9NM='Cn10000.d'

!      Date and time for timing of program by calling Date_And_Time
       CHARACTER(LEN=10):: CDate, CTime

!  Cell and face array files.
       CHARACTER(Len=256) :: CelPath='./'
       CHARACTER(LEN=26)::  CelFile='Med36125Cel0.dat', &
        &                   ISdFile='Med325GISide.dat', &
        &                   JSdFile='Med325GJSide.dat', &  
        &                   ArcFile, AISFile, AJSFile       

!      CHARACTER(LEN=20):: RUNDATE=' SMC24816  2Sep2014 '
       CHARACTER(LEN=20):: RUNDATE=' Med36125  9Jul2019 '

      END MODULE Constants

