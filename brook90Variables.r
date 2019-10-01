NPINT <- 1         
SUBDAYDATA <- FALSE
RRD <- 0.55
UWD <- 3.0

# rmarkdown::render(file.path(projectpath,"Rmd_files",'GLOBDECL.Rmd'),output_dir=output_html)
## ------------------------------------------------------------------------
ML <- 25    
gctMaxoutvars <- 60 
gctGenoutv <- 19 
gctIntoutv <- 22 
gctLayoutv <- 12 
gvals <- 8 
maxgraphs <- 100
#
#A <- numeric(ML) # parsed values from INSTRNG
#aparsed  # parsed values from gsbparsestring, first subscript is 0
#bad      # error indicator passed back to calling routine
#B90path$  # path for running program, used for .INI and .TMP
#cancelsave As Boolean #true if save or print should be cancelled
#chkout(1 To gctMaxoutvars, 1 To 5) As Integer #1 for output, 5 intervals
#        indicates which variables/interval selected
#        stored in .INI
#commadec As Boolean
#dfiletype    #data file type, 1-csv,2-decimal period, 3-decimal comma
#
dum   <-0  
errnum <-0   
#
#EV(1 To 25) As New frmevalout # EVAL windows
#
evalnumber <-0 
#
#FG(1 To maxgraphs) As New frmgraph  # multiple instances of graph
#gnumber   # graph number
#graphon    # true if graph has been initialized
#
INIDAYS<-0  
#
#inoutdir$  # latest input-output directory
#intrvl     #output interval index, 1 ANN, 2 MON, 3 DAY, 4 PRE, 5 ITR
#ivar(1 To gctMaxoutvars), pvar(1 To gctMaxoutvars), dvar(1 To gctMaxoutvars), mvar(1 To gctMaxoutvars), yvar(1 To gctMaxoutvars) # list of var numbers for output, does not count mulltiple soil layers
#
lstart <-0
lend <-0 
lstep <-0 
#
#msg As String#
#noenderror As Integer  # 0 if no END in parameter file, -1 if END is found
#noruncontinue #1 to prevent Run-Continue after soil parameter changes
#NLN As String  # new line
#outfilestatus(1 To 5) As Integer  #0 if no output, 1 if output, subscript 1 ann, 2 mon, 3 day, 4 pre, 5 itr
##        indicates presence of .TMP output files, controls mnuview???.Enabled
#outselect(1 To 5)  As Integer # True(-1) if output selected for that time interval, False(0) if not
#        controls checking in Select Output menu
#
parserr <-0  
prfiletype  <-0  
rerunnumber<-0 
rstop <-0 
RUNDAYS<-0  
runflag <-0
#
#strng$
#Title$(1 To gctMaxoutvars)   # column titles for gctMaxoutvars variables
#txt$(1 To 30) # text strings for help
#userfont As String
#varno # output variable number, from 1 to gctMaxoutvars, always locally used
#yvars, mvars, dvars, pvars, ivars  # number of vars/cols in output


## ------------------------------------------------------------------------
DT <- 1
WTOMJ <- 0.0864
ETOM <- 0.4085
CPRHO <- 1240
GAMMA <- 0.067
CVLQ <- 0.00418
CVICE <- 0.00192
LF <- 0.335
LS <- 2.824
RHOWG <- 0.00981
SIGMA <- 0.0000000567
SC <- 1367
K <- 0.4
PI <- 3.1416


## ------------------------------------------------------------------------
AA  <-0       
ADEF<-0   
ALB   <-0.07       
ALBEDO<-0    
ALBSN <-0.3           
ALPHA<-numeric(ML) 
ASPECT  <-0         
ASPECTD <-0       
ASUBS   <-0          
ATR<-numeric(2)          
ATRANI<-numeric(ML)  
ATRI<-matrix(0,2,ML) 
AWAT    <-0          
BALERD<-0 
BALERM<-0  
BALERY<-0 
BEXP<-rep(0,ML)  
BEXP[1]<-5.37433
BEXP[2]<-4.03320
BEXP[3]<-5.64096
BEXP[4]<-6
BEXP[5]<-5
BYFL   <-0       
BYFLI<-numeric(ML)  
BYFLPI<-numeric(ML)
BYFLDI<-numeric(ML)
BYFLMI<-numeric(ML)
BYFLYI<-numeric(ML)
BYFLP<-rep(0,ML)
BYFLD<-rep(0,ML) 
BYFLM<-rep(0,ML)  
BYFLY<-rep(0,ML)
BYFRAC<-numeric(ML)  
BYPAR<-1            
#C    <-0             
C1    <-0.25
C2    <-0.5     
C3    <-0.2            
CC      <-0   
CCFAC  <-0.3          
CHM<-numeric(ML)     
CHN<-numeric(ML)     
CINTRL<-0.15            
CINTRS<-0.15            
CINTSL<-0.6         
CINTSS<-0.6             
CR   <-0.5              
CS  <-0.035           
CVPD <-2           
CZR  <-0.05              
CZS <-0.13          
DAYLEN <-0           
DAYMO<-numeric(12)     
DD  <-0            
DELTA <-0           
DENSEF <-0.8580571           
DISP  <-0        
DISPC  <-0           
DOM  <-0        
DOY  <-0            
DPSIDW<-numeric(ML) 
DPSIMX  <-0.01          
DRAIN  <-0.626259     
DSFL  <-0         
DSFLI<-numeric(ML)   
DSFLP<-0 
DSFLD<-0 
DSFLM<-0 
DSFLY<-0   
DSFLPI<-numeric(ML)
DSFLDI<-numeric(ML)
DSFLMI<-numeric(ML)
DSFLYI<-numeric(ML)
DSLOPE  <-0         
DSLOPED  <-0      
DSWMAX   <-2        
DTI   <-0          
DTIMAX  <-0.5          
DTINEW  <-0.0          
DTP    <-1           
DTRI  <-0             
DUMM<-rep(0,gctMaxoutvars)  
dummy   <-0          
DURATN<-c(4,4,4,4,4,4,4,4,4,4,4,4)   
EA    <-0             
ES   <-0              
ESLOPE  <-0           
ESLOPED  <-2         
EVAPP<-0 
EVAPD<-0  
EVAPM<-0  
EVAPY <-0
FARR<-rep(0,366)    
FETCH   <-5000           
FLOWP<-0 
FLOWD<-0 
FLOWM<-0 
FLOWY<-0  
FRINTL  <-0.06           
FRINTS  <-0.06        
FSINTL  <-0.04           
FSINTS  <-0.04          
FXYLEM  <-0.5           
GER<-numeric(2)            
GEVP  <-0       
GIR<-numeric(2)          
GIVP <-0       
GLMAX  <-0           
GLMAXC  <-0.494377          
GLMIN <-0           
GLMINC  <-0.03        
GRAPH <-0            
GRDMLT <-0.35           
GSC   <-0            
GSP  <-0.085              
GWAT   <-0            
GWATIN <-20            
GWFL  <-0             
GWFLP<-0 
GWFLD<-0 
GWFLM<-0 
GWFLY<-0   
mnuhalfiter<-FALSE
HEIGHT  <-0          
HR   <-10              
HS    <-1            
#I    <-0            
I0HDAY <-0            
IDAY<-0             
II <-0               
ILAYER <-0          
IDEPTH <-1000            
IMPERV  <-0.025 
INFEXP  <-0.8797636          
INFRAC<-rep(0,ML)  
INFLI<-rep(0,ML) 
INFLP<-rep(0,ML)
INFLD<-rep(0,ML)
INFLM<-rep(0,ML)
INFLY<-rep(0,ML) 
INFLPI<-rep(0,ML)
INFLDI<-rep(0,ML)
INFLMI<-rep(0,ML)
INFLYI<-rep(0,ML)
INTR   <-0           
INTRIN  <-0           
INTS   <-0            
INTSIN <-0            
IRVP <-0              
IRVPD<-0 
IRVPM<-0 
IRVPY <-0             
ISVP  <-0           
ISVPD<-0 
ISVPM<-0 
ISVPY<-0  
J   <-0              
KF<-rep(0,ML)      
KF[1]<-6.9
KF[2]<-2.7
KF[3]<-2.9
KF[4]<-1
KF[5]<-3.5
KK<-rep(0,ML)      
KSAT<-rep(0,ML)    
KSNVP <-0.009734            
L1  <-0               
L2  <-0               
LAI  <-0              
LAIMLT <-0.2           
LAT   <-0             
LATD <-50.5              
LENGTH  <-0           
LPC  <-4              
LWIDTH    <-0.004      
MARR<-c(seq(1,366,1))    
MAXHT<-25             
MAXLAI  <-7.693270           
MAXLQF <-0.05            
MELFAC <-1.728930            
MESFL <-0            
MESFLD<-0
MESFLM<-0
MESFLY<-0 
MESFLP <-0           
MM    <-0     
MONTHN <-0        
MXKPL  <-7.03463      
MXRTLN <-3000.001    
N   <-0            
NN  <-2.5             
NDAYS<-0            
NITS<-0             
NITSD<-0            
NITSM<-0           
NITSY<-0            
NITSR<-0            
NLAYER  <-5         
NOOUTF <-1          
NTFLI<-rep(0,ML)  
NTFLPI<-rep(0,ML)
NTFLDI<-rep(0,ML)
NTFLMI<-rep(0,ML)
NTFLYI<-rep(0,ML)
PINT <-0        
PINTD<-0 
PINTM<-0 
PINTY <-0 
PIR<-c(0,0)          
PREC  <-0      
PRECD<-0 
PRECM<-0 
PRECY<-0             
PREINT  <-0          
PRECIN  <-0          
PSICR  <--2           
PSIF<-rep(0,ML)   
PSIF[1]<--11.818
PSIF[2]<--11.516 
PSIF[3]<--10.22
PSIF[4]<--10
PSIF[5]<--10
PSIG<-rep(0,ML)     
PSIM<-rep(0,ML)      
PSIMIN<-rep(-10,ML)  
PSITI<-rep(0,ML)    
PSNVP <-0          
PTR<-c(0,0)          
PTRAN  <-0          
PTRAND<-0 
PTRANM<-0 
PTRANY<-0          
QFFC<-0.00104        
QFPAR <-0.834524        
QLAYER <-1        
QDEPTH <-0.1          
R5    <-100           
RAA    <-0           
RAC   <-0            
RAS  <-0           
RELHT<-c(1,1,366,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)    
RELLAI<-c(1,1,366,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)   
ROOTDEN<-c(50,1,50,1,50,1,50,1,50,1,50,.3,50,.2,50,.1,50,.1,50,.1,50,0.0,50,0.0,50,0,50,0,50,0,50,0,50,0,50,0,50,0,50,0,50,0,50,0,50,0,50,0,50,0)  
RELDEN<-rep(0,ML)  
RFAL   <-0          
RFALD<-0
RFALM<-0
RFALY <-0            
RHOTP <-2            
RINT  <-0           
RINTD<-0
RINTM<-0
RINTY<-0           
RM <-1000               
RNET   <-0          
RNETD<-0
RNETM<-0 
RNETY <-0            
RPLANT<-0            
RROOTI<-rep(0,ML)    
RTHR   <-0           
RTHRD<-0
RTHRM<-0
RTHRY<-0            
RTLEN <-0           
RTRAD <-0.35           
RRD<-0.55           
RSC <-0              
RSNO <-0           
RSNOD<-0
RSNOM<-0
RSNOY<-0         
RSS<-0              
RSSA <-500             
RSSB <-1             
RSTEMP<- -1.29978         
RXYLEM <-0           
SAFRAC <-0           
SAI <-0              
SAIMLT <-0.5           
SEEP <-0             
SEEPP<-0
SEEPD<-0
SEEPM<-0
SEEPY<-0             
SFAL <-0     
SFALD<-0
SFALM<-0
SFALY<-0             
SHEAT<-0        
SINT <-0           
SINTD<-0 
SINTM<-0
SINTY <-0
SLFDAY  <-0          
SLFL  <-0            
SLFLP<-0
SLFLD<-0
SLFLM<-0
SLFLY<-0
SLFLI<-rep(0,ML)  
SLFLPI<-rep(0,ML)
SLFLDI<-rep(0,ML)
SLFLMI<-rep(0,ML)
SLFLYI<-rep(0,ML)
SLRAD <-0            
SLRADd<-0
SLVP <-0            
SLVPD<-0
SLVPM<-0
SLVPY<-0          
SMLT  <-0            
SMLTD<-0
SMLTM<-0 
SMLTY <-0           
SNODEN <-0.3         
SNOEN <-0           
SNOFRC <-0        
SNOW <-0             
SNOWIN  <-20          
SNOWLQ<-0            
SNVP  <-0            
SNVPD<-0
SNVPM<-0
SNVPY <-0            
SOLRAD <-0           
SOLRADC<-0           
SRFL <-0             
SRFLP<-0
SRFLD<-0
SRFLM<-0
SRFLY <-0
STHR <-0             
STHRD<-0
STHRM<-0
STHRY<-0 
STONEF<-rep(0.00,ML)
STONEF[1]<-0.02
STONEF[2]<-0.2
STONEF[3]<-0.25
STONEF[4]<-0.25
STONEF[5]<-0.7
STORD<-0
STORM<-0
STORY<-0 
STRES<-0            
STRX<-0             
SWAT<-0    
SWATI<-rep(0,ML)  
SWATMX<-rep(0,ML)  
SWATQF   <-0         
SWATQX <-0           
T1    <-10           
T2    <-30           
TA   <-40            
TADTM <-0           
TAJ  <-0             
TANTM <-0            
TEMP  <-0           
TH   <-40            
THETA<-rep(0,ML)   
THETAF<-rep(0,ML) 
THETAF[1]<-0.34062341
THETAF[2]<-0.39705807  
THETAF[3]<-0.24359704
THETAF[4]<-0.35
THETAF[5]<-0.23
THICK<-rep(0,ML)   
THICK[1]<-50
THICK[2]<-250
THICK[3]<-300
THICK[4]<-300
THICK[5]<-100
THSAT<-rep(0,ML)  
THSAT[1]<-0.6313342
THSAT[2]<-0.6189386
THSAT[3]<-0.4930716
THSAT[4]<-0.680
THSAT[5]<-0.600
TL    <-0           
TMAX  <-0            
TMIN <-0             
TRANI<-rep(0,ML)   
TRANP<-0
TRAND<-0
TRANM<-0
TRANY<-0             
TRANPI<-rep(0,ML)
TRANDI<-rep(0,ML)
TRANMI<-rep(0,ML) 
TRANYI<-rep(0,ML)
TSNOW  <-0          
UA     <-0           
UADTM  <-0           
UAJ   <-0           
UANTM <-0          
UW   <-0             
VPD  <-0             
VRFLI<-rep(0,ML)  
VRFLPI<-rep(0,ML)
VRFLDI<-rep(0,ML)
VRFLMI<-rep(0,ML)
VRFLYI<-rep(0,ML)
VV<-rep(0,ML)      
WETC<-rep(0,ML)    
WETF<-rep(0,ML)   
WETFR  <-0         
WETINF<-rep(0,ML)  
WETINF[1]<-0.92
WETINF[2]<-0.92
WETINF[3]<-0.92
WETINF[4]<-0.92
WETINF[5]<-0.92
WETNES<-rep(0,ML)  
WNDRAT <-.3           
XMAX  <-0            
XMIN  <-0           
YEARN <-0           
YMAX<-numeric(gvals)      
YMIN<-numeric(gvals)      
YNAME<-numeric(gvals)    
YVAL<-numeric(gvals)      
YY   <-0            
Z0   <-0             
Z0C  <-0             
Z0G  <-0.02             
Z0GS  <-0            
Z0S  <-0.001             
Z0W  <-2.8             
ZA   <-0            
ZMINH <-2           
ZW  <-14              
