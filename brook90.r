## ----chunkpath-----------------------------------------------------------
#projectpath<-"."
#output_html<-file.path(projectpath,"Documentation\\HTML_Files")


## ----chunkmain, message=FALSE, warning=FALSE, results='hide'-------------
# rmarkdown::render('MainProg.Rmd',output_dir = output_html)
## ----chunk1--------------------------------------------------------------
#rm(list = setdiff(ls(), lsf.str()))


## ----chunk2--------------------------------------------------------------
NPINT <- 1         
SUBDAYDATA <- FALSE
RRD <- 0.55
UWD <- 3.0


## ----chunk3, message=FALSE, results='hide'-------------------------------
#projectpath<-"."
#output_html<-file.path(projectpath,"Documentation/HTML_Files")

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


# rmarkdown::render(file.path(projectpath,"Rmd_files",'AXX.Rmd'),output_dir = output_html)
## ------------------------------------------------------------------------
ACCUM <-function(A1, A2, A3, A4, A5, B1, B2, B3, B4, B5){
B1 <- B1 + A1
B2 <- B2 + A2
B3 <- B3 + A3
B4 <- B4 + A4
B5 <- B5 + A5
return(list(B1, B2, B3, B4, B5))
}


## ------------------------------------------------------------------------
ACCUMI<-function (N, A1, A2, A3, A4, A5, B1, B2, B3, B4, B5){
  B1 <- B1 + A1
  B2 <- B2 + A2
  B3 <- B3 + A3
  B4 <- B4 + A4
  B5 <- B5 + A5
#for( i in 1:N){
#  B1[i] <- B1[i] + A1[i]
#  B2[i] <- B2[i] + A2[i]
#  B3[i] <- B3[i] + A3[i]
#  B4[i] <- B4[i] + A4[i]
#  B5[i] <- B5[i] + A5[i]
#}
return(list(B1, B2, B3, B4, B5))
}


## ------------------------------------------------------------------------
ACOSF<-function (T){
TA<-0
AC<-0
TA <- abs(T)
if (TA > 1) {

}
if (TA < 0.7) {
  AC <- 1.570796 - atan(TA / (1 - TA * TA)^(1/2))
}else{
  AC <- atan((1 - TA * TA)^(1/2) / TA)
}
if (T < 0) {
  ACOS <- 3.141593 - AC
}else{
  ACOS <- AC
}
return(ACOS) #acos(T))
}


## ------------------------------------------------------------------------
ASINF<-function (temp){
TA<-0
TA <- abs(temp)
if(TA > 1){
#  
}
if (TA < 0.7) {
  ASIN <- sign(temp) * (atan(TA / (1 - TA * TA)^(1/2)))
}else{
  ASIN <- sign(temp) * (1.570796 - atan((1 - TA * TA)^(1/2) / TA))
}
return(ASIN)#asin(temp))
}


## ------------------------------------------------------------------------
RMAXF <-function(T1, T2){
if (T1 < T2) {
  RMAX <- T2
}else{
  RMAX <- T1
}
return(RMAX)
}


## ------------------------------------------------------------------------
RMINF<-function (T1, T2){
if (T1 > T2) {
RMIN <- T2
}else{
RMIN <- T1
}
return(RMIN)
}


## ------------------------------------------------------------------------
SUMI<-function (N, A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6){
B<-rep(0,N)
for( i in 1:N){
B[1] <- B[1] + A1[i]
B[2] <- B[2] + A2[i]
B[3] <- B[3] + A3[i]
B[4] <- B[4] + A4[i]
B[5] <- B[5] + A5[i]
B[6] <- B[6] + A6[i]
}
return(list(B[1], B[2], B[3], B[4], B[5], B[6]))
}


## ------------------------------------------------------------------------
ZERO<-function (V1, V2, V3, V4, V5, V6){
  V1 <- 0
  V2 <- 0
  V3 <- 0
  V4 <- 0
  V5 <- 0
  V6 <- 0
    return(list(V1, V2, V3, V4, V5, V6))
}


## ------------------------------------------------------------------------
ZEROA <-function(N,A1,A2,A3,A4){
#zeroes arrays
#for( i in 1:N){
  A1<-rep(0,ML)
  A2<-rep(0,ML)
  A3<-rep(0,ML)
  A4<-rep(0,ML)
#}
return(list(A1,A2,A3,A4))
}


# rmarkdown::render(file.path(projectpath,"Rmd_files",'KPT.Rmd'),output_dir = output_html)
## ------------------------------------------------------------------------
FDPSIDWF<-function(i){
  if (WETNES[i] < WETINF[i]){ 
    FDPSIDW <- (-BEXP[i] * PSIF[i] / WETF[i]) * (WETNES[i] / WETF[i]) ^ (-BEXP[i] - 1)
  }else if (WETNES[i] < 1){
# in near-saturated range
    FDPSIDW <- CHM[i] * (2 * WETNES[i] - CHN[i] - 1)
  }else{
# saturated
    FDPSIDW <- 0
  }
  return(FDPSIDW)
}


## ------------------------------------------------------------------------
FPSIMF<-function(WETNESi, PSIFi, BEXPi, WETINFi, WETFi, CHMi, CHNi){
  if (WETNESi <= 0){ 
    FPSIM <- -10000000000
  }else if (WETNESi < WETINFi) {
    FPSIM <- PSIFi * (WETNESi / WETFi) ^ (-BEXPi)
  }else if (WETNESi < 1) {
# in near-saturated range
    FPSIM <- CHMi* (WETNESi - CHNi) * (WETNESi - 1)
  }else{
# saturated
  FPSIM <- 0
  }
  return(FPSIM)
}


## ------------------------------------------------------------------------
SOILPAR<-function(){
#local
  #Dim i%
  PSIINF<-rep(1,50)
#
  wetff<-WETF
  psigg<-PSIG
  thickk<-THICK
  SWATMx<-SWATMX
  PSIINf<-  PSIINF
  CHm<-CHM
  CHn<-CHN
  WETNEs<-WETNES
  SWATi<-SWATI
  KSAt<-KSAT
  WETc<-WETC
  
  for(i in 1:NLAYER){
    if (i == 1) {
      psigg[1] <- -RHOWG * thickk[1] / 2
    }else{
      psigg[i] <- psigg[i - 1] - RHOWG * ((thickk[i - 1] + thickk[i]) / 2)
    }
    SWATMx[i] <- thickk[i] * THSAT[i] * (1 - STONEF[i])
    wetff[i] <- THETAF[i] / THSAT[i]
    PSIINf[i] <- PSIF[i] * (WETINF[i] / wetff[i]) ^ -BEXP[i]
    CHm[i] <- (-PSIINf[i] / (1 - WETINF[i]) ^ 2) - BEXP[i] * (-PSIINf[i]) / (WETINF[i] * (1 - WETINF[i]))
    CHn[i] <- 2 * WETINF[i] - 1 - (-PSIINf[i] * BEXP[i] / (CHm[i] * WETINF[i]))
    if (PSIM[i] > 0) {
     # Stop
    }else if (PSIM[i] == 0) {
      WETNEs[i] <- 1
    }else{
      WETNEs[i]<- wetff[i] * (PSIM[i] / PSIF[i]) ^ (-1 / BEXP[i])
      if (WETNEs[i] > WETINF[i]){
        WETNEs[i] <- (1 + CHn[i]) / 2 + 0.5 * (CHn[i] ^ 2 - 2 * CHn[i] + 1 + 4 * PSIM[i] / CHm[i])^(1/2)
      }
    }
    SWATi[i] <- WETNEs[i] * SWATMx[i]
    KSAt[i] <- KF[i] * (1 / wetff[i]) ^ (2 * BEXP[i] + 3)
    WETc[i] <- wetff[i] * (1000 * PSICR / PSIF[i]) ^ (-1 / BEXP[i])
  }
return(list(PSICR, psigg, SWATMx, wetff, WETc, CHm, CHn, WETNEs, SWATi, KSAt))
}


## ------------------------------------------------------------------------
SOILVAR<-function(){
  PSITi<-PSITI
  THETa<-THETA
  Kk<-KK
  SWAt <- 0
  for (i in 1:NLAYER){
    PSITi[i] <- PSIM[i] + PSIG[i]
    THETa[i] <- WETNES[i] * THSAT[i]
    if(WETNES[i] > 0.0001){
      Kk[i] <- KF[i] * (WETNES[i] / WETF[i]) ^ (2 * BEXP[i] + 3)
    }else{
      Kk[i] <- 0.0000000001
    }
    SWAt <- SWAt + SWATI[i]
  }
return(c(PSITi, THETa, Kk, SWAt))
}


# rmarkdown::render(file.path(projectpath,"Rmd_files",'EVP.Rmd'),output_dir = output_html)
## ------------------------------------------------------------------------
INTER<-function(RFAL, PINT, LAI, SAI, FRINTL, FRINTS, CINTRL, CINTRS, DTP, INTR, RINT, IRVP){
#local
  INTRMX<-0
  CATCH<-0
  NEWINT<-0
  #
  CATCH <- (FRINTL * LAI + FRINTS * SAI) * RFAL
  INTRMX <- CINTRL * LAI + CINTRS * SAI
  NEWINT <- INTR + (CATCH - PINT) * DTP

  if(NEWINT > 0){
    IRVP <- PINT
    if(NEWINT > INTRMX){ 
      RINT <- PINT + (INTRMX - INTR) / DTP
    }else{
      RINT <- CATCH
    }
  }else{
    RINT <- CATCH
    IRVP <- (INTR / DTP) + CATCH
  }
  return(list(RINT,IRVP));
}


## ------------------------------------------------------------------------
INTER24<-function(RFAL, PINT, LAI, SAI, FRINTL, FRINTS, CINTRL, CINTRS, DURATN, INTR, RINT, IRVP, MONTHN){
#local
  INTRMX<-0 
  INTRNU<-0  
  NEWINT<-0
  RINTHR<-0 
  CATCH<-0   
  IRVPHR<-0 
  SMINT<-0    
  SMVP<-0     
  IHD<-0     
  #hh<<-0     
  DTH<-0
  #
  IHD <- as.integer((DURATN[MONTHN] + 0.1) / 2)
  INTRMX <- CINTRL * LAI + CINTRS * SAI
  INTRNU <- INTR
  SMINT <- 0
  SMVP <- 0
  DTH <- 1
  for(i in seq(0,23,1)){
    if((i < (12 - IHD)) || (i >= (12 + IHD))){
# before or after rain
      CATCH <- 0
    }else{
      CATCH <- (FRINTL * LAI + FRINTS * SAI) * RFAL / (2 * IHD)
    }
    NEWINT <- INTRNU + (CATCH - PINT / 24) * DTH
    if (NEWINT > 0.0001) {
      IRVPHR <- PINT / 24
      if (NEWINT > INTRMX) {
        RINTHR <- IRVPHR + (INTRMX - INTRNU) / DTH
      }else{
        RINTHR <- CATCH
      }
    }else{
        RINTHR <- CATCH
        IRVPHR <- INTRNU / DTH + CATCH
    }
  INTRNU <- INTRNU + (RINTHR - IRVPHR) * DTH
  SMVP <- SMVP + IRVPHR * DTH
  SMINT <- SMINT + RINTHR * DTH
  }
  IRVP <- SMVP
# / 1 d
  RINT <- SMINT
# / 1 d
  return(list(RINT,IRVP))
}


## ------------------------------------------------------------------------
PLNTRES<-function(NLAYER, THICK, STONEF, RTLEN, RELDEN, RTRAD, RPLANT, FXYLEM, RXYLEM, RROOTI, ALPHA){
#local
  #Dim I As Integer
  Dic<-c(seq(1,50,1))
  SUM<-0 
  RTFRAC<-0 
  RTDENI<-0 
  DELT<-0   
  RXYLEm<-0
  RROOTi<-rep(0,ML)
  ALPHa<-rep(0,ML)
  #
  RXYLEm <- FXYLEM * RPLANT
  for( i in seq( 1,NLAYER, 1)){
    Dic[i] <- THICK[i] * (1 - STONEF[i])
    SUM <- SUM + RELDEN[i] * Dic[i]
  }
  for( i in seq( 1,NLAYER,1)){
    if ((RELDEN[i] < 0.00001) || (RTLEN < 0.1)){
      RROOTi[i] <- 1E+20
      ALPHa[i] <- 1E+20
    }else{
      RTFRAC <- RELDEN[i] * Dic[i] / SUM
      RROOTi[i] <- (RPLANT - RXYLEm) / RTFRAC
      RTDENI <- RTFRAC * 0.001 * RTLEN / Dic[i]
# .001 is (mm/mm2)/(m/m2) conversion
      DELT <- PI * RTRAD ^ 2 * RTDENI
      ALPHa[i] <- (1 / (8 * PI * RTDENI)) * (DELT - 3 - 2 * (log(DELT)) / (1 - DELT))
      ALPHa[i] <- ALPHa[i] * 0.001 * RHOWG / Dic[i]
# .001 is MPa/kPa conversion
    }
  }
  return(c(RXYLEm,RROOTi,ALPHa))
}


## ------------------------------------------------------------------------
TBYLAYER<-function(J, PTR, DISPC, ALPHA, KK, RROOTI, RXYLEM, PSITI, NLAYER, PSICR, NOOUTF){
#local
  #Dim i%  
  RI<-rep(0,50) 
  RT<-0    
  SUM<-0  
  TRMIN<-0    
  PSIT<-0  
  R<-0     
  SUPPLY<-0         
  IDEL<-0   
  FLAG<-rep(0,50) 
  NEGFLAG<-0      
  ATr<-0
  ATRANi<-rep(0,ML)
  #
  for (i in seq( 1,NLAYER,1)){
    if (RROOTI[i] > 1E+15){
      FLAG[i] <- 1
    } else if((NOOUTF == 1) && (PSITI[i] / 1000 <= PSICR)) {
      FLAG[i] <- 1
    } else {
      FLAG[i]<- 0 # this layer has roots
    }
  }
  dfw<-0
# top of loop for recalculation of transpiration if more layers get flagged
 repeat{
    NEGFLAG <- 0
    SUM <- 0
    for(i in 1:NLAYER){
      if (FLAG[i]== 0){
        RI[i] <- RROOTI[i] + ALPHA[i] / KK[i]
        SUM <- SUM + 1 / RI[i]
      }else{
        ATRANi[i] <- 0
      }
    }
    if (SUM < 1E-20){
        ATr <- 0
        PSIT <- -10000000000
        return(list(ATr,ATRANi))
    }else{
        RT <- 1 / SUM
    }
# weighted mean soil water potential
    PSIT <- 0
    for (i in 1:NLAYER){
        if (FLAG[i] == 0){
          PSIT <- PSIT + RT * PSITI[i] / RI[i]
        }
    }
    SUPPLY <- (PSIT / 1000 - PSICR - RHOWG * DISPC) / (RT + RXYLEM)
    if (J == 1){
# daytime
      R <- (2 / PI) * (SUPPLY / PTR)
      if (R <= 0){ 
        ATr <- 0
      }else if (R < 1){
        ATr <- PTR * (1 + R * ACOSF(R) - sin(ACOSF(R)))
      }else{
        ATr <- PTR
      }
    }else{
# nighttime
      if ((SUPPLY <= 0) || (PTR <= 0)){
        ATr <- 0
      }else {
        ATr <- RMINF(SUPPLY, PTR)
      }
    }
    for (i in 1:NLAYER){
      if (FLAG[i] == 1){
        ATRANi[i] <- 0
      }else{
        ATRANi[i] <- ((PSITI[i]- PSIT) / 1000 + RT * ATr) / RI[i]
        if (ATRANi[i] < -0.000001){
          NEGFLAG <- 1
        }
      }  
    }
    dfw<-dfw+1
    if (NOOUTF == 1 && NEGFLAG == 1){ 
           IDEL <- 0
           TRMIN <- 0
           for(i in 1:NLAYER){
            if (ATRANi[i] < TRMIN) {
              TRMIN <- ATRANi[i]
              IDEL <- i
            }
           }
           FLAG[IDEL] <- 1
# repeat main loop with flagged layers excluded
    }else{
# done
       return(list(ATr,ATRANi))
    }
  }
}


# rmarkdown::render(file.path(projectpath,"Rmd_files",'PET.Rmd'),output_dir = output_html)
## ------------------------------------------------------------------------
CANOPY<-function(DOY, MAXHT, RELHT, MAXLAI, RELLAI, SNOW, SNODEN, MXRTLN, MXKPL, CS, DENSEF){
  #local
  SNODEP<-0 
  HNOSNO<-0 
  HSNO<-0   
  RATIO<-0 
  RELHIT<-0 
  KPL<-0
  #
  RELHIT <- INTERP(10, RELHT, DOY)
  SNODEP <- 0.001 * SNOW / SNODEN
  HNOSNO <- RMAXF(0.01, RELHIT * MAXHT)
  HSNO <- RMAXF(0, HNOSNO - SNODEP)
  RATIO <- HSNO / HNOSNO
  HEIGHT <- RMAXF(0.01, HSNO)
  #
  LAI <- RATIO * DENSEF * INTERP(10, RELLAI, DOY) * MAXLAI
  SAI <- DENSEF * CS * HEIGHT
  if(LAI < 0.00001)  LAI <- 0.00001
  #
  RTLEN <- DENSEF * RELHIT * MXRTLN
  KPL <- DENSEF * RELHIT * MXKPL
  if (KPL < 0.00000001)  KPL <- 0.00000001
  RPLANT <- 1 / KPL
  #
return(c(HEIGHT, LAI, SAI, RTLEN, RPLANT))
}


## ------------------------------------------------------------------------
ESAT<-function(TA, ES, DELTA){
  Es <- 0.61078 * exp(17.26939 * TA / (TA + 237.3))
  DELTa <- 4098 * Es / (TA + 237.3) ^ 2
  if (TA < 0) {
    Es <- 0.61078 * exp(21.87456 * TA / (TA + 265.5))
    DELTa <- 5808 * Es / (TA + 265.5) ^ 2
  }
return(c(Es, DELTa))
}


## ------------------------------------------------------------------------
FRSS<-function(RSSA, RSSB, PSIF, PSIM){
if (RSSA < 0.0001) {
  FRSs <- 10000000000
}else{
  FRSs <- RSSA * (PSIM / PSIF) ^ RSSB
}
return(FRSs)
}


## ------------------------------------------------------------------------
INTERP<-function(NPAIRS, FUNCT, XVALUE){
#local
#Dim I%, J%
 XX<-c(seq(1,10,1)) 
 YY<-c(seq(1,10,1))
# put FUNCT into XX and YY
  i <- 0
  for (J in seq(1,(2 * NPAIRS - 1),2)){
    i <- i + 1
    XX[i] <- FUNCT[J]
    YY[i] <- FUNCT[J + 1]
  }
# interpolate using XX and YY
  for (J in 1:NPAIRS){
    if (XVALUE == XX[J]){
      INTERp <- YY[J]
      return(INTERp)
    }else if (XVALUE < XX[J]){
      INTERp <- YY[J - 1] + (XVALUE - XX[J - 1]) * (YY[J] - YY[J - 1]) / (XX[J] - XX[J - 1])
      return(INTERp)
    }
  }
return(INTERp)
}


## ------------------------------------------------------------------------
PM<-function(AA, VPD, DELTA, RA, RC){
Pm <- (RA * DELTA * AA + CPRHO * VPD) / ((DELTA + GAMMA) * RA + GAMMA * RC)
return(Pm)
}


## ------------------------------------------------------------------------
ROUGH<-function(HEIGHT, ZMINH, LAI, SAI, CZS, CZR, HS, HR, LPC, CS, Z0GS){
#local
RATIO<-0
XX<-0

if (HEIGHT >= HR) {
  Z0C <- CZR * HEIGHT
}else if (HEIGHT <= HS){
  Z0C <- CZS * HEIGHT
}else{
  Z0C <- CZS * HS + (CZR * HR - CZS * HS) * (HEIGHT - HS) / (HR - HS)
}
  DISPC <- HEIGHT - Z0C / 0.3
  if (Z0GS > Z0C)  Z0GS <- Z0C
  RATIO <- (LAI + SAI) / (LPC + CS * HEIGHT)
  if (RATIO >= 1) {
    Z0 <- Z0C
    DISP <- DISPC
  }else{
    XX <- RATIO * (-1 + exp(0.909 - 3.03 * Z0C / HEIGHT)) ^ 4
    DISP <- 1.1 * HEIGHT * log(1 + XX ^ 0.25)
    Z0 <- RMINF(0.3 * (HEIGHT - DISP), Z0GS + 0.3 * HEIGHT * XX ^ 0.5)
  }
ZA <- HEIGHT + ZMINH
return(list(Z0GS, Z0C, DISPC, Z0, DISP, ZA))
}


## ------------------------------------------------------------------------
SRSC<-function(RAD, TA, VPD, LAI, SAI, GLMIN, GLMAX, R5, CVPD, RM, CR, TL, T1, T2, TH){
#local
FS<-0 
R0 <-0 
FRINT<-0 
FD <-0 
FT <-0    
GSC <-0  
#solar radiation limitation
FS <- (2 * LAI + SAI) / (2 * LAI)
if (RAD <= 0.0000000001){
    FRINT <- 0
}else{
  R0 <- RM * R5 / (RM - 2 * R5)
  FRINT <- ((RM + R0) / (RM * CR * FS)) * log((R0 + CR * RAD) / (R0 + CR * RAD * exp(-CR * FS * LAI)))
}
#vapor deficit limitation
FD <- 1 / (1 + VPD / CVPD)
#temperature limitation
if (TA <= TL) {
  FT <- 0
}else if (TA > TL && TA < T1) {
  FT <- 1 - ((T1 - TA) / (T1 - TL)) ^ 2
}else if (TA >= T1 && TA <= T2) {
  FT <- 1
}else if (TA > T2 && TA < TH) {
  FT <- 1 - ((TA - T2) / (TH - T2)) ^ 2
}else{
  FT <- 0
}
GSC <- FD * FT * FRINT * (GLMAX - GLMIN) + LAI * GLMIN
RSC <- 1 / GSC
return(RSC)
}


## ------------------------------------------------------------------------
SWGE<-function(AA, ASUBS, VPD, RAA, RAS, RSS, DELTA, ARATE, ERATE){
#local
RS<-0
RA<-0 
LE<-0  
LEC<-0 
#
LEC <- ARATE / (ETOM * WTOMJ)
RS <- (DELTA + GAMMA) * RAS + GAMMA * RSS
RA <- (DELTA + GAMMA) * RAA
LE <- (RS / (RS + RA)) * LEC + (CPRHO * VPD + DELTA * RAS * ASUBS + DELTA * RAA * AA) / (RS + RA)
ERATE <- ETOM * WTOMJ * (LE - LEC)
return(ERATE)
}


## ------------------------------------------------------------------------
SWGRA<-function(UA, ZA, HEIGHT, Z0, DISP, Z0C, DISPC, Z0GS, LWIDTH, RHOTP, NN, LAI, SAI, RAA, RAC, RAS){
#local
USTAR<-0
KH<-0
UH<-0
RB<-0
#
USTAR <- K * UA / (log((ZA - DISP) / Z0))
KH <- K * USTAR * (HEIGHT - DISP)
RAS <- (HEIGHT * exp(NN) / (NN * KH)) * (exp(-NN * Z0GS / HEIGHT) - exp(-NN * (Z0C + DISPC) / HEIGHT))
if (RAS < 1) RAS <- 1
RAA <- log((ZA - DISP) / (HEIGHT - DISP)) / (K * USTAR) + (HEIGHT / (NN * KH)) * (-1 + exp(NN * (HEIGHT - DISPC - Z0C) / HEIGHT))
UH <- (USTAR / K) * log((HEIGHT - DISP) / Z0)
RB <- (100 * NN) * (LWIDTH / UH) ^ 0.5 / (1 - exp(-NN / 2))
RAC <- RB / (RHOTP * LAI + PI * SAI)
return(list(RAA, RAC, RAS))
}


## ------------------------------------------------------------------------
SWPE<-function(AA, ASUBS, VPD, RAA, RAC, RAS, RSC, RSS, DELTA){
#local
RS<-0
RC<-0
RA<-0
PMS<-0
PMC<-0
D0<-0
CCS<-0
CCC<-0
LE<-0
#
RS <- (DELTA + GAMMA) * RAS + GAMMA * RSS
RC <- (DELTA + GAMMA) * RAC + GAMMA * RSC
RA <- (DELTA + GAMMA) * RAA
CCS <- 1 / (1 + RS * RA / (RC * (RS + RA)))
CCC <- 1 / (1 + RC * RA / (RS * (RC + RA)))
PMS <- PM(AA, VPD - DELTA * RAS * (AA - ASUBS) / CPRHO, DELTA, RAA + RAS, RSS)
PMC <- PM(AA, VPD - DELTA * RAC * ASUBS / CPRHO, DELTA, RAA + RAC, RSC)
LE <- (CCC * PMC + CCS * PMS)
D0 <- VPD + RAA * (DELTA * AA - (DELTA + GAMMA) * LE) / CPRHO
PRATE <- ETOM * WTOMJ * PM(AA - ASUBS, D0, DELTA, RAC, RSC)
ERATE <- ETOM * WTOMJ * PM(ASUBS, D0, DELTA, RAS, RSS)
return(list(PRATE, ERATE))
}


## ------------------------------------------------------------------------
WEATHER<-function(TMAX, TMIN, DAYLEN, I0HDAY, EA, UW, ZA, DISP, Z0, WNDRAT, FETCH, Z0W, ZW, SOLRAD, SOLRADC, TA, TADTM, TANTM, UA, UADTM, UANTM){
#local
dummy<-0
#
if (SOLRAD < 0.001) {
  SOLRADC <<- RRD * I0HDAY
}else if (SOLRAD > I0HDAY) {
  SOLRADC <<- 0.99 * I0HDAY
}else{
  SOLRADC <<- SOLRAD
}
TA <<- (TMAX + TMIN) / 2
TADTM <<- TA + ((TMAX - TMIN) / (2 * PI * DAYLEN)) * sin(PI * DAYLEN)
TANTM <<- TA - ((TMAX - TMIN) / (2 * PI * (1 - DAYLEN))) * sin(PI * DAYLEN)
if (EA == 0) {esat<-ESAT(TMIN, EA, dummy)
  EA<<-unlist(esat[1])
}
if (UW == 0) UW <<- UWD  #[28022018]  changed after Federer:be <<- UWD, where UWD is specified as 3.0 or some other value
if (UW < 0.2) UW <<- 0.2
if (Z0W < 0.000001) {
UA <<- UW
}else{
UA <<- UW * WNDADJ(ZA, DISP, Z0, FETCH, ZW, Z0W)
}
UADTM <<- UA / (DAYLEN + (1 - DAYLEN) * WNDRAT)
UANTM <<- WNDRAT * UADTM
return(list(SOLRADC, TA, TADTM, TANTM, UA, UADTM, UANTM))
}


## ------------------------------------------------------------------------
WNDADJ<-function(ZA, DISP, Z0, FETCH, ZW, Z0W){
#local
HIBL<-0
HIBL <- 0.334 * FETCH ^ 0.875 * Z0W ^ 0.125
WNDADj <- log(HIBL / Z0W) * log((ZA - DISP) / Z0) / (log(HIBL / Z0) * log(ZW / Z0W))
return(WNDADj)
}


# rmarkdown::render(file.path(projectpath,"Rmd_files",'SUN.Rmd'),output_dir = output_html)
## ------------------------------------------------------------------------
AVAILEN<-function (SLRAD, ALBEDO, C1, C2, C3, TA, EA, RATIO, SHEAT, CR, LAI, SAI){
#local
 SOLNET <-0
 EFFEM<-0   
 NOVERN<-0  
 CLDCOR<-0  
 LNGNET<-0   
 RN  <-0    
#
SOLNET <- (1 - ALBEDO) * SLRAD
EFFEM <- 1.24 * (EA * 10 / (TA + 273.15)) ^ (1 / 7)
NOVERN <- (RATIO - C1) / C2
if (NOVERN > 1)  NOVERN <- 1
if (NOVERN < 0)  NOVERN <- 0
CLDCOR <- C3 + (1 - C3) * NOVERN
# emissivity of the surface taken as 1.0 to also account for reflected
LNGNET <- (EFFEM - 1) * CLDCOR * SIGMA * (TA + 273.15) ^ 4
RN <- SOLNET + LNGNET
AA <- RN - SHEAT
ASUBS <- RN * exp(-CR * (LAI + SAI)) - SHEAT
return(list(RN, AA, ASUBS))
}


## ------------------------------------------------------------------------
EQUIVSLP<-function (LAT, SLOPE, ASPECT){
#Swift#s L1 and L2, Lee (3.31, 3.32)
#local
D1<-0
#
L1 <- ASINF(cos(SLOPE) * sin(LAT) + sin(SLOPE) * cos(LAT) * cos(ASPECT))
D1 <- cos(SLOPE) * cos(LAT) - sin(SLOPE) * sin(LAT) * cos(ASPECT)
if (D1 == 0)  D1 <- .0000000001
L2 <- atan(sin(SLOPE) * sin(ASPECT) / D1)
if (D1 < 0) L2 <- L2 + PI
return(list(L1, L2))
}



## ------------------------------------------------------------------------
SUNDS<-function (LAT, SLOPE, DOY, L1, L2, DAYLEN, I0HDAY, SLFDAY){
#local
I0SDAY <-0
SCD <-0  
DEC <-0 
TWORIS <-0 
Temp  <-0 
T0 <-0  
T1 <-0  
T2 <-0 
T3 <-0  
T6 <-0 
T7 <-0
T8 <-0  
T9 <-0   
#
SCD <- SC / (1 - .0167 * cos(.0172 * (DOY - 3))) ^ 2
DEC <- ASINF(.39785 * sin(4.868961 + .017203 * DOY + .033446 * sin(6.224111 + .017202 * DOY)))
Temp <- HAFDAY(LAT, DEC)
DAYLEN <- RMAXF(.0001, RMINF(.9999, Temp / PI))
# to avoid zero divides for 0 and 1
T1 <- Temp
T0 <- -Temp
Temp <- HAFDAY(L1, DEC)
T7 <- Temp - L2
T6 <- -Temp - L2
T3 <- RMINF(T1, T7)
T2 <- RMAXF(T0, T6)
  if (T3 < T2) {
    T2 <- 0
     T3 <- 0
  }
T6 <- T6 + 2 * PI
  if (T6 < T1) {
  T8 <- T6
  T9 <- T1
  TWORIS <- 1
  }
  T7 <- T7 - 2 * PI
  if (T7 > T0) {
    T8 <- T0
    T9 <- T7
    TWORIS <- 1
  }else{
  TWORIS <- 0
  }
if (TWORIS == 1) {   # two sunrises
  I0SDAY <- WTOMJ * SCD * (FUNC3(DEC, L2, L1, T3, T2) + FUNC3(DEC, L2, L1, T9, T8)) / cos(SLOPE)
# "daylength" on the slope = ((T3 - T2) + (T9 - T8)) / (2. * PI)
}else{    #  one sunrise
  I0SDAY <- WTOMJ * SCD * FUNC3(DEC, L2, L1, T3, T2) / cos(SLOPE)
# COS(SLOPE) adjusts from slope area to map area
# "daylength" on the slope = (T3 - T2) / (2. * PI)
}
  I0HDAY <- WTOMJ * SCD * FUNC3(DEC, 0, LAT, T1, T0)
  if (I0HDAY <= 0){
    SLFDAY <- 0
  }else{
    SLFDAY <- I0SDAY / I0HDAY
  }
  return(list(DAYLEN, I0HDAY, SLFDAY))
}


## ------------------------------------------------------------------------
FUNC3<-function(DEC, L2, L1, T3, T2){
#
FUnC3 <- (1 / (2 * 3.14159)) * (sin(DEC) * sin(L1) * (T3 - T2) + cos(DEC) * cos(L1) * (sin(T3 + L2) - sin(T2 + L2)))
return(FUnC3)
}


## ------------------------------------------------------------------------
HAFDAY<-function (LAT, DEC){
#local
ARG<-0
#
if (abs(LAT) >= PI / 2)  LAT <- sign(LAT) * (PI / 2 - 0.01)
ARG <- -tan(DEC) * tan(LAT)
if (ARG >= 1) {
HAFDAy <- 0
}else if (ARG <= -1) {
HAFDAy <- PI
}else{
HAFDAy <- ACOSF(ARG)
}
return(HAFDAy)
}


# rmarkdown::render(file.path(projectpath,"Rmd_files",'WAT.Rmd'),output_dir = output_html)
## ------------------------------------------------------------------------
BYFLFR<-function(){
#
for(i in 1:NLAYER){
  if (BYPAR == 1) {
    if (QFPAR > 0.01){
      BYFRAC[i] <- QFFC ^ (1 - (1 / QFPAR) * (WETNES[i] - WETF[i]) / (1 - WETF[i]))
      if (BYFRAC[i] > 1)  BYFRAC[i] <- 1
    }else{ # bucket for the layer
      if(WETNES[i] >= WETF[i]){
        BYFRAC[i] <- 1
      }else{
        BYFRAC[i] <- 0
      }
    }
  }else{
    BYFRAC[i] <- 0
  }
}
return(BYFRAC)
}


## ------------------------------------------------------------------------
DSLOP<-function(i){
    #local
    LL<-0 
    GRAD <-0
    ARATIO  <-0
    #
    LL <- 1000 * LENGTH
    GRAD <- RHOWG * sin(DSLOPE) + (2 * PSIM[i] / LL) * cos(DSLOPE)
    ARATIO <- THICK[i] * (1 - STONEF[i]) * cos(DSLOPE) / LL
    DSFLi <- KK[i] * ARATIO * GRAD / RHOWG
    if (DSFLi < 0)  DSFLi <- 0
return(DSFLi)
}


## ------------------------------------------------------------------------
GWATER<-function(GWAT, GSC, GSP, DT, VRFLN){
    if (GSC < 0.00000001){
      SEEP <- GSP * VRFLN
      GWFL <- VRFLN - SEEP
    }else{
      SEEP <- GWAT * GSC * GSP
      GWFL <- GWAT * GSC * (1 - GSP)
      if (GWAT / DT - (GWFL + SEEP) < 0){
        SEEP <- GSP * GWAT / DT
        GWFL <- (1 - GSP) * GWAT / DT
      }
    }
return(list(GWFL, SEEP))
}


## ------------------------------------------------------------------------
INFLOW<-function(){
    #local
    #Dim i% 
    INFIL<-0 
    MAXIN<-0
    INFLi<-INFLI
    #
    for(i in seq(NLAYER,1,-1)){
      INFIL <- SLFL * INFRAC[i]
      BYFLI[i] <- BYFRAC[i] * INFIL
      INFLi[i] <- INFIL - BYFLI[i]
      if (i == NLAYER)
        {VV[i] <- VRFLI[i]
           }
      if (i > 1) {
       MAXIN <- (SWATMX[i] - SWATI[i]) / DTI + VV[i] + DSFLI[i] + TRANI[i]
      
        if (VRFLI[i - 1] + INFLi[i] > MAXIN) {
          if (BYFRAC[i] > 0) {
            if (VRFLI[i - 1] < MAXIN) {
              BYFLI[i] <- BYFLI[i] + INFLi[i] - (MAXIN - VRFLI[i - 1])
              INFLi[i] <- MAXIN - VRFLI[i - 1]
              VV[i-1] <- VRFLI[i-1]
            }else{
              BYFLI[i] <- BYFLI[i] + INFLi[i]
              INFLi[i] <- 0
              VV[i-1] <- MAXIN
            }
          }else{
            VV[i-1] <- MAXIN - INFLi[i]
          }
        }else{
          VV[i-1] <- VRFLI[i-1]
        }
       NTFLI[i] <- VV[i-1] + INFLi[i] - VV[i] - DSFLI[i] - TRANI[i]

      }else{
      # i% = 1
        MAXIN <- (SWATMX[1] - SWATI[1]) / DTI + VV[1] + DSFLI[1] + TRANI[1] + SLVP
        if (INFLi[1] > MAXIN) {
          BYFLI[1] <- BYFLI[1] + INFLi[1] - MAXIN
          INFLi[1] <- MAXIN
        }
        NTFLI[1] <- INFLi[1] - VV[1] - DSFLI[1] - TRANI[1] - SLVP
      }
    }
   INFLI<- INFLi
return(list(VV, INFLI, BYFLI, NTFLI))
}


## ------------------------------------------------------------------------
INFPAR<-function(INFEXP, IDEPTH, NLAYER, THICK){
#local
THICKT<-0 
THICKA<-rep(0,ML)
#
if (INFEXP <= 0 || IDEPTH == 0) {
  ILAYER <- 1  # probably not used
  INFRAC[1] <- 1
  for (i in seq(2,NLAYER,1)){
    INFRAC[i] <- 0
  }
}else{
# must have at least one layer
  THICKT <- THICK[1]
  ILAYER <- 1
  for (i in 2:NLAYER){
    if (THICKT + 0.5 * THICK[i] <= IDEPTH){
      ILAYER <- ILAYER + 1
      THICKT <- THICKT + THICK[i]
    }else{
      i<-NLAYER
    }
  }
  THICKA[1] <- 0
  for(i in 1:NLAYER){
    if (i <= ILAYER) {
      
      if(i==1){
        THICKA[i] <- THICK[i]
        INFRAC[i] <- (THICKA[i] / THICKT) ^ INFEXP - (0 / THICKT) ^ INFEXP
      }else{
        THICKA[i] <- THICKA[i - 1] + THICK[i]
        INFRAC[i] <- (THICKA[i] / THICKT) ^ INFEXP - (THICKA[i - 1] / THICKT) ^ INFEXP
      }
    }else{
      INFRAC[i] <- 0
    }
  }  
}
return(list(ILAYER, INFRAC))
}



## ------------------------------------------------------------------------
ITER<-function(NLAYER, DTI, DPSIDW, NTFLI, SWATMX, PSITI, DSWMAX, DPSIMX){
#local
A<- rep(0,50)
temp<-rep(0,50)
TT  <-0 
PP<-0
#
for ( i in 1:NLAYER){
    A[i] <- NTFLI[i] * DPSIDW[i] / SWATMX[i]
    temp[i] <- PSITI[i] + A[i] * DTI
}
DTINEW <- DTI
for( i in 1:NLAYER){
  DTINEW <- RMINF(DTINEW, 0.01 * DSWMAX * SWATMX[i] / RMAXF(0.000001, abs(NTFLI[i])))
  if (i < NLAYER) {
    PP <- PSITI[i] - PSITI[i + 1]
    TT <- temp[i] - temp[i + 1]
    if ((abs(TT) > DPSIMX) && (abs(PP) > DPSIMX) && (sign(TT) != sign(PP))){
      DTINEW <- RMINF(DTINEW, -PP / (A[i] - A[i + 1]))
    }
  }
}
return(DTINEW)
}


## ------------------------------------------------------------------------
RTDEN<-function(ROOTDEN, NLAYER, THICK){
  
RTHICK<-rep(0,ML) 
RDEN<-rep(0,ML)
RREMAIN<-0
TREMAIN<-0

for( J in 1:ML){
  RTHICK[J] <- ROOTDEN[2 * J - 1]
  RDEN[J] <- ROOTDEN[2 * J]
}
DONE <- FALSE
j <- 1
RREMAIN <- RTHICK[j]
for( i in 1:NLAYER){ # new soil layer
# accumulate RELDEN as total root length in soil layer
  RELDEN[i] <- 0
  if (!DONE){
    TREMAIN <- THICK[i]
    while(RREMAIN < TREMAIN && j < ML-1){
# remaining root layer thickness < remaining soil layer thickness
      RELDEN[i] <- RELDEN[i] + RDEN[j] * RREMAIN
      TREMAIN <- TREMAIN - RREMAIN
      j <- j + 1
      if( j == ML){
        DONE <- TRUE
      }
      RREMAIN <- RTHICK[j]
    }
# remaining root layer thickness >= remaining soil layer thickness
    if(!DONE){
       RELDEN[i] <- RELDEN[i] + RDEN[j] * TREMAIN
       RREMAIN <- RREMAIN - TREMAIN
    }
  }
# convert back to unit volume basis
RELDEN[i] <- RELDEN[i] / THICK[i]
}
return(RELDEN)
}


## ------------------------------------------------------------------------
SRFLFR<-function(){
#local
SUM<-0
safra<-0
#
for( i in  1:QLAYER){
  SUM <- SUM + SWATI[i]
}
if( QFPAR > 0.01){
  safra <- QFFC ^ (1 - (1 / QFPAR) * (SUM - SWATQF) / (SWATQX - SWATQF))
  if (safra > 1) {safra <- 1}
}else{ # bucket over QLAYERs
  if (SUM >= SWATQF){
      safra <- 1
  }else{
      safra <- 0
  }
}
return(safra)
}


## ------------------------------------------------------------------------
SRFPAR<-function(QDEPTH, NLAYER, THETAF, THICK, STONEF, SWATMX){
#local
THICKT<-0
#
if( QDEPTH == 0 ){
  QLAYER <- 0
  SWATQX <- 0
  SWATQF <- 0
  return(list( QLAYER, SWATQX, SWATQF))
}
QLAYER <- 1
THICKT <- THICK[1]
for( i in  2:NLAYER){
  if( THICKT + 0.5 * THICK[i] <= QDEPTH){
    THICKT <- THICKT + THICK[i]
    QLAYER <- QLAYER+ 1
  }else{
    i<-NLAYER
  }
}
  SWATQX <- 0
  SWATQF <- 0
  for( i in  1:QLAYER){
    SWATQX <- SWATQX + SWATMX[i]
    SWATQF <- SWATQF + THETAF[i] * THICK[i] * (1 - STONEF[i])
  }
return(list( QLAYER, SWATQX, SWATQF))
}


## ------------------------------------------------------------------------
VERT<-function(i){
#local
GRAD <-0
KKMEAN<-0
#
KKMEAN <- exp((log(KK[i]) + log(KK[i+1])) / 2)
if (KKMEAN > KSAT[i])  KKMEAN <- KSAT[i]
if (KKMEAN > KSAT[i+1])  KKMEAN <- KSAT[i+1]
GRAD <- (PSITI[i] - PSITI[i+1]) / RMINF(THICK[i], THICK[i+1])
VRFLi <- (GRAD * KKMEAN / RHOWG) * (1 - (STONEF[i] + STONEF[i+1]) / 2)

return(VRFLi)
}


# rmarkdown::render(file.path(projectpath,"Rmd_files",'SNO.Rmd'),output_dir = output_html)
## ------------------------------------------------------------------------
SNOENRGY<-function(TSNOW, TA, DAYLEN, CCFAC, MELFAC, SLFDAY, LAI, SAI, LAIMLT, SAIMLT){
if (TA <= 0) {
  SNOEN <- CCFAC * 2 * DAYLEN * (TA - TSNOW)
}else{
  SNOEN <- MELFAC * 2 * DAYLEN * TA * exp(-SAIMLT * SAI) * exp(-LAIMLT * LAI) * SLFDAY
}
return(SNOEN)
}


## ------------------------------------------------------------------------
SNOFRAC<-function (TMAX, TMIN, RSTEMP){
if (TMIN >= RSTEMP) {
  SNOFRC <- 0
}else if (TMAX < RSTEMP){
  SNOFRC <- 1
}else{
  SNOFRC <- 1 - (TMAX - RSTEMP) / (TMAX - TMIN)
}
return(SNOFRC)
}


## ------------------------------------------------------------------------
SNOVAP<-function (TSNOW, TA, EA, UA, ZA, HEIGHT, Z0, DISP, Z0C, DISPC, Z0GS, LWIDTH, RHOTP, NN, LAI, SAI, KSNVP){
#local
ESNOW <-0 
RAA   <-0 
RAS   <-0  
# ignores effect of interception on PSNVP or of PSNVP on PTRAN
if (TSNOW > -0.1) {
  ESNOW <- 0.61
}else{
# snow surface vapor pressure saturated at lower of TA and TSNOW
  esatt<-ESAT(RMINF(TA, TSNOW), ESNOW, dummy)
   ESNOW<-unlist(esatt[1])
}
swgra<-SWGRA(UA, ZA, HEIGHT, Z0, DISP, Z0C, DISPC, Z0GS, LWIDTH, RHOTP, NN, LAI, SAI, RAA, RAC, RAS)
  RAA<-unlist(swgra[1])
  RAC<-unlist(swgra[2])
  RAS<-unlist(swgra[3])

PSNVP <- (WTOMJ / LS) * (CPRHO / GAMMA) * (ESNOW - EA) / (RAA + RAS)
PSNVP <- KSNVP * PSNVP
return(PSNVP)
}


## ------------------------------------------------------------------------
SNOWPACK<-function(RTHR, STHR, PSNVP, SNOEN, CC, SNOW, SNOWLq, DTP, TA, MAXLQF, GRDMLT){
#local
SNOWLQ<-SNOWLq
FRAC <-0
EQEN <-0
NMLT<-0
ALQ <-0
RIN <-0 
#
# snow throughfall and its cold content, SNOWLQ unchanged
SNOW <- SNOW + STHR * DTP
CC <- CC + CVICE * RMAXF(-TA, 0) * STHR * DTP

if (CC > 0 && SNOWLQ > 0) {
  if (CC > SNOWLQ * LF) {
# refreeze all liquid
    CC <- CC - SNOWLQ * LF
    SNOWLQ <- 0
  }else{
# refreeze part
    SNOWLQ <- SNOWLQ - CC / LF
    CC <- 0
  }
}
# groundmelt and evaporation loss as fraction of SNOW
FRAC <- (GRDMLT + PSNVP) * DTP / SNOW
# FRAC can be negative if condensation exceeds groundmelt
if (FRAC < 1) {
  SMLT <- GRDMLT
  SNVP <- PSNVP
# reduce CC, SNOWLQ, and SNOW proportionally for groundmelt and evaporation
# increase them proportionally if condensation exceeds groundmelt
  CC <- CC * (1 - FRAC)
  SNOWLQ <- SNOWLQ * (1 - FRAC)
  SNOW <- SNOW * (1 - FRAC)
}else{
# all SNOW disappears from groundmelt and/or evaporation
  SMLT <- GRDMLT / FRAC
  SNVP <- PSNVP / FRAC
  RSNO <- 0
  CC <- 0
  SNOWLQ <- 0
  SNOW <- 0
}
# snowpack cooling or warming
if (SNOW > 0) {
# equivalent ice melted by energy input including warm rain (mm)
  EQEN <- DTP * (SNOEN + RTHR * RMAXF(TA, 0) * CVLQ) / LF
  if (EQEN <= 0) {
# snowpack cooling
    NMLT <- -EQEN
    if (NMLT < SNOWLQ) {
# only part of SNOWLQ refreezes
      CC <- 0
# should be 0 already because SNOWLQ is positive
      SNOWLQ <- SNOWLQ - NMLT
    }else{
# all SNOWLQ (if any) refreezes, remaining NMLT increases CC
      NMLT <- NMLT - SNOWLQ
      SNOWLQ <- 0
      CC <- CC + NMLT * LF
# do not allow TSNOW to cool below TA
      CC <- RMINF(CC, -TA * SNOW * CVICE)
    }
  }else{
# snowpack warming  (cant have both CC and SNOWLQ)
    if (EQEN * LF < CC || TA < 0) {
# reduce but dont eliminate CC
      if (TA < 0) {
# do not allow TSNOW to warm above TA when TA < 0
        CC <- RMAXF(CC - EQEN * LF, -TA * SNOW * CVICE)
      }else{
        CC <- CC - EQEN * LF
      }
      SNOWLQ <- 0
    }else{
# CC eliminated
      EQEN <- EQEN - CC / LF
      CC <- 0
      if (EQEN <= MAXLQF * SNOW - SNOWLQ){
# remaining energy increases liquid water
        SNOWLQ <- SNOWLQ + EQEN
# SMLT and SNOW unchanged
      }else{
# liquid water capacity reached, SNOW melt produced
        EQEN <- EQEN - (MAXLQF * SNOW - SNOWLQ)
        if (SNOW * (1 - MAXLQF) > EQEN) {
# melt is ice plus the liquid included in it
          SMLT <- SMLT + (EQEN / DTP) / (1 - MAXLQF)
          SNOW <- SNOW - EQEN / (1 - MAXLQF)
          SNOWLQ <- MAXLQF * SNOW
        }else{
# all SNOW melts
          SMLT <- SMLT + SNOW / DTP
          SNOW <- 0
          SNOWLQ <- 0
        }
      }
    }
  }
# add rain to snowpack,
if (RTHR == 0 || SNOW == 0) {
  RSNO <- 0
}else{
# rain on SNOW
  RIN <- RTHR * DTP
  if (CC > 0) {
# use CC to refreeze rain
      if (CC > RIN * LF) {
# refreezes all rain
        CC <- CC - RIN * LF
        RSNO <- RTHR
        SNOW <- SNOW + RIN
      }else{
# CC refreezes part of rain
        SNOW <- SNOW + CC / LF
        RSNO <- (CC / LF) / DTP
        CC <- 0
# remaining rain
        RIN <- RIN - RSNO * DTP
# increase liquid water, SNOWLQ initially zero
        if (RIN < MAXLQF * SNOW / (1 - MAXLQF)) {
# remaining RIN all to SNOWLQ
          SNOWLQ <- RIN
          RSNO <- RSNO + RIN / DTP
          SNOW <- SNOW + RIN
        }else{
          SNOWLQ <- MAXLQF * SNOW / (1 - MAXLQF)
          RSNO <- RSNO + SNOWLQ / DTP
          SNOW <- SNOW + SNOWLQ
        }
      }
    }else{
# CC = 0.
        if (SNOWLQ >= MAXLQF * SNOW) {
# SNOW already holding maximum liquid
          RSNO <- 0
        }else{
          ALQ <- MAXLQF * SNOW - SNOWLQ
          if (RIN < ALQ) {
# all RIN to SNOW
            RSNO <- RTHR
            SNOWLQ <- SNOWLQ + RIN
            SNOW <- SNOW + RIN
          }else{
# maximum liquid reached
            RSNO <- (ALQ / (1 - MAXLQF)) / DTP
            SNOW <- SNOW + RSNO * DTP
            SNOWLQ <- MAXLQF * SNOW
          }
        }
      }
    }
}    
return (list(CC,SNOW,SNOWLQ,RSNO, SNVP, SMLT))
}


# rmarkdown::render(file.path(projectpath,"Rmd_files",'B90V4_sub.Rmd'),output_dir = output_html)
## ------------------------------------------------------------------------
subdatafileline<-function(row){

YY<<- MData[[1]][row]
if( YY < 100){
  if (YY > 20){
    YY <<- YY + 1900
  }else{
    YY <<- YY + 2000
  }
}
MM<<- MData[[2]][row]
DD<<- MData[[3]][row]
SOLRAD<<- MData[[4]][row]
TMAX <<- MData[[5]][row]
TMIN <<- MData[[6]][row]
EA <<- MData[[7]][row]
UW <<- MData[[8]][row]
PRECIN <<- MData[[9]][row]
MESFL <<- MData[[10]][row]
}


## ------------------------------------------------------------------------
subprfileline<-function(row){
YY <<- MhhData[[1]][row]
if( YY < 100){
  if( YY > 20 ){
YY <<- YY + 1900
}else{
YY <<- YY + 2000
}
}
MM <<- MhhData[[2]][row]
DD <<- MhhData[[3]][row]
II <<- MhhData[[4]][row]
PREINT <<- MhhData[[5]][row]
MESFLP <<- MhhData[[6]][row]
}


## ------------------------------------------------------------------------
subprfilelineStatPrecip<-function(row,precip){
  YY <<- MhhData[[1]][row]
   if( YY < 100){
    if( YY > 20 ){
      YY <<- YY + 1900
    }else{
      YY <<- YY + 2000
    }
  }
  MM <<- MhhData[[2]][row]
  DD <<- MhhData[[3]][row]
  II <<- MhhData[[4]][row]
  PREINT <<- precip
  MESFLP <<- MhhData[[6]][row]
}


## ------------------------------------------------------------------------
MSBSETVARS<-function(){
# solar parameters depending on DOY%
sundss<-SUNDS(LAT, ESLOPE, DOY, L1, L2)
  
  DAYLEN<<-unlist(sundss[1])
  I0HDAY<<-unlist(sundss[2])
  SLFDAY<<-unlist(sundss[3])
# canopy parameters depending on DOY%
cano<-CANOPY(DOY, MAXHT, RELHT, MAXLAI, RELLAI, SNOW, SNODEN, MXRTLN, MXKPL, CS, DENSEF)
    HEIGHT<<-unlist(cano[1])
    LAI<<-unlist(cano[2])
    SAI<<-unlist(cano[3])
    RTLEN<<-unlist(cano[4])
    RPLANT<<-unlist(cano[5])
# roughness parameters
if (SNOW > 0) {
  Z0GS <<- Z0S
}else{
  Z0GS <<- Z0G
}
rough<-ROUGH(HEIGHT, ZMINH, LAI, SAI, CZS, CZR, HS, HR, LPC, CS, Z0GS)
  Z0GS<<-unlist(rough[1])
  Z0C<<-unlist(rough[2])
  DISPC<<-unlist(rough[3])
  Z0<<-unlist(rough[4])
  DISP<<-unlist(rough[5])
  ZA<<-unlist(rough[6])
# plant resistance components
plnt<-PLNTRES(NLAYER, THICK, STONEF, RTLEN, RELDEN, RTRAD, RPLANT, FXYLEM)
   RXYLEM<<-plnt[1]
   RROOTI<<-plnt[2:(ML+1)]
   ALPHA<<-plnt[(ML+2):(ML*2+1)]
# calculated weather data
SHEAT <<- 0
WEATHER(TMAX, TMIN, DAYLEN, I0HDAY, EA, UW, ZA, DISP, Z0, WNDRAT, FETCH, Z0W, ZW, SOLRAD, SOLRADC, TA, TADTM, TANTM, UA, UADTM, UANTM)
# fraction of precipitation as SFAL
SNOFRC<<- SNOFRAC(TMAX, TMIN, RSTEMP)
if (SNOW > 0) {
# snowpack temperature at beginning of day
  TSNOW <<- -CC / (CVICE * SNOW)
# potential snow evaporation
  PSNVP<<-SNOVAP(TSNOW, TA, EA, UA, ZA, HEIGHT, Z0, DISP, Z0C, DISPC, Z0GS, LWIDTH, RHOTP, NN, LAI, SAI, KSNVP)
  ALBEDO <<- ALBSN
  RSS <<- 0
}else{
  TSNOW <<- 0
  PSNVP <<- 0
  ALBEDO <<- ALB
# soil evaporation resistance
  RSS <<- FRSS(RSSA, RSSB, PSIF[1], PSIM[1])
# check for zero or negative RSS
  if (RSS < 0.000001) {
# MsgBox ("RSS is very small or negative. Run ends. Check RSSA and RSSB values.")
    rstop <<- 3
  }
}
# snow surface energy balance
SNOEN<<-SNOENRGY(TSNOW, TA, DAYLEN, CCFAC, MELFAC, SLFDAY, LAI, SAI, LAIMLT, SAIMLT)
}


## ------------------------------------------------------------------------
MSBPREINT<-function(){
 PREC <<- PREINT / DTP
 SFAL <<- SNOFRC * PREC
 RFAL <<- PREC - SFAL
if (NPINT > 1) {
# more than one precip interval in day
# snow interception
  if (PINT < 0 && TA > 0) {
    # prevent frost when too warm, carry negative PINT to rain
    temppp<-INTER(SFAL, 0, LAI, SAI, FSINTL, FSINTS, CINTSL, CINTSS, DTP, INTS, SINT, ISVP)
    SINT<<-unlist(temppp[1])
    ISVP<<-unlist(temppp[2])
  }else{
    temppp<-INTER(SFAL, PINT, LAI, SAI, FSINTL, FSINTS, CINTSL, CINTSS, DTP, INTS, SINT, ISVP)
      SINT<<-unlist(temppp[1])
      ISVP<<-unlist(temppp[2])
  }
# rain interception,  note potential interception rate is PID/DT-ISVP
  temppp<-INTER(RFAL, PINT - ISVP, LAI, SAI, FRINTL, FRINTS, CINTRL, CINTRS, DTP, INTR, RINT, IRVP)
    RINT<<-unlist(temppp[1])
    IRVP<<-unlist(temppp[2])
}else{
# one precip interval in day, use storm DURATN and INTER24
# snow interception
  if (PINT < 0 && TA > 0) {
  # prevent frost when too warm, carry negative PINT to rain
    temm<-INTER24(SFAL, 0, LAI, SAI, FSINTL, FSINTS, CINTSL, CINTSS, DURATN, INTS, SINT, ISVP, MONTHN)
      SINT<<-unlist(temm[1])
      ISVP<<-unlist(temm[2])
  }else{
    temm<-INTER24(SFAL, PINT, LAI, SAI, FSINTL, FSINTS, CINTSL, CINTSS, DURATN, INTS, SINT, ISVP, MONTHN)
      SINT<<-unlist(temm[1])
      ISVP<<-unlist(temm[2])
  }
# rain interception,  note potential interception rate is PID/DT-ISVP
   temm<-INTER24(RFAL, PINT - ISVP, LAI, SAI, FRINTL, FRINTS, CINTRL, CINTRS, DURATN, INTR, RINT, IRVP, MONTHN)
     RINT<<-unlist(temm[1])
     IRVP<<-unlist(temm[2])
}
# throughfall
RTHR <<- RFAL - RINT
STHR <<- SFAL - SINT
# reduce transpiration for fraction of precip interval that canopy is wet
WETFR <<- RMINF(1, (IRVP + ISVP) / PINT)
PTRAN <<- (1 - WETFR) * PTRAN
for( i in 1:NLAYER){
  TRANI[i] <<- (1 - WETFR) * TRANI[i]
}
if (SNOW <= 0 && STHR <= 0) {
# no snow, soil evaporation weighted for WETFR
  SLVP <<- WETFR * GIVP + (1 - WETFR) * GEVP
  RNET <<- RTHR
  RSNO <<- 0
  SNVP <<- 0
  SMLT <<- 0
}else{
  if (SNOW <= 0 && STHR > 0){
    CC <<- 0
    SNOWLQ <<- 0
  }
# snow accumulation and melt
  spa<-SNOWPACK(RTHR, STHR, PSNVP, SNOEN, CC, SNOW, SNOWLQ, DTP, TA, MAXLQF, GRDMLT)
    CC      <<-unlist(spa[1])
    SNOW    <<-unlist(spa[2])
    SNOWLQ  <<-unlist(spa[3])
    RSNO    <<-unlist(spa[4])
    SNVP    <<-unlist(spa[5])
    SMLT    <<-unlist(spa[6])
  RNET <<- RTHR - RSNO
  SLVP <<- 0
}
}


## ------------------------------------------------------------------------
MSBITERATE<-function(){
# source area flow rate
if (QLAYER > 0) {
  SAFRAC<<-SRFLFR()
}else{
  SAFRAC <<- 0
}
SRFL <<- RMINF(1, (IMPERV + SAFRAC)) * (RNET + SMLT)
# water supply rate to soil surface
SLFL <<- RNET + SMLT - SRFL
# bypass fraction of infiltration to each layer
BYFRAC<<-BYFLFR()
#
for( i in  seq(NLAYER,1,-1)){
    # downslope flow rates
    if( LENGTH == 0 || DSLOPE == 0){  
    # added in Version 4
      DSFLI[i]<<- 0
    }else{
      DSFLI[i]<<-DSLOP(i)
    }
    # vertical flow rates
    if (i < NLAYER) {
      if (abs(PSITI[i] - PSITI[i+1]) < DPSIMX) {
        VRFLI[i] <<- 0
       
      }else{
        VRFLI[i]<<-VERT(i)
      }
    }else{
    # bottom layer
      if( DRAIN > 0.0001){
      # gravity drainage only
        VRFLI[NLAYER] <<- DRAIN * KK[NLAYER] * (1 - STONEF[NLAYER])
      }else{
      # bottom of profile sealed
        VRFLI[NLAYER] <<- 0
      }
    }  
    if (IDAY >= 6 && i==NLAYER) {
      DRAIN<-DRAIN
    }
}
DTI <<- RMINF(DTRI, DTIMAX)
inflo<-INFLOW()
   VV<<-unlist(inflo[1]) 
   INFLI<<-unlist(inflo[2])
   BYFLI<<-unlist(inflo[3])
   NTFLI<<-unlist(inflo[4])
for( i in 1:NLAYER){
  DPSIDW[i] <<- FDPSIDWF(i)
}
DTINEW<<-ITER(NLAYER, DTI, DPSIDW, NTFLI, SWATMX, PSITI, DSWMAX, DPSIMX)
  if (DTINEW < DTI) {
    # recalculate flow rates with new DTI
    if (mnuhalfiter == FALSE) {
      DTI <<- DTINEW
    }else{
      DTI <<- 0.5 * DTINEW
    }
    inflo<-INFLOW()
       VV<<-unlist(inflo[1]) 
       INFLI<<-unlist(inflo[2])
       BYFLI<<-unlist(inflo[3])
       NTFLI<<-unlist(inflo[4])
  }
for( i in 1:NLAYER){
  VRFLI[i] <<- VV[i]
}
# groundwater flow and seepage loss
gwa<-GWATER(GWAT, GSC, GSP, DT, VRFLI[NLAYER])
    GWFL<<-unlist(gwa[1])
    SEEP<<-unlist(gwa[2])
}


## ------------------------------------------------------------------------
MSBDAYNIGHT<-function(){
SOVERI<<-0
for( J in  1:2){
# net radiation
  if (J ==1){
    SLRAD <<- SLFDAY * SOLRADC / (WTOMJ * DAYLEN)
    SLRADd<<-SLRAD
    TAJ <<- TADTM
    UAJ <<- UADTM
  }else{
    SLRAD <<- 0
    TAJ <<- TANTM
    UAJ <<- UANTM
  }
  if (I0HDAY <= 0.01){
  # no sunrise, assume 50% clouds for longwave
    SOVERI <<- 0.5
  }else{
    SOVERI <<- SOLRADC / I0HDAY
  }
  avai<-AVAILEN(SLRAD, ALBEDO, C1, C2, C3, TAJ, EA, SOVERI, SHEAT, CR, LAI, SAI)
    AA<<-unlist(avai[2])
    ASUBS<<-unlist(avai[3])
# vapor pressure deficit
  esat<-ESAT(TAJ, ES, DELTA)
    ES<<-unlist(esat[1])
    DELTA<<-unlist(esat[2])
    VPD <<- ES - EA
# S.-W. resistances
  swgra<-SWGRA(UAJ, ZA, HEIGHT, Z0, DISP, Z0C, DISPC, Z0GS, LWIDTH, RHOTP, NN, LAI, SAI, RAA, RAC, RAS)
    RAA<<-unlist(swgra[1])
    RAC<<-unlist(swgra[2])
    RAS<<-unlist(swgra[3])
  if (J == 1) {
    RSC<<-SRSC(SLRAD, TA, VPD, LAI, SAI, GLMIN, GLMAX, R5, CVPD, RM, CR, TL, T1, T2, TH)
  }else{
    RSC <<- 1 / (GLMIN * LAI)
  }
# S.-W. potential transpiration and ground evaporation rates
  swpe<-  SWPE(AA, ASUBS, VPD, RAA, RAC, RAS, RSC, RSS, DELTA)
    PTR[J]<<-unlist(swpe[1])
    GER[J]<<-unlist(swpe[2])
# S.-W. potential interception and ground evap. rates
# RSC = 0, RSS not changed
  swpe<-  SWPE(AA, ASUBS, VPD, RAA, RAC, RAS, 0, RSS, DELTA)
    PIR[J]<<-unlist(swpe[1])
    GIR[J]<<-unlist(swpe[2])
# actual transpiration and ground evaporation rates
  if (PTR[J] > 0.001) {
      rbl<-TBYLAYER(J, PTR[J], DISPC, ALPHA, KK, RROOTI, RXYLEM, PSITI, NLAYER, PSICR, NOOUTF)
      ATR[J]<<-unlist(rbl[1])
      ATRANI<<-unlist(rbl[2])
    for (i in 1:NLAYER){
      ATRI[J,i] <<- ATRANI[i]
    }
    if (ATR[J] < PTR[J]){
    # soil water limitation
      GER[J]<<-SWGE(AA, ASUBS, VPD, RAA, RAS, RSS, DELTA, ATR[J], GER[J])
    }
  }else{
    # no transpiration, condensation ignored
    PTR[J] <<- 0
    ATR[J] <<- 0
    for( i in 1:NLAYER){
      ATRI[J,i] <<- 0
    }
    GER[J]<<-SWGE(AA, ASUBS, VPD, RAA, RAS, RSS, DELTA, 0, GER[J])
  }
}
}


## ------------------------------------------------------------------------
psum<-function(){
 EVAPP <<- (ISVP + IRVP + SNVP + SLVP) * DTP + TRANP
 FLOWP <<- SRFLP + BYFLP + DSFLP + GWFLP
}


## ------------------------------------------------------------------------
ysum<-function(){
PRECY <<- RFALY + SFALY
STHRY <<- SFALY - SINTY
RTHRY <<- RFALY - RINTY
RNETY <<- RTHRY - RSNOY
EVAPY <<- IRVPY + ISVPY + SNVPY + SLVPY + TRANY
FLOWY <<- SRFLY + BYFLY + DSFLY + GWFLY
}


## ------------------------------------------------------------------------
msum<-function(){
PRECM <<- RFALM + SFALM
STHRM <<- SFALM - SINTM
RTHRM <<- RFALM - RINTM
RNETM <<- RTHRM - RSNOM
EVAPM <<- IRVPM + ISVPM + SNVPM + SLVPM + TRANM
FLOWM <<- SRFLM + BYFLM + DSFLM + GWFLM
}


## ------------------------------------------------------------------------
dsum<-function(){
PRECD <<- RFALD + SFALD
STHRD <<- SFALD - SINTD
RTHRD <<- RFALD - RINTD
RNETD <<- RTHRD - RSNOD
EVAPD <<- IRVPD + ISVPD + SNVPD + SLVPD + TRAND
FLOWD <<- SRFLD + BYFLD + DSFLD + GWFLD
}


## ------------------------------------------------------------------------
paccum<-function(){
  VRFLPI <<- VRFLPI + VRFLI * DTI
  SLFLPI <<- SLFLPI + SLFLI * DTI
  INFLPI <<- INFLPI + INFLI * DTI
  BYFLPI <<- BYFLPI + BYFLI * DTI
  DSFLPI <<- DSFLPI + DSFLI * DTI
  NTFLPI <<- NTFLPI + NTFLI * DTI
  TRANPI <<- TRANPI + TRANI * DTI

SRFLP <<- SRFLP + SRFL * DTI
SLFLP <<- SLFLP + SLFL * DTI
GWFLP <<- GWFLP + GWFL * DTI
SEEPP <<- SEEPP + SEEP * DTI

# sum flows for precip interval from components
sumii<-SUMI(NLAYER, BYFLPI, INFLPI, DSFLPI, TRANPI, DUMM, DUMM, BYFLP, INFLP, DSFLP, TRANP, dummy, dummy)
BYFLP<<-unlist(sumii[1])
INFLP<<-unlist(sumii[2])
DSFLP<<-unlist(sumii[3])
TRANP<<-unlist(sumii[4])
}


## ------------------------------------------------------------------------
yaccum<-function(){
ACCUMI(NLAYER, VRFLMI, INFLMI, BYFLMI, DSFLMI, NTFLMI, VRFLYI, INFLYI, BYFLYI, DSFLYI, NTFLYI)
ACCUMI(NLAYER, TRANMI, SLFLMI, DUMM, DUMM, DUMM, TRANYI, SLFLYI, DUMM, DUMM, DUMM)
ACCUM(SRFLM, SLFLM, GWFLM, SEEPM, dummy, SRFLY, SLFLY, GWFLY, SEEPY, dummy)
ACCUM(ISVPM, IRVPM, SNVPM, SLVPM, SFALM, ISVPY, IRVPY, SNVPY, SLVPY, SFALY)
ACCUM(RFALM, SINTM, RINTM, RSNOM, SMLTM, RFALY, SINTY, RINTY, RSNOY, SMLTY)
ACCUM(MESFLM, PTRANM, PINTM, dummy, dummy, MESFLY, PTRANY, PINTY, dummy, dummy)
SUMI(NLAYER, BYFLYI, INFLYI, DSFLYI, TRANYI, DUMM, DUMM, BYFLY, INFLY, DSFLY, TRANY, dummy, dummy)
}


## ------------------------------------------------------------------------
maccum<-function(){
accumi<-ACCUMI(NLAYER, VRFLDI, INFLDI, BYFLDI, DSFLDI, NTFLDI, VRFLMI, INFLMI, BYFLMI, DSFLMI, NTFLMI)
  VRFLMI<<-unlist(accumi[1])
  INFLMI<<-unlist(accumi[2])
  BYFLMI<<-unlist(accumi[3])
  DSFLMI<<-unlist(accumi[4])
  NTFLMI<<-unlist(accumi[5])

accumi<-ACCUMI(NLAYER, TRANDI, SLFLDI, DUMM, DUMM, DUMM, TRANMI, SLFLMI, DUMM, DUMM, DUMM)
TRANMI<<-unlist(accumi[1])
SLFLMI<<-unlist(accumi[2])

accumii<-ACCUM(SRFLD, SLFLD, GWFLD, SEEPD, dummy, SRFLM, SLFLM, GWFLM, SEEPM, dummy)
SRFLM<<-unlist(accumii[1])
SLFLM<<-unlist(accumii[2])
GWFLM<<-unlist(accumii[3])
SEEPM<<-unlist(accumii[4])

accumii<-ACCUM(ISVPD, IRVPD, SNVPD, SLVPD, SFALD, ISVPM, IRVPM, SNVPM, SLVPM, SFALM)
ISVPM<<-unlist(accumii[1])
IRVPM<<-unlist(accumii[2])
SNVPM<<-unlist(accumii[3] )
SLVPM<<-unlist(accumii[4])
SFALM<<-unlist(accumii[5])

acumii<-ACCUM(RFALD, SINTD, RINTD, RSNOD, SMLTD, RFALM, SINTM, RINTM, RSNOM, SMLTM)
RFALM<<-unlist(acumii[1])
SINTM<<-unlist(acumii[2])
RINTM<<-unlist(acumii[3])
RSNOM<<-unlist(acumii[4])
SMLTM<<-unlist(acumii[5])

acumii<-ACCUM(MESFLD, PTRAND, PINTD, dummy, dummy, MESFLM, PTRANM, PINTM, dummy, dummy)
MESFLM<<-unlist(acumii[1] )
PTRANM<<-unlist(acumii[2])
PINTM<<-unlist(acumii[3])

# sum flows for month from components
summi<-SUMI(NLAYER, BYFLMI, INFLMI, DSFLMI, TRANMI, DUMM, DUMM, BYFLM, INFLM, DSFLM, TRANM, dummy, dummy)
BYFLM<<-unlist(summi[1])
INFLM<<-unlist(summi[2])
DSFLM<<-unlist(summi[3])
TRANM<<-unlist(summi[4])
}


## ------------------------------------------------------------------------
daccum<-function(){
# accumulate above ground flows over day
ISVPD <<- ISVPD + ISVP * DTP
IRVPD <<- IRVPD + IRVP * DTP
SNVPD <<- SNVPD + SNVP * DTP
SLVPD <<- SLVPD + SLVP * DTP
SFALD <<- SFALD + SFAL * DTP
RFALD <<- RFALD + RFAL * DTP
SINTD <<- SINTD + SINT * DTP
RINTD <<- RINTD + RINT * DTP
RSNOD <<- RSNOD + RSNO * DTP
SMLTD <<- SMLTD + SMLT * DTP
MESFLD <<- MESFLD + MESFLP * DTP
PTRAND <<- PTRAND + PTRAN * DTP
PINTD <<- PINTD + PINT * DTP

# accumulate below ground flows over day
accumi<-ACCUMI(NLAYER, VRFLPI, INFLPI, BYFLPI, DSFLPI, NTFLPI, VRFLDI, INFLDI, BYFLDI, DSFLDI, NTFLDI)
VRFLDI<<-unlist(accumi[1]) 
INFLDI<<-unlist(accumi[2])
BYFLDI<<-unlist(accumi[3]) 
DSFLDI<<-unlist(accumi[4]) 
NTFLDI<<-unlist(accumi[5])

accumi<-ACCUMI(NLAYER, TRANPI, SLFLPI, DUMM, DUMM, DUMM, TRANDI, SLFLDI, DUMM, DUMM, DUMM)
TRANDI<<-unlist(accumi[1])
SLFLDI<<-unlist(accumi[2])

accum<-ACCUM(SRFLP, SLFLP, GWFLP, SEEPP, dummy, SRFLD, SLFLD, GWFLD, SEEPD, dummy)
SRFLD<<-unlist(accum[1])
SLFLD<<-unlist(accum[2])
GWFLD<<-unlist(accum[3])
SEEPD<<-unlist(accum[4])

# sum flows for day from components
summii<-SUMI(NLAYER, BYFLDI, INFLDI, DSFLDI, TRANDI, DUMM, DUMM, BYFLD, INFLD, DSFLD, TRAND, dummy, dummy)
BYFLD<<-unlist(summii[1])
INFLD<<-unlist(summii[2])
DSFLD<<-unlist(summii[3])
TRAND<<-unlist(summii[4])
}


## ------------------------------------------------------------------------
zpint<-function(){
VRFLPI<<-rep(0,ML)
INFLPI<<-rep(0,ML)
BYFLPI<<-rep(0,ML)
DSFLPI<<-rep(0,ML)
 NTFLPI<<-rep(0,ML)
 TRANPI<<-rep(0,ML)
 SLFLPI<<-rep(0,ML)
SRFLP<<-0
SLFLP<<-0
GWFLP<<-0
SEEPP<<-0
}


## ------------------------------------------------------------------------
zyear<-function(){
  VRFLYI<<-rep(0,ML)
  INFLYI<<-rep(0,ML)
  BYFLYI<<-rep(0,ML)
  DSFLYI<<-rep(0,ML)
 NTFLYI<<-rep(0,ML)
 TRANYI<<-rep(0,ML)
 SLFLYI<<-rep(0,ML) 
SRFLY<<-0
GWFLY<<-0
SEEPY<<-0
SLFLY<<-0
IRVPY<<-0
ISVPY<<-0
SLVPY<<-0
SNVPY<<-0
SFALY<<-0
RFALY<<-0
SINTY<<-0
RINTY<<-0
RSNOY<<-0
SMLTY<<-0
MESFLY<<-0
PTRANY<<-0
PINTY<<-0
}


## ------------------------------------------------------------------------
zmonth<-function(){
  VRFLMI<<-rep(0,ML) 
  INFLMI<<-rep(0,ML) 
  BYFLMI<<-rep(0,ML) 
  DSFLMI<<-rep(0,ML) 
  NTFLMI<<-rep(0,ML) 
  TRANMI<<-rep(0,ML) 
  SLFLMI<<-rep(0,ML)

  SRFLM<<-0
  GWFLM<<-0
  SEEPM<<-0
  SLFLM<<-0
  IRVPM<<-0 
  ISVPM<<-0 
  SLVPM<<-0 
  SNVPM<<-0
  SFALM<<-0
  RFALM<<-0 
  SINTM<<-0 
  RINTM<<-0 
  RSNOM<<-0 
  SMLTM<<-0
  MESFLM<<-0 
  PTRANM<<-0 
  PINTM<<-0
}


## ------------------------------------------------------------------------
zday<-function(){
VRFLDI<<-rep(0,ML)
INFLDI<<-rep(0,ML)
BYFLDI<<-rep(0,ML)
DSFLDI<<-rep(0,ML)
NTFLDI<<-rep(0,ML)
TRANDI<<-rep(0,ML) 
SLFLDI<<-rep(0,ML)
SRFLD<<-0
GWFLD<<-0
SEEPD<<-0
SLFLD<<-0

IRVPD<<-0
ISVPD<<-0
SLVPD<<-0
SNVPD<<-0
SFALD<<-0
RFALD<<-0
SINTD<<-0
RINTD<<-0
RSNOD<<-0
SMLTD<<-0
MESFLD<<-0
PTRAND<<-0
PINTD<<-0
}


## ------------------------------------------------------------------------
fnleap<-function(){
if ((YEARN %% 4 == 0)  && ((YEARN %% 100 != 0) || (YEARN %% 400 == 0))) {
return(TRUE)
}else{
return(FALSE)
}
}


## ------------------------------------------------------------------------
swchek<-function(i){
  if (SWATI[i] <= 0) {
      if(swatproblem >0){}
  }else if (SWATI[i] > SWATMX[i]){
  if (SWATI[i] > SWATMX[i] + 0.00001) {
    if(swatproblem >0){}
  }else{
# rounding error only
    SWATI[i] <<- SWATMX[i]
  }
}
}


## ------------------------------------------------------------------------
DOYF<-function(day,month, daymo){
  doyy<-0
  if(fnleap()){
    daymo[2]<-29
  }else{
    daymo[2]<-28
  }
  
  if(month>1)
    doyy<-daymo[1]+doyy
  if(month>2)
    doyy<-daymo[2]+doyy
  if(month>3)
    doyy<-daymo[3]+doyy
  if(month>4)
    doyy<-daymo[4]+doyy
  if(month>5)
    doyy<-daymo[5]+doyy
  if(month>6)
    doyy<-daymo[6]+doyy
  if(month>7)
    doyy<-daymo[7]+doyy
  if(month>8)
    doyy<-daymo[8]+doyy
  if(month>9)
    doyy<-daymo[9]+doyy
  if(month>10)
    doyy<-daymo[10]+doyy 
  if(month>11)
    doyy<-daymo[11]+doyy
  if(month>12)
    doyy<-daymo[12]+doyy 
  
  doyy<-doyy+day
  return(doyy)
}



## ----chunk4--------------------------------------------------------------
meteoFile <-read.table(file.path(projectpath,"Input_data","data_input_WB_1999_point.txt"))
MData<-matrix(meteoFile)

precFile <-read.table(file.path(projectpath,"Input_data","P_WBG_0060_Minuten_1999_point.txt"))
MhhData<-matrix(precFile)


## ----chunk5, warning=FALSE-----------------------------------------------
canopy_file <- read.table(file.path(projectpath,"Input_data","canopy.txt"), header = FALSE, fill=TRUE , dec = "." , sep = ",", strip.white = TRUE)

#canopy parameters
canopy_varname <- as.character(canopy_file[,1])
canopy_varvalue <- as.character(canopy_file [,2])
quant_canopy <- 18

for (i in 1:quant_canopy){
  a <- unlist(canopy_varname[i])
  b <- as.numeric(unlist(canopy_varvalue[i]))
  assign(a,b)
  rm(a,b)
}

#ROOTDEN
roota <- strsplit(canopy_varname[20],"[ ]+")
rootb <- unlist(roota)
rootc <- as.numeric(rootb)
rootaa <- strsplit(canopy_varname[22],"[ ]+")
rootbb <- unlist(rootaa)
rootcc <- as.numeric(rootbb)
root_a <- c(rootc,rootcc)
root_a[is.na(root_a)] <- 0

rootd <- strsplit(canopy_varname[21],"[ ]+")
roote <- unlist(rootd)
rootf <- as.numeric(roote)
rootdd <- strsplit(canopy_varname[23],"[ ]+")
rootee <- unlist(rootdd)
rootff <- as.numeric(rootee)
root_b <- c(rootf,rootff)
root_b[is.na(root_b)] <- 0

rootden <- as.vector(rbind(root_a,root_b)) 
assign(canopy_varname[19],rootden)

rm(roota,rootb,rootc,rootd,roote,rootf,rootaa,rootbb,rootcc,root_a,rootdd,rootee,rootff,root_b,quant_canopy,canopy_varvalue,canopy_varname,rootden)


## ----chunk6--------------------------------------------------------------
fixed_file <-read.table(file.path(projectpath,"Input_data","fixed.txt"), header = FALSE, fill=TRUE, sep = "," , dec = "." , row.names = 1)

fixed_varname <- rownames(fixed_file)
fixed_varvalue <- fixed_file [,c(1)]

quant_fixed = length(fixed_varname) - 1

for (i in 1:quant_fixed){
  a <- unlist(fixed_varname[i])
  b <- unlist(fixed_varvalue[i])
  assign(a,b)
  rm(a,b)
}

rm(fixed_varname,fixed_varvalue,quant_fixed)


## ----chunk7--------------------------------------------------------------
flow_file <-read.table(file.path(projectpath,"Input_data","flow_standard.txt"), header = FALSE, fill=TRUE, sep = "," , dec = "." , row.names = 1)

flow_varname <- rownames(flow_file)
flow_varvalue <- flow_file [,c(1)]

quant_flow = length(flow_varname) - 1

for (i in 1:quant_flow){
  a <- unlist(flow_varname[i])
  b <- unlist(flow_varvalue[i])
  assign(a,b)
  rm(a,b)
}

rm(flow_varname,flow_varvalue,quant_flow)


## ----chunk8--------------------------------------------------------------
initial_file <-read.table(file.path(projectpath,"Input_data","initial.txt"),header = FALSE,fill=TRUE,sep = ",",dec = ".",row.names = 1)

initial_varname <- rownames(initial_file)
initial_varvalue <- as.character(initial_file [,1])
quant_initial <- 4

for (i in 1:quant_initial){
  a <- unlist(initial_varname[i])
  b <- as.numeric(unlist(initial_varvalue[i]))
  assign(a,b)
  rm(a,b)
}

#PSIMIN
PSIMIN <- as.numeric(initial_varvalue[6:30])

rm(initial_varname,initial_varvalue,quant_initial)


## ----chunk9--------------------------------------------------------------
location_file <-read.table(file.path(projectpath,"Input_data","location.txt"), header = FALSE, fill=TRUE, sep = "," , dec = ".", strip.white = TRUE)

location_varname <- as.character(location_file[,1])
location_varvalue <- as.character(location_file [,2])
quant_location <- 5

for (i in 1:quant_location){
  a <- unlist(location_varname[i])
  b <- as.numeric(unlist(location_varvalue[i]))
  assign(a,b)
  rm(a,b)
}

#DURATN
duratna <- strsplit(location_varvalue[6],"[ ]+")
duratnb <- unlist(duratna)
duratn<- as.numeric(duratnb)
assign(location_varname[6],duratn)

rm(duratna,duratnb,duratn)

#RELHT
hta <- strsplit(location_varvalue[7],"[ ]+")
htb<- unlist(hta)
ht_a <- suppressWarnings(as.numeric(htb))
ht_a[is.na(ht_a)] <- 0

htc <- strsplit(location_varname[8],"[ ]+")
htd <- unlist(htc)
ht_b <- suppressWarnings(as.numeric(htd))
ht_b[is.na(ht_b)] <- 0

relht <- as.vector(rbind(ht_a,ht_b)) 
assign(location_varname[7],relht)

rm(hta,htb,htc,htd,ht_a,ht_b,relht)

#RELLAI
laia <- strsplit(location_varvalue[9],"[ ]+")
laib<- unlist(laia)
lai_a <- suppressWarnings(as.numeric(laib))
lai_a[is.na(lai_a)] <- 0

laic <- strsplit(location_varname[10],"[ ]+")
laid <- unlist(laic)
lai_b <- suppressWarnings(as.numeric(laid))
lai_b[is.na(lai_b)] <- 0

rellai <- as.vector(rbind(lai_a,lai_b)) 
assign(location_varname[9],rellai)

rm(laia,laib,laic,laid,lai_a,lai_b,rellai)
rm(location_varname,location_varvalue,quant_location)


## ----chunk10-------------------------------------------------------------
soil_file <-read.table(file.path(projectpath,"Input_data","soil.txt"), header = FALSE, fill=TRUE , dec = "." , sep = "", strip.white = TRUE)

soil_var1 <- as.character(soil_file[,1])
soil_var2 <- as.character(soil_file[,2])
soil_var3 <- as.character(soil_file[,3])
soil_var4 <- as.character(soil_file[,4])
soil_var5 <- as.character(soil_file[,5])
soil_var6 <- as.character(soil_file[,6])
soil_var7 <- as.character(soil_file[,7])
soil_var8 <- as.character(soil_file[,8])

# NLAYER
a_nlayer <- unlist(soil_var1[1])
a_nlayer <- as.character(gsub(",","",a_nlayer))
b_nlayer <- as.numeric(unlist(soil_var2[1]))
assign(a_nlayer,b_nlayer)
rm(a_nlayer,b_nlayer)

# THICK
a_thick <- unlist(soil_var1[2])
b_thick <- as.numeric(unlist(soil_var1[3:27]))
assign(a_thick,b_thick)
rm(a_thick,b_thick)

# STONEF
a_stonef <- unlist(soil_var2[2])
b_stonef <- as.numeric(unlist(soil_var2[3:27]))
assign(a_stonef,b_stonef)
rm(a_stonef,b_stonef)

# PSIF
a_psif <- unlist(soil_var3[2])
b_psif <- as.numeric(unlist(soil_var3[3:27]))
assign(a_psif,b_psif)
rm(a_psif,b_psif)

# THETAF
a_thetaf <- unlist(soil_var4[2])
b_thetaf <- as.numeric(unlist(soil_var4[3:27]))
assign(a_thetaf,b_thetaf)
rm(a_thetaf,b_thetaf)

# THSAT
a_thsat <- unlist(soil_var5[2])
b_thsat <- as.numeric(unlist(soil_var5[3:27]))
assign(a_thsat,b_thsat)
rm(a_thsat,b_thsat)

# BEXP
a_bexp <- unlist(soil_var6[2])
b_bexp <- as.numeric(unlist(soil_var6[3:27]))
assign(a_bexp,b_bexp)
rm(a_bexp,b_bexp)

# KF
a_kf <- unlist(soil_var7[2])
b_kf <- as.numeric(unlist(soil_var7[3:27]))
assign(a_kf,b_kf)
rm(a_kf,b_kf)

# WETINF
a_wetinf <- unlist(soil_var8[2])
b_wetinf <- as.numeric(unlist(soil_var8[3:27]))
assign(a_wetinf,b_wetinf)
rm(a_wetinf,b_wetinf)

rm(soil_var1,soil_var2,soil_var3,soil_var4,soil_var5,soil_var6,soil_var7,soil_var8)



## ----chunkini------------------------------------------------------------
if((runflag == 0) || (runflag == 1)){
  DAYMO[1] = 31
  DAYMO[2] = 28
  DAYMO[3] = 31
  DAYMO[4] = 30 
  DAYMO[5] = 31
  DAYMO[6] = 30
  DAYMO[7] = 31
  DAYMO[8] = 31
  DAYMO[9] = 30
  DAYMO[10] = 31
  DAYMO[11] = 30
  DAYMO[12] = 31
	IDAY =1
	IInterValDay=1
	NDAYS=length(MData[[1]])
	NITSR = 0
	NITSY = 0
	NITSM = 0
	YEARN = as.numeric(MData[[1]][IDAY])
	
	daymax=NDAYS-IDAY+1
	maxF=0
	timeseries_prec=rep(0,daymax)
	timeseries_evp=rep(0,daymax)
	timeseries_flow=rep(0,daymax)
	timeseries_rnet=rep(0,daymax)
	timeseries_ptran=rep(0,daymax)
	timeseries_irvp=rep(0,daymax)
	timeseries_isvp=rep(0,daymax)
	timeseries_snow=rep(0,daymax)
	timeseries_swat=rep(0,daymax)
	timeseries_pint=rep(0,daymax)
	timeseries_snvp=rep(0,daymax)
	timeseries_slvp=rep(0,daymax)
	timeseries_trand=rep(0,daymax)
	timeseries_mesfld=rep(0,daymax)
	timeseries_smltd=rep(0,daymax)
	timeseries_slfld=rep(0,daymax)
	timeseries_rfald=rep(0,daymax)
	timeseries_sfald=rep(0,daymax)
	timeseries_awat=rep(0,daymax)
	timeseries_adef=rep(0,daymax)
	timeseries_sintd=rep(0,daymax)
	timeseries_rintd=rep(0,daymax)
	timeseries_rthrd=rep(0,daymax)
	timeseries_sthrd=rep(0,daymax)
	timeseries_rsnod=rep(0,daymax)
	
	if( YEARN < 100){
		if(YEARN > 20){
			YEARN = YEARN + 1900
		}else{
			YEARN = YEARN + 2000
		}
	}
  MONTHN = as.numeric(MData[[2]][IDAY])
	DOM = as.numeric(MData[[3]][IDAY])
	DOY = DOY=DOYF(DOM,MONTHN,DAYMO)
		
	  if (fnleap()) {
			DAYMO[2] = 29
		}else{
			DAYMO[2] = 28
		}
	if (SUBDAYDATA) {
  	DTP = DT / NPINT
	}else{
		DTP = DT
	}
# zero accumulators
	zyear()
	zmonth()
# initial values
	SNOW = SNOWIN
	GWAT = GWATIN
	INTR = INTRIN
	INTS = INTSIN
	for( i in 1:NLAYER){
		PSIM[i] = PSIMIN[i]
	}
# soil water parameters and initial variables
	soilp<-SOILPAR()
	  PSIG<-unlist(soilp[2])
	  SWATMX<-unlist(soilp[3])
	  WETF<-unlist(soilp[4])
	  WETC<-unlist(soilp[5])
	  CHM<-unlist(soilp[6])
	  CHN<-unlist(soilp[7]) 
	  WETNES<-unlist(soilp[8])
	  SWATI<-unlist(soilp[9])
	  KSAT<-unlist(soilp[10])
# ^^
# initial soil water variables
	soil<-SOILVAR()
	  PSITI<-soil[1:ML]
	  THETA<-soil[(ML+1):(2*ML)]
	  KK<-soil[(2*ML+1):(3*ML)]
	  SWAT<-soil[(3*ML+1)]
# ^^
# initial total water in system
	STORD = INTR + INTS + SNOW + SWAT + GWAT
	STORM = STORD
	STORY = STORD
# any initial snow has zero liquid water and cold content
	CC = 0
	SNOWLQ = 0
}


## ----chunkpara-----------------------------------------------------------
# parameter conversions
GLMAX = GLMAXC / 100
GLMIN = GLMINC / 100
LAT = LATD / 57.296
ESLOPE = ESLOPED / 57.296
DSLOPE = DSLOPED / 57.296
ASPECT = ASPECTD / 57.296
# equivalent slope for radiation calculations
equi<-EQUIVSLP(LAT, ESLOPE, ASPECT)
  L1<-unlist(equi[1])
  L2<-unlist(equi[2])
# ^^
# infiltration parameters
infpa<-INFPAR(INFEXP, IDEPTH, NLAYER, THICK)
  ILAYER<-unlist(infpa[1])
  INFRAC<-unlist(infpa[2])
# ^^
# source area parameters
srfp<-SRFPAR(QDEPTH, NLAYER, THETAF, THICK, STONEF, SWATMX)
  QLAYER<-unlist(srfp[1]) 
  SWATQX<-unlist(srfp[2])
  SWATQF<-unlist(srfp[3])
# ^^
# root density parameters
RELDEN<-RTDEN(ROOTDEN, NLAYER, THICK)


## ----chunkmodel----------------------------------------------------------
while( IDAY <= NDAYS){  
	NITSD = 0
		subdatafileline(IDAY)
  if( IDAY == INIDAYS + 1){
# end of initialization, reinitialize year and month accumulators
		STORD = INTR + INTS + SNOW + SWAT + GWAT
		STORM = STORD
		STORY = STORD
		NITSY = 0
		NITSM = 0
		zyear()
		zmonth()
  }
# calculate derived variables
	MSBSETVARS()
#
#* * * * *  B E G I N   D A Y - N I G H T   E T   L O O P  * * * * * * * * *
#potential and actual interception, evaporation, and transpiration
	MSBDAYNIGHT()
#
#* * * * * * * *  E N D   D A Y - N I G H T   L O O P  * * * * * * * * * *
# average rates over day
	PTRAN = (PTR[1] * DAYLEN + PTR[2] * (1 - DAYLEN)) / DT
	GEVP = (GER[1] * DAYLEN + GER[2] * (1 - DAYLEN)) / DT
	PINT = (PIR[1] * DAYLEN + PIR[2] * (1 - DAYLEN)) / DT
	GIVP = (GIR[1] * DAYLEN + GIR[2] * (1 - DAYLEN)) / DT
	for(i in 1:NLAYER){
		TRANI[i] = (ATRI[1, i] * DAYLEN + ATRI[2, i] * (1 - DAYLEN)) / DT
	}
# zero daily integrators
	zday()
#
#* * * * * * * * B E G I N   P R E C I P   I N T E R V A L * * * * * * * * *
for( N in 1:NPINT){  
	if (SUBDAYDATA){
    subprfileline(IInterValDay)
	  if (MESFLP <= -0.01) {MESFLP = MESFL / DT}
	}else{
# precip data from data file
		PREINT = PRECIN / DT
		MESFLP = MESFL / DT
	}
# interception and snow accumulation/melt
	MSBPREINT()
# initialize for iterations
# initial time remaining in iteration time step = precip time step
	DTRI = DTP
# initialize iteration counter
	NITS = 0
# zero precip interval integrators
	zpint()
#
#  *  *  *  *  *  *  B E G I N   I T E R A T I O N   *  *  *  *  *  *  *  *
while(!(DTRI <= 0)){  
		NITS = NITS + 1
# check for events
		if (NITS %% 100 == 0) {}
# water movement through soil
		MSBITERATE() 
# iteration calculations
# calculate SLFLI vertical macropore infiltration out of layer
		SLFLI[1] = SLFL - INFLI[1] - BYFLI[1]
		if (ILAYER >= 2){
		  if (NLAYER >= ILAYER +1){
		    for (i in 2:ILAYER){ 
		    # does not execute if ILAYER% = 1 or 0
			    SLFLI[i] = SLFLI[i - 1] - INFLI[i] - BYFLI[i]
		    }
		    for( i in (ILAYER + 1):NLAYER){ 
		    # does not execute if NLAYER% < ILAYER% + 1
			    SLFLI[i] = 0
		    }
		  }
		}
# integrate below ground storages over iteration interval
		for( i in 1:NLAYER){
			SWATI[i] = SWATI[i] + NTFLI[i] * DTI
		}
		GWAT = GWAT + (VRFLI[NLAYER] - GWFL - SEEP) * DTI
# new soil water variables and test for errors
		for (i in 1:NLAYER){
			swchek(i)
			WETNES[i] = SWATI[i] / SWATMX[i]
			PSIM[i] = FPSIMF(WETNES[i], PSIF[i], BEXP[i], WETINF[i], WETF[i], CHM[i], CHN[i])
		}
		soil<-SOILVAR()
		   PSITI<-soil[1:ML]
		   THETA<-soil[(ML+1):(2*ML)]
		   KK<-soil[(2*ML+1):(3*ML)]
		   SWAT<-soil[(3*ML+1)]
# ^^
# iteration output
# flows accumulated over precip interval
	paccum()
# time remaining in precipitation time-step
	DTRI = DTRI - DTI
	NITSR = NITSR + 1  # for visible display of iterations
}
#
#  *  *  *  *   E N D   i T E R A T I O N    L O O P  *  *  *  *  *  *  *  *
# display iterations
# integrate interception storages over precip interval
INTS = INTS + (SINT - ISVP) * DTP
INTR = INTR + (RINT - IRVP) * DTP
#  flows for precip interval summed from components
psum()
# precipitation interval output
# flows accumulated over day
daccum()
# accumulate iterations
  NITSD = NITSD + NITS
  NITSM = NITSM + NITS
  NITSY = NITSY + NITS
  IInterValDay<-IInterValDay+1
}
#
#* * * * *  E N D   P R E C I P   I N T E R V A L   L O O P  * * * * * * * *
# flows for day summed from components
dsum()
# check for water balance error
BALERD = STORD - (INTR + INTS + SNOW + SWAT + GWAT) + PRECD - EVAPD - FLOWD - SEEPD
STORD = INTR + INTS + SNOW + SWAT + GWAT
# flows accumulated over month
maccum()
# date checking on
if(DOM == DAYMO[MONTHN]){
# set up for next month
zmonth()
MONTHN = MONTHN + 1
DOM = 0
NITSM = 0
}  # for end of month
if (MONTHN == 13) {
# end of year
# set up for next year
  MONTHN = 1
  DOM = 1
  DOY = 1
  YEARN = YEARN + 1
  zyear()
  if (fnleap() ){
    DAYMO[2] = 29
  }else{
    DAYMO[2] = 28
  }
NITSY = 0
NITSM = 0
} 
#set up for next day
IDAY = IDAY + 1
  MONTHN = as.numeric(MData[[2]][IDAY])
  DOM = as.numeric(MData[[3]][IDAY])
  YEARN = as.numeric(MData[[1]][IDAY])
  if(IDAY <= NDAYS)
  DOY=DOYF(DOM,MONTHN,DAYMO)
 
#* * * I N P U T   W E A T H E R   L I N E   F R O M   D F I L E * * *
#subdatafileline()
#
# ***************   E N D    D A Y   L O O P    **************************
	timeseries_prec[daymax-NDAYS+IDAY-1]<-PRECD
	timeseries_evp[daymax-NDAYS+IDAY-1]<-EVAPD
	timeseries_flow[daymax-NDAYS+IDAY-1]<-FLOWD
	timeseries_rnet[daymax-NDAYS+IDAY-1]<-RNET
	timeseries_irvp[daymax-NDAYS+IDAY-1]<-IRVPD
	timeseries_isvp[daymax-NDAYS+IDAY-1]<-ISVPD
	timeseries_ptran[daymax-NDAYS+IDAY-1]<-PTRAND
	timeseries_snow[daymax-NDAYS+IDAY-1]<-SNOW
	timeseries_swat[daymax-NDAYS+IDAY-1]<-SWAT
	timeseries_pint[daymax-NDAYS+IDAY-1]<-PINTD
	timeseries_snvp[daymax-NDAYS+IDAY-1]<-SNVPD 
	timeseries_slvp[daymax-NDAYS+IDAY-1]<-SLVPD 
	timeseries_trand[daymax-NDAYS+IDAY-1]<-TRAND
	timeseries_mesfld[daymax-NDAYS+IDAY-1]<-MESFLD
	timeseries_smltd[daymax-NDAYS+IDAY-1]<-SMLTD
	timeseries_slfld[daymax-NDAYS+IDAY-1]<-SLFLD
	timeseries_rfald[daymax-NDAYS+IDAY-1]<-RFALD
	timeseries_awat[daymax-NDAYS+IDAY-1]<-AWAT
	timeseries_adef[daymax-NDAYS+IDAY-1]<-ADEF
	timeseries_sintd[daymax-NDAYS+IDAY-1]<-SINTD
	timeseries_rintd[daymax-NDAYS+IDAY-1]<-RINTD
	timeseries_sfald[daymax-NDAYS+IDAY-1]<-SFALD
	timeseries_rthrd[daymax-NDAYS+IDAY-1]<-RTHRD
	timeseries_sthrd[daymax-NDAYS+IDAY-1]<-STHRD
	timeseries_rsnod[daymax-NDAYS+IDAY-1]<-RSNOD
}
