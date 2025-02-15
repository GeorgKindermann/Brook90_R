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
