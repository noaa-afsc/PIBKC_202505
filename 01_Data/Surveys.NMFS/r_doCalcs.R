#--do calculations for NMFS survey data
#--NOTE: output units for biomass are 1000's t, for abundance are millions

rstudio = rstudioapi::isAvailable();
if (!rstudio) 
  stop("Must run/source r_doCalcs.R from the RStudio console.");

dirPrj = rstudioapi::getActiveProject();
setup = wtsUtilities::getObj(file.path(dirPrj,"rda_AssessmentSetup.RData"));

#--if file is being run from RStudio, set dirOut to path to file for output
if (rstudio){
  dirInp = dirname(rstudioapi::getActiveDocumentContext()$path);
  dirOut = dirname(rstudioapi::getActiveDocumentContext()$path);
}
#--otherwise dirOut should have been already set
if (!exists("dirInp")) stop("RStudio not running and dirInp does not exist");
if (!exists("dirOut")) stop("RStudio not running and dirOut does not exist");

#--set options
species    = 'BKC';
strataType = '2015';
minYr = 1975;
maxYr = setup$asmtYr;
downloadDate = "August 15, 2023";#--TODO: move to assessment setup(?)
dirData = "~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current";
fn.SD   = "PIBKC_SurveyStrata.csv"
fn.CH   = "PIBKC_HaulData.csv"

#--load required packages
require(tcsamSurveyData);

#--required preliminary calculations
options(stringsAsFactors=FALSE);
codes.TS<-Codes.TrawlSurvey();
Sum<-wtsUtilities::Sum;#sum function with na.rm=TRUE

outDir<-dirOut;
if (!dir.exists(outDir)) dir.create(outDir);


HaulTypes = NULL;#get ALL haul types for new time series
minSize   = 0;#minimum size of individuals to extract
cutpts    = seq(from=minSize,to=200,by=5);

verbosity = 0;##0=off,1=minimal,2=full
  
#--define size groups for males
#----immature/mature groups for females are based on abdomen morphology and egg condition
dfrZGT<-rbind(
           data.frame(sex=  'male',`size range`="< 120 mm CL",category="immature male"),
           data.frame(sex=  'male',`size range`="> 119 mm CL",category="mature male"),
           data.frame(sex=  'male',`size range`="< 135 mm CL",category="sublegal male"),
           data.frame(sex=  'male',`size range`="> 134 mm CL",category="legal male"),
           data.frame(sex='female',`size range`="all",        category="all females"),
           data.frame(sex=  'male',`size range`="all",        category="all males"));

#--select strata definitions
tmp0<-readr::read_csv(file.path(dirData,fn.SD),guess_max=1000000)[,1:8];
names(tmp0)[8]<-"TOTAL_AREA";
dfr.SD<-selectStrata.TrawlSurvey(tmp0,
                                 species=species,
                                 strataType=strataType,
                                 export=FALSE,
                                 verbosity=verbosity);

#--load crabhaul files
dfr.crabhauls<-NULL;
for (f in fn.CH){
  tmp1<-readr::read_csv(file.path(dirData,f),skip=5,guess_max=1000000);
  dfr.crabhauls<-rbind(dfr.crabhauls,tmp1);
}

rm(tmp0,tmp1);


#--select haul data
dfr.HD<-selectHauls.TrawlSurvey(dfr.SD,
                                tbl=dfr.crabhauls,
                                YearRange=c(minYr,maxYr),
                                export=FALSE,
                                verbosity=verbosity);

#--select individuals
dfrID.ImmF<-selectIndivs.TrawlSurvey(dfr.HD,
                                     tbl=dfr.crabhauls,
                                     col.Size="LENGTH",
                                     sex='FEMALE',
                                     shell_condition='ALL',
                                     maturity='IMMATURE',
                                     calcMaleMaturity=FALSE,
                                     minSize=0,
                                     maxSize=Inf,
                                     export=FALSE,
                                     verbosity=verbosity);
dfrID.ImmF$SEX<-"1. immature females";
dfrID.ImmF$SEXp<-"females";

dfrID.ImmM<-selectIndivs.TrawlSurvey(dfr.HD,
                                     tbl=dfr.crabhauls,
                                     col.Size="LENGTH",
                                     sex='MALE',
                                     shell_condition='ALL',
                                     maturity='ALL',
                                     calcMaleMaturity=FALSE,
                                     minSize=0,
                                     maxSize=setup$maleMatZ-1,
                                     export=FALSE,
                                     verbosity=verbosity);
dfrID.ImmM$SEX<-"2. immature males";
dfrID.ImmM$SEXp<-"males";

dfrID.MatF<-selectIndivs.TrawlSurvey(dfr.HD,
                                     tbl=dfr.crabhauls,
                                     col.Size="LENGTH",
                                     sex='FEMALE',
                                     shell_condition='ALL',
                                     maturity='MATURE',
                                     calcMaleMaturity=FALSE,
                                     minSize=0,
                                     maxSize=Inf,
                                     export=FALSE,
                                     verbosity=verbosity);
dfrID.MatF$SEX<-"3. mature females";
dfrID.MatF$SEXp<-"females";

dfrID.MatM<-selectIndivs.TrawlSurvey(dfr.HD,
                                     tbl=dfr.crabhauls,
                                     col.Size="LENGTH",
                                     sex='MALE',
                                     shell_condition='ALL',
                                     maturity='ALL',
                                     calcMaleMaturity=FALSE,
                                     minSize=setup$maleMatZ,
                                     maxSize=Inf,
                                     export=FALSE,
                                     verbosity=verbosity);
dfrID.MatM$SEX<-"4. mature males";
dfrID.MatM$SEXp<-"males";

dfrID.SubL<-selectIndivs.TrawlSurvey(dfr.HD,
                                     tbl=dfr.crabhauls,
                                     col.Size="LENGTH",
                                     sex='MALE',
                                     shell_condition='ALL',
                                     maturity='ALL',
                                     calcMaleMaturity=FALSE,
                                     minSize=0,
                                     maxSize=setup$maleLglZ-1,
                                     export=FALSE,
                                     verbosity=verbosity);
dfrID.SubL$SEX<-"5. sublegal males";
dfrID.SubL$SEXp<-"males";

dfrID.LglM<-selectIndivs.TrawlSurvey(dfr.HD,
                                     tbl=dfr.crabhauls,
                                     col.Size="LENGTH",
                                     sex='MALE',
                                     shell_condition='ALL',
                                     maturity='ALL',
                                     calcMaleMaturity=FALSE,
                                     minSize=setup$maleLglZ,
                                     maxSize=Inf,
                                     export=FALSE,
                                     verbosity=verbosity);
dfrID.LglM$SEX<-"6. legal males";
dfrID.LglM$SEXp<-"males";

dfrID.AllF<-selectIndivs.TrawlSurvey(dfr.HD,
                                     tbl=dfr.crabhauls,
                                     col.Size="LENGTH",
                                     sex='FEMALE',
                                     shell_condition='ALL',
                                     maturity='ALL',
                                     calcMaleMaturity=FALSE,
                                     minSize=0,
                                     maxSize=Inf,
                                     export=FALSE,
                                     verbosity=verbosity);
dfrID.AllF$SEX<-"7. all females";
dfrID.AllF$SEXp<-"females";

dfrID.AllM<-selectIndivs.TrawlSurvey(dfr.HD,
                                     tbl=dfr.crabhauls,
                                     col.Size="LENGTH",
                                     sex='MALE',
                                     shell_condition='ALL',
                                     maturity='ALL',
                                     calcMaleMaturity=FALSE,
                                     minSize=0,
                                     maxSize=Inf,
                                     export=FALSE,
                                     verbosity=verbosity);
dfrID.AllM$SEX<-"8. all males";
dfrID.AllM$SEXp<-"males";

dfr.ID<-rbind(dfrID.ImmF,dfrID.MatF,dfrID.AllF,
              dfrID.ImmM,dfrID.MatM,dfrID.SubL,dfrID.LglM,dfrID.AllM);

#--Calculate annual survey abundance and biomass by sex/size group for the Pribilof District
#----calculate CPUE by haul
dfrCPUE.ByX <-calcCPUE.ByHaul(dfr.HD,
                              dfr.ID,
                              bySex=TRUE,
                              byMaturity=FALSE,
                              byShellCondition=FALSE,
                              bySize=FALSE,
                              export=FALSE,
                              verbosity=verbosity);

#----calculate CPUE by station
dfrCPUE.ByX <-calcCPUE.ByStation(dfr.SD,
                                 dfrCPUE.ByX,
                                 export=FALSE,
                                 verbosity=verbosity);

#----calculate biomass/abundance by stratum
dfrACD.ByX.ByStrata <-calcAB.ByStratum(tbl_strata=dfr.SD,
                                       tbl_cpue=dfrCPUE.ByX,
                                       export=FALSE,
                                       verbosity=verbosity);

#----calculate biomass/abundance for the Pribilof District
tmp<-calcAB.EBS(dfrACD.ByX.ByStrata,
                export=FALSE,
                verbosity=verbosity);
tmp$STRATUM<-"Pribilof District";
tmp$SEXp<-"males";
idx<-tmp$SEX %in% c("1. immature females","3. mature females","7. all females");
tmp$SEXp[idx]<-"females";
write.csv(tmp,file.path(outDir,"dfrACD.BySexGroup.PribDistrict.csv"),row.names=FALSE);

lvls<-c("1. immature females","2. immature males",
        "3. mature females","4. mature males",
        "5. sublegal males","6. legal males",
        "7. all females","8. all males");
lbls<-c("immature females","immature males",
        "mature females","mature males",
        "sublegal males","legal males", 
        "all females","all males");
tmp$SEX<-factor(tmp$SEX,levels=lvls,labels=lbls);

dfrACD.ByX.PD<-tmp;
rm(tmp,dfrACD.ByX.ByStrata);

#--get stats on biomass/abundance
tmp = dfrACD.ByX.PD |> 
        dplyr::select(SEXp,SEX,YEAR,totABUNDANCE,totBIOMASS) |>
        dplyr::mutate(decade=10*floor(YEAR/10));
dfrStats = tmp |> dplyr::group_by(SEXp,SEX,decade) |>
                   dplyr::summarize(meanAbd=mean(totABUNDANCE),
                                    maxAbd=max(totABUNDANCE),
                                    meanBio=mean(totBIOMASS),
                                    maxBio=max(totBIOMASS)) |>
                  dplyr::ungroup() |> 
                  dplyr::arrange(SEXp,SEX,decade);

#--Size compositions
##replace SEX categories with SEXp categories
dfrID.AllX<-rbind(dfrID.AllF,dfrID.AllM);
dfrID.AllX$SEX<-dfrID.AllX$SEXp;
##calculate size comps by stratum for all individuals
dfrZCs.ByXS.ByS<-calcSizeComps.ByStratum(dfr.SD,
	                                  tbl_hauls=dfr.HD,
	                                  tbl_indivs=dfrID.AllX,
	                                  avgHaulsByStation=TRUE,
	                                  bySex=TRUE,
	                                  byShellCondition=TRUE,
	                                  byMaturity=FALSE,
	                                  cutpts=cutpts,
	                                  truncate.low=TRUE,
	                                  truncate.high=FALSE,
	                                  export=FALSE,
	    							                verbosity=verbosity);

dfrZCs.ByXS.PD<-calcSizeComps.EBS(dfrZCs.ByXS.ByS,
                                  export=FALSE,
                                  verbosity=verbosity);
dfrZCs.ByXS.PD$STRATUM<-"Pribilof District";
rm(dfrZCs.ByXS.ByS);

#--calculate survey CPUE
#----calculate CPUE by haul
dfrCPUE.ByX<-calcCPUE.ByHaul(dfr.HD,
                            dfrID.AllX,
                            bySex=TRUE,
                            byMaturity=FALSE,
                            byShellCondition=FALSE,
                            bySize=FALSE,
                            export=FALSE,
                            verbosity=verbosity);

#----calculate CPUE by station
dfrCPUE.ByX<-calcCPUE.ByStation(dfr.SD,
                                dfrCPUE.ByX,
                                export=FALSE,
                                verbosity=verbosity);

#--stations in latest survey
surveyGridLayers<-tcsamSurveyData::gisGetSurveyGridLayers();
dfr.SDp = dfr.SD |> subset(YEAR==maxYr);
sd.grid<-wtsGIS::mergeDataframeWithLayer(dfr.SDp,surveyGridLayers$grid,
                                          dataID="GIS_STATION",geomsID="STATION_ID",
                                          sfJoinType="left join",
                                          spAllData=FALSE,spDuplicateGeoms=TRUE);
sd.stns<-wtsGIS::mergeDataframeWithLayer(dfr.SDp,surveyGridLayers$stations,
                                          dataID="GIS_STATION",geomsID="ID",
                                          sfJoinType="left join",
                                          spAllData=FALSE,spDuplicateGeoms=TRUE);

#--environmental data
evs.csv  <-tcsamSurveyData::calcEnvData.ByStation(dfr.SD,dfr.HD,verbosity=verbosity);
evs.grid<-wtsGIS::mergeDataframeWithLayer(evs.csv,surveyGridLayers$grid,
                                          dataID="GIS_STATION",geomsID="STATION_ID",
                                          sfJoinType="left join",
                                          spAllData=FALSE,spDuplicateGeoms=TRUE);
evs.stns<-wtsGIS::mergeDataframeWithLayer(evs.csv,surveyGridLayers$stations,
                                          dataID="GIS_STATION",geomsID="ID",
                                          sfJoinType="left join",
                                          spAllData=FALSE,spDuplicateGeoms=TRUE);

rm(dfr.crabhauls,dfr.HD,dfr.ID,dfr.SD,
   dfrID.ImmF,dfrID.MatF,dfrID.AllF,dfrID.ImmM,dfrID.MatM,dfrID.SubL,dfrID.LglM,dfrID.AllM,
   idx,lbls,lvls,Sum,verbosity);

tsr=list(dirData=dirData,
         units_biomass="1000s t",
         units_abundance="millions",
         downloadDate=downloadDate,
         fn.SD=fn.SD,
         fn.CH=fn.CH,
         species=species,
         strataType=strataType,
         minYr=minYr,
         maxYr=maxYr,
         minSize=minSize,
         cutpts=cutpts,
         outDir=outDir,
         dfrZGT=dfrZGT,
         dfrACD.ByX.PD=dfrACD.ByX.PD,
         dfrStats=dfrStats,
         dfrCPUE.ByX=dfrCPUE.ByX,
         dfrZCs.ByXS.PD=dfrZCs.ByXS.PD,
         sd.stns=sd.stns,
         sd.grid=sd.grid,
         evs.stns=evs.stns,
         evs.grid=evs.grid);
wtsUtilities::saveObj(tsr,fn=setup$fnData_TSs);

