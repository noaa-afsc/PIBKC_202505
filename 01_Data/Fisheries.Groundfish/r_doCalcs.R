#--calculate bycatch and discard mortality in the groundfish fisheries
#----NOTE: input bycatch biomass is in KG!
Sum = wtsUtilities::Sum;
options(stringsAsFactors=FALSE);
verbosity<-0;##0=off,1=minimal,2=full
gisDir = "~/Work/StockAssessments-Crab/Data/GIS/Shapefiles";

rstudio = rstudioapi::isAvailable();
if (!rstudio) 
  stop("Must run/source r_doCalcs.R from the RStudio console.");

dirPrj = rstudioapi::getActiveProject();
setup = wtsUtilities::getObj(file.path(dirPrj,"rda_AssessmentSetup.RData"));

#--if file is being run from RStudio, set dirOut to path to file for output
if (rstudio){
  dirOut = dirname(rstudioapi::getActiveDocumentContext()$path);
}
#--otherwise dirOut should have been already set
if (!exists("dirOut")) stop("RStudio not running and dirOut does not exist");

assmtYr = setup$asmtYr; #assessment year
strYr = paste0(assmtYr-1,"/",stringr::str_sub(paste0(assmtYr),3,4));

#--discard mortality rates
hm.fxd = setup$hm.fxd; 
hm.trl = setup$hm.trl;

#--path to "current" data
dirCurr = "~/Work/StockAssessments-Crab/Assessments/PIBKC/Data/Current";

#--extract historical bycatch in the groundfish fisheries and calculate discard mortality
#----1. data from CAS database (1990/91-2008/09)
#----2. assumes bycatch of PIBKC occurs only in NMFS stat area 513
#----3. catch biomass in kg
fnHist = file.path(dirCurr,"FromAKFIN.PIBKC.GroundfishBycatchEstimates.CAS.1991-2009.csv");
tfHist = "1991/92-2008/09";
dfrHist = readr::read_csv(fnHist) |> 
            dplyr::select(year=`Crab Year`,gear=Gear,
                          `NMFS stat area`=`Reporting Area Code`,
                          biomass=`Estimate Wt (Sum)`) |> 
            dplyr::mutate(year=as.numeric(stringr::str_sub(year,1,4)),
                          gear=tolower(gear)) |> 
            dplyr::filter(year<2009,`NMFS stat area`==513) |>
            dplyr::select(!`NMFS stat area`) |> 
            dplyr::mutate(mortality=ifelse(gear=="fixed",hm.fxd,hm.trl)*biomass);

#--calculate recent bycatch and discard mortality in the groundfish fisheries
#----1. data from the Catch-In-Areas database (2009/2010+)
#----2. CIA data is expanded to the Pribilof District and may include bycatch outside NMFS area 513
fnCurr     = file.path(dirCurr,"FromAKFIN.PIBKC.GroundfishBycatchEstimates.CIA.2009+.csv");
downloadDate = "July 15, 2024";
tmp<-readr::read_csv(file = fnCurr,skip = 6);
names(tmp)<-tolower(names(tmp));
#----drop columns, revise names
tbl = tmp |> 
        dplyr::select(year=`crab year`,
                      `vessel id`,
                      target=`trip target name`,
                      gearcode=`agency gear code`,
                      `stat area`=`adfg stat area code`,
                      `haul count`=`haul count`,
                      biomass=`estimated crab weight (kg)`,
                      number=`estimated number`) |>
        dplyr::filter(year<assmtYr);#--drop too-recent data
tfCurr = paste0("2009/10-",wtsQMD::crabYear(max(tbl$year)));
rm(tmp);
#substitute gear names
tbl$gear = "";
tbl$gear[tbl$gearcode %in% c("PTR","NPT")]<-"trawl";
tbl$gear[tbl$gearcode %in% c("HAL","POT")]<-"fixed";
#replace target "" with "unknown"
tbl$target[tbl$target==""]<-"unknown";
#----drop crab years before 2009
tbl = tbl |> dplyr::filter(year>=2009);
#----find targets with non-zero bycatch
qry<-"select target,sum(biomass) as biomass 
      from tbl
      group by target";
tmp0<-sqldf::sqldf(qry);
targets<-tmp0$target[tmp0$biomass>0.0];#--targets with PIBKC bycatch
rm(tmp0);
#--add in discard mortality (in t)
tbl$mortality = 0.0;
tbl = rbind(tbl |> dplyr::filter(gear=="fixed") |> dplyr::mutate(mortality=hm.fxd*biomass),
            tbl |> dplyr::filter(gear=="trawl") |> dplyr::mutate(mortality=hm.trl*biomass));
#--save tbl
dfrBycatchDetails = tbl;

#--calculate totals by year and gear to match dfrHist
dfrCurr = dfrBycatchDetails |> 
           dplyr::group_by(year,gear) |>
           dplyr::summarize(biomass=Sum(biomass)) |> 
           dplyr::ungroup() |> 
           dplyr::mutate(mortality=ifelse(gear=="fixed",hm.fxd,hm.trl)*biomass);
tmp = dfrCurr |> 
        dplyr::group_by(year) |> 
        dplyr::summarize(biomass=Sum(biomass),
                         mortality=Sum(mortality));
maxByc = max(tmp$biomass);
maxBycYr = tmp$year[tmp$biomass==maxByc];
maxMrt = max(tmp$mortality);
maxMrtYr = tmp$year[tmp$mortality==maxMrt];
dfrLast5Avg = dfrCurr |>
                dplyr::filter(year %in% wtsQMD::last(sort(unique(dfrCurr$year)),5)) |> 
                dplyr::group_by(year) |>                        #--need to sum over gear here before averaging
                dplyr::summarize(catch=sum(biomass),
                                 mortality=sum(mortality)) |> 
                dplyr::ungroup() |> 
                dplyr::summarize(catch=mean(catch),
                                 mortality=mean(mortality));

#--combine historical (from CAS) and current (from CIA) tables for annual mortality
dfrGFM = dplyr::bind_rows(dfrHist,dfrCurr) |> 
          dplyr::group_by(year) |> 
          dplyr::summarize(mortality=Sum(mortality)) |>
          dplyr::ungroup();

#----calculate bycatch by gear type for tables
cap<-"Bycatch of PIBKC in the groundfish fisheries, by gear type. Biomass is in kilograms.";
tmp0 = dfrBycatchDetails |> 
         dplyr::group_by(year,gear) |> 
         dplyr::summarize(vessels=length(unique(`vessel id`)),
                          conf=vessels<3,
                          number=Sum(number),
                          biomass=Sum(biomass),
                          mortality=Sum(mortality)) |> 
         dplyr::ungroup();
dfrLast5AvgByGearType = tmp0 |> 
                          dplyr::filter(year %in% wtsQMD::last(sort(unique(tmp0$year)),5)) |> 
                          dplyr::group_by(gear) |> 
                          dplyr::summarize(catch=mean(biomass),
                                           mortality=mean(mortality)) |>
                          dplyr::ungroup();
dfrBycatchByGearType = dplyr::bind_rows(
                         tmp0,
                         dfrHist |> dplyr::mutate(vessels=NA,
                                                  conf=FALSE,
                                                  number=NA)) |>
                         dplyr::arrange(year,gear);
rm(tmp0);

#--calculate bycatch by target fishery for tables
#--keep targets with overall bycatch at least > 10 kg
tmp0 = tbl |> 
         dplyr::group_by(year,target) |> 
         dplyr::summarize(vessels=length(unique(`vessel id`)),
                          conf=vessels<3,
                          number=Sum(number),
                          biomass=Sum(biomass),
                          mortality=Sum(mortality)) |> 
         dplyr::ungroup();
tmp1 = tmp0 |> 
         dplyr::group_by(target) |> 
         dplyr::summarize(biomass=Sum(biomass)) |> 
         dplyr::ungroup() |> 
         dplyr::filter(biomass>=10);
tmp0 = tmp0 |> dplyr::filter(target %in% tmp1$target);
rm(tmp1);

dfrBycatchByTargetType  = tmp0;
dfrLast5AvgByTargetType = tmp0 |> 
                           dplyr::filter(year %in% wtsQMD::last(sort(unique(tmp0$year)),5)) |> 
                           dplyr::group_by(target) |> 
                           dplyr::summarize(biomass=mean(biomass),
                                            mortality=mean(mortality)) |>
                           dplyr::ungroup();
rm(tmp0);

#--list for output results
gfr = list(
  units_biomass="kg",
  units_abundance="ones",
  assmtYr=assmtYr,
  strYr=strYr,
  tfHist=tfHist,
  tfCurr=tfCurr,
  timeFrame=paste0(stringr::str_sub(tfHist,1,7),"-",stringr::str_sub(tfCurr,-7,-1)),
  fnHist=fnHist,
  fnCurr=fnCurr,
  downloadDate=downloadDate,
  gisDir = gisDir,
  hm.fxd=hm.fxd,
  hm.trl=hm.trl,
  targets=targets,
  dfrBycatchDetails=dfrBycatchDetails,
  dfrGFM=dfrGFM,
  maxByc=maxByc,
  maxBycYr=maxBycYr,
  maxMrt=maxMrt,
  maxMrtYr=maxMrtYr,
  avgByc5=dfrLast5Avg$catch,
  avgMrt5=dfrLast5Avg$mortality,
  dfrBycatchByGearType=dfrBycatchByGearType,
  dfrLast5AvgByGearType=dfrLast5AvgByGearType,
  dfrBycatchByTargetType=dfrBycatchByTargetType,
  dfrLast5AvgByTargetType=dfrLast5AvgByTargetType
);

wtsUtilities::saveObj(gfr,fn=file.path(dirOut,"rda_GroundfishFisheriesResults.RData"));

#--clean up
rm(assmtYr,cap,dfrBycatchByGearType,dfrBycatchByTargetType,
   dfrBycatchDetails,dfrCurr,dfrGFM,dfrHist,
   dfrLast5Avg,dfrLast5AvgByGearType,dfrLast5AvgByTargetType,
   dirCurr,dirOut,dirPrj,downloadDate,fnCurr,fnHist,
   gisDir,hm.fxd,hm.trl,maxByc,maxBycYr,maxMrt,maxMrtYr,
   qry,strYr,Sum,targets,tbl,tfCurr,tfHist,tmp,verbosity);
if (!rstudio) rm(gfr); rm(rstudio);

