#--get retained catch and bycatch in the crab fisheries
#----1. retained catch data starts in 1973/74
#----2. bycatch data starts in 1996/97 (females, legal males, sublegal males)
#--NOTE: input units for biomass are t, for abundance are ones
#--if file is being run from RStudio, set dirOut to path to file for output
rstudio = rstudioapi::isAvailable();
if (!rstudio) 
  stop("Must run/source r_doCalcs.R from the RStudio console.");

dirPrj = rstudioapi::getActiveProject();
setup = wtsUtilities::getObj(file.path(dirPrj,"rda_AssessmentSetup.RData"));

if (rstudio){
  dirOut = dirname(rstudioapi::getActiveDocumentContext()$path);
}
#--otherwise dirOut should have been already set
if (!exists("dirOut")) stop("RStudio not running and dirOut does not exist");

assYr = setup$asmtYr; #assessment year

#--get input files (must update these beforehand with info from previous fishery year)
dirInp = "~/Work/StockAssessments-Crab/Assessments/PIBKC/Data/Current";
dataFileRetained = file.path(dirInp,"data.RetainedCatch.csv");
dataFileBycatch  = file.path(dirInp,"data.BycatchInCrabFisheries.csv");
report_date      = "July 6, 2023";#--report from Ben Daly: 0 crab taken as retained or bycatch (TODO: move to r_doAssessmentSetup.R)
#--copy files to dirOut
file.copy(from=dataFileRetained,to=file.path(dirOut,basename(dataFileRetained)),overwrite=TRUE);
file.copy(from=dataFileBycatch, to=file.path(dirOut,basename(dataFileBycatch)), overwrite=TRUE);

strYr = wtsQMD::crabYear(assYr-1);#--final fishery year

hm.pot = setup$hm.pot;

outDir = dirOut;
if (!dir.exists(outDir)) dir.create(outDir);

#--read in retained catch and average CPUE (num. legal crabs/pot)
#----assumes no retained catch 1999/00+
tmp    = readr::read_csv(file=dataFileRetained,skip=2);
dfrRet = tmp |> 
           dplyr::mutate(`crab year`=year,year=as.numeric(substr(`crab year`,1,4))) |>
           dplyr::filter(!dplyr::if_all(.fns=is.na));
retMax   = max(dfrRet$biomass);
retMaxYr = dfrRet$`crab year`[dfrRet$biomass==retMax];
rm(tmp);

#--read in bycatch data
tmp    = readr::read_csv(file=dataFileBycatch,skip=1);
dfrByc = tmp |> 
           dplyr::mutate(`crab year`=year,year=as.numeric(substr(`crab year`,1,4))) |> 
           dplyr::filter(year<assYr) |> 
           dplyr::mutate(total=females+`legal males`+`sublegal males`,
                         mortality = hm.pot*total);
nyrs   = nrow(dfrByc);
bycMaxM  = max(dfrByc$mortality);
bycMaxYr = dfrByc$`crab year`[dfrByc$mortality==bycMaxM];
bycAvgM  = mean(dfrByc$mortality[nyrs+(-4:0)])
bycCurM  = dfrByc$mortality[nyrs];
rm(tmp);

#--save
cfr = list(assYr=assYr,
           units_biomass="t",
           units_abundance="ones",
           hm.pot=hm.pot,
           dataFileRetained=dataFileRetained,
           dataFileBycatch=dataFileBycatch,
           report_date=report_date,
           strYr=strYr,
           outDir=outDir,
           dfrRet=dfrRet,
           retMax=retMax,
           retMaxYr=retMaxYr,
           dfrByc=dfrByc,
           bycMaxM=bycMaxM,
           bycMaxYr=bycMaxYr,
           bycAvgM=bycAvgM,
           bycCurM=bycCurM
           );
wtsUtilities::saveObj(cfr,fn=setup$fnData_CFs);

#--clean up
rm(setup,assYr,hm.pot,dataFileRetained,dataFileBycatch,report_date,strYr,outDir,
   dfrRet,retMax,retMaxYr,dfrByc,bycMaxM,bycMaxYr,bycAvgM,bycCurM,
   dirInp,dirOut,nyrs);
if (!rstudio) rm(cfr); rm(rstudio);
