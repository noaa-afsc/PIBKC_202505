#--set up parameters for the assessment (full assessment in 2025)
#----must source this as active document from RStudio
dirPrj = rstudioapi::getActiveProject();
if (dirname(rstudioapi::getActiveDocumentContext()$path)=="") 
  stop("must source r_doAssessment.R as active document!",call.=FALSE);

if (FALSE){
  #--copy latex includes file for qmd (can't figure out how to reference package file in yaml)
  fn = system.file("files/ltx_ExtraLatexIncludesNew.tex",package="wtsQMD");
  file.copy(from=fn,to=file.path(dirPrj,"Text","ltx_ExtraLatexIncludes.tex"),overwrite=TRUE);
}

s = list();
s[["asmtYr"]]     = 2024;#--assessment year
s[["lastAsmtYr"]] = 2023;#--last assessment year
s[["nextAsmtYr"]] = 2025;#--next assessment year
s[["prvOFL"]]   = 1.16;#--OFL (in t) for the year prior to the assessment
s[["prvABC"]]   = 0.87;#--ABC (in t) for the year prior to the assessment
s[["prvBuf"]]   = 0.25;#--ABC buffer for the year prior to the assessment
s[["ci"]]       = 0.80;#--confidence interval for plots, etc.
s[["maleMatZ"]] = 120; #--min mature male size
s[["maleLglZ"]] = 135; #--min legal male size
s[["hm.pot"]]   = 0.2; #--assumed discard mortality in crab pot fisheries
s[["hm.fxd"]]   = 0.2; #--assumed discard mortality in groundfish fixed gear fisheries
s[["hm.trl"]]   = 0.8; #--assumed discard mortality in groundfish trawl gear fisheries
s[["M"]]        = 0.18;#--assumed rate of natural mortality
s[["t.sf"]]     = 3/12;#--assumed time from survey to fishing mortality
s[["t.fm"]]     = 4/12;#--assumed time from fishing mortality to mating
s[["pct.male"]] = 0.5; #--fraction of discard mortality in groundfish fisheries assumed to be MMB
s[["timeFrameForBmsy"]] = "1980:1984,1990:1997"; #--time frame for calculating Bmsy (as string)
s[["yrsForBmsy"]]       = c(1980:1984,1990:1997);#--time frame for calculating Bmsy
s[["nYrsTheta"]]  = 3;   #--number of years to average for theta calculation
s[["gamma"]]      = 1.0; #--control rule parameter
s[["alpha"]]      = 0.1; #--control rule parameter
s[["beta"]]       = 0.25; #--control rule parameter
s[["abcBuf"]]     = 0.25; #--ABC buffer

#--paths to various folders----
s["pthHstRes"]  = file.path(dirPrj,"../HistoricalAssessmentResults_PIBKC");
s["fnLastAsmt.RWMforSrvMMB"] = file.path(dirPrj,"00_LastAsmt","rda_02_RWModelForSurveyMMB.RData");
s["fnLastAsmt.MMBatMating"]  = file.path(dirPrj,"00_LastAsmt","rda_03_MMBatMating.RData");
s["fnLastAsmt.Tier4MQs"]     = file.path(dirPrj,"00_LastAsmt","rda_04_Tier4MQs.RData");
s["fnLastAsmt.Res4SAFE"]     = file.path(dirPrj,"00_LastAsmt","rda_ResultsForSAFE.RData");
s["fnData_CFs"] = file.path(dirPrj,"01_Data/Fisheries.Crab",        "rda_CrabFisheriesResults.RData");
s["fnData_GFs"] = file.path(dirPrj,"01_Data/Fisheries.Groundfish",  "rda_GroundfishFisheriesResults.RData");
s["fnData_TSs"] = file.path(dirPrj,"01_Data/Surveys.NMFS",          "rda_TrawlSurveyResults.RData");
s["fnSrvMMB"]   = file.path(dirPrj,"02_ModelForSurveyMMB",          "rda_02_RWModelForSurveyMMB.RData");
s["fnMMB"]      = file.path(dirPrj,"03_MMB-at-mating",              "rda_03_MMBatMating.RData");
s["fnT4MQs"]    = file.path(dirPrj,"04_Tier4ManagementQuantities",  "rda_04_Tier4MQs.RData");
s["fnRes4SAFE"] = file.path(dirPrj,"Text",                          "rda_ResultsForSAFE.RData");

#--paths to external images----
s[["fnFoflControlRule"]] = file.path(dirPrj,"Text/extern_images/PNGs","FoflControlRule.png");

#--save list----
fn = file.path(dirname(rstudioapi::getActiveDocumentContext()$path),"rda_AssessmentSetup.RData");
wtsUtilities::saveObj(s,fn);

#--clean up
rm(dirPrj,fn,s);
