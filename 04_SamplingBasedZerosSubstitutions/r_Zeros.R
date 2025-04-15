#--dealing with 0's in mature male biomass time series
require(ggplot2);

rstudio = rstudioapi::isAvailable();#--is file run from RStudio?
if (!rstudio) stop("r_doCalcs.R for survey MMB random walk model must be run/sourced from RStudio.");

dirPrj = rstudioapi::getActiveProject();
setup  = wtsUtilities::getObj(file.path(dirPrj,"rda_AssessmentSetup.RData"));

#--set up file paths----
##--if file is being run from RStudio, set dirOut to path to file for output
if (rstudio){
  testing = TRUE;
  dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
  dirInp = dirThs;
  dirOut = dirThs;
} else {
  testing = FALSE;
}
##--otherwise dirOut should have been already set
if (!exists("dirInp")) stop("RStudio not running and dirInp does not exist");
if (!exists("dirOut")) stop("RStudio not running and dirOut does not exist");

#--load survey data (input biomass units: 1000's t)----
lstSrv = wtsUtilities::getObj(setup$fnData_TSs);
dfr = lstSrv$dfrACD.ByX.PD |>
        dplyr::filter(SEX=="mature males") |>
        dplyr::mutate(year=YEAR,
                      area=STRATUM_AREA,
                      numHls=numHauls,
                      numNon0=numNonZeroHauls,
                      totAbd=1000*totABUNDANCE,       #--thousands of individuals
                      cvAbd=cvABUNDANCE,
                      totBio=1000*totBIOMASS,         #--in t
                      cvBio=cvBIOMASS,
                      mnWgt=(totBIOMASS/totABUNDANCE)/1000,  #--in t
                      .keep="none");#--keep only columns created here

##--plot number of non-0 hauls by year----
ggplot(dfr,aes(x=year)) +
  geom_line(aes(y=numNon0),colour="blue") + geom_point(aes(y=numNon0),colour="blue") +
  scale_y_log10() +
  labs(x="year",y="number of non-0 hauls") +
  wtsPlots::getStdTheme();

##--plot number of non-0 hauls by total abundance----
ggplot(dfr,aes(x=numNon0)) +
  geom_path(aes(y=totAbd,colour=year)) + geom_point(aes(y=totAbd,colour=year)) +
  geom_smooth(aes(y=totAbd),alpha=0.2) +
  scale_y_log10() +
  scale_colour_viridis_c(option="magma") +
  labs(x="number of non-0 hauls",y="total abundance (1,000's)") +
  wtsPlots::getStdTheme();

#--calculate median mean mature male weight----
mmWgt = median(dfr$mnWgt,na.rm=TRUE);#--in t

#--assumed area swept----
As = 0.012; # square nmi

#--define function to calculate density that yields Pr(N=0) = p given area swept and K hauls
calcDensity<-function(p,K,As){
  d = -(1/(K*As))*log(p);
  return(d);
}
ds = tidyr::expand_grid(p=seq(0.01,0.99,0.01),K=1:100) |>
     dplyr::mutate(d=calcDensity(p,K,As));
ggplot(ds,aes(x=d,y=p,colour=K,group=K)) +
  geom_line() +
  scale_color_viridis_c(option="magma") +
  labs(x="density",y="Pr(N=0|K hauls)",colour="number\nhauls");

#--calculate densities and equivalent total abundance and biomass----
##--use sampled hauls for density calculation----
dfp = dfr |> dplyr::mutate(d20=calcDensity(0.2,numHls,As), #--density in #/sq nmi
                           d50=calcDensity(0.5,numHls,As),
                           d80=calcDensity(0.8,numHls,As),
                           a20=area*d20/1000,  #--abundance (in 1,000's) expanded to total area
                           a50=area*d50/1000,
                           a80=area*d80/1000,
                           b20=area*d20*mmWgt, #--biomass (in t) expanded to total area
                           b50=area*d50*mmWgt,
                           b80=area*d80*mmWgt
                           );

ggplot(dfp,aes(x=year)) +
  geom_line(aes(y=totAbd),colour="blue") + geom_point(aes(y=totAbd),colour="blue") +
  geom_ribbon(aes(ymin=a80,ymax=a20),alpha=0.2) +
  geom_line(aes(y=a50)) +
  scale_y_log10() +
  labs(x="year",y="total abundance (1,000's)") +
  wtsPlots::getStdTheme();

ggplot(dfp,aes(x=year)) +
  geom_line(aes(y=totBio),colour="blue") + geom_point(aes(y=totBio),colour="blue") +
  geom_ribbon(aes(ymin=b80,ymax=b20),alpha=0.2) +
  geom_line(aes(y=b50)) +
  scale_y_log10() +
  labs(x="year",y="total biomass (t)") +
  wtsPlots::getStdTheme();

##--use median number of non-zero hauls for density calculation----
mdHls = median(dfr$numNon0);
dfq = dfr |> dplyr::mutate(d20=calcDensity(0.2,mdHls,As), #--density in #/sq nmi
                           d50=calcDensity(0.5,mdHls,As),
                           d80=calcDensity(0.8,mdHls,As),
                           a20=area*d20/1000,  #--abundance (in 1,000's) expanded to total area
                           a50=area*d50/1000,
                           a80=area*d80/1000,
                           b20=area*d20*mmWgt, #--biomass (in t) expanded to total area
                           b50=area*d50*mmWgt,
                           b80=area*d80*mmWgt
                           );

ggplot(dfq,aes(x=year)) +
  geom_line(aes(y=totAbd),colour="blue") + geom_point(aes(y=totAbd),colour="blue") +
  geom_ribbon(aes(ymin=a80,ymax=a20),alpha=0.2) +
  geom_line(aes(y=a50)) +
  scale_y_log10() +
  labs(x="year",y="total abundance (1,000's)") +
  wtsPlots::getStdTheme();

ggplot(dfq,aes(x=year)) +
  geom_line(aes(y=totBio),colour="blue") + geom_point(aes(y=totBio),colour="blue") +
  geom_ribbon(aes(ymin=b80,ymax=b20),alpha=0.2) +
  geom_line(aes(y=b50)) +
  scale_y_log10() +
  labs(x="year",y="total biomass (t)") +
  wtsPlots::getStdTheme();

