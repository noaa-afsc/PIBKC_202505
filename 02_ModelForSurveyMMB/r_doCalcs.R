#--Run rema model on mature male  biomass time series
#----to reduce variance. Results in "rda_REMA_Results.RData".

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

#--get assessment year----
assmt_yr = setup$asmtYr;

#--create list for output----
lstAll = list();
lstAll$assmt_yr = assmt_yr;

#--read in results from previous (2023) assessment (will need to change in future)----
##--TODO: need to redefine this for 2023 results, not 2021 results----
ssr_last = wtsUtilities::getObj(setup$fnLastAsmt.RWMforSrvMMB);
lstAll$dfrSrvMMB_last = ssr_last$dfrPrds |>
                          dplyr::filter(option=="0's as NAs") |>     #--selected model name
                          dplyr::select(year,pred,pred_lci,pred_uci);

#--load survey data (input biomass units: 1000's t)----
lstSrv = wtsUtilities::getObj(setup$fnData_TSs);
dfr = lstSrv$dfrACD.ByX.PD |>
        dplyr::filter(SEX=="mature males") |>
        dplyr::mutate(totBIOMASS=1000*totBIOMASS) |>  #--convert to t
        dplyr::select(year=YEAR,
                      value=totBIOMASS,
                      cv=cvBIOMASS);

#--2023, 2024: survey MMB time series now has observed zeros!----
##----3 options considered:----
###----1. 0's as NA's (GPT's current choice)----
###----2. small value, large cv----
###----3. use Tweedie distribution----

#--define function to plot one-step-ahead residuals----
plotOSAs<-function(osas){
    old_thm = ggplot2::theme_set(cowplot::theme_cowplot(font_size = 10) +
                                 cowplot::background_grid() +
                                 cowplot::panel_border());
    thm = wtsPlots::getStdTheme();
    p1 = osas$plots$biomass_resids + thm;
    p2 = osas$plots$biomass_fitted + thm;
    p3 = osas$plots$biomass_hist + thm;
    p4 = osas$plots$biomass_qqplot + thm;
    pg1 = cowplot::plot_grid(p1,p2,nrow=1);
    pg2 = cowplot::plot_grid(p3,p4,nrow=1);
    pg  = cowplot::plot_grid(pg1,pg2,ncol=1);
    ggplot2::theme_set(old_thm);
    return(pg);
}
lstAll$plotOSAs = plotOSAs;

#--Define function to run rema model with zeros option and generate results and diagnostics as a list----
runREMA<-function(model_name,dfr,zeros,doMCMC=TRUE){
  if (testing) cat("Running rema model",model_name,"\n")
  resLst = list();
  #--set cv's to NA where value = 0
  dfr = dfr |> dplyr::mutate(cv=ifelse(value==0,NA,cv));
  #--prepare rema input
  inp = rema::prepare_rema_input(
                model_name="rema",
                end_year=assmt_yr,
                biomass_dat=tibble::tibble(strata="Pribilof Islands",
                                           year=dfr$year,
                                           biomass=dfr$value,
                                           cv=dfr$cv),
                zeros = zeros)
  resLst$rema_inp = inp;

  #--fit the data
  mdl = rema::fit_rema(inp,
                       n.newton=1,
                       do.sdrep=TRUE,
                       do.check=TRUE,
                       do.fit=TRUE,
                       save.sdrep=TRUE,
                       MakeADFun.silent=FALSE);
  resLst$rema_mdl = mdl;

  #--tidy up the model results
  tdy = rema::tidy_rema(mdl,alpha=1-setup$ci);
  resLst$tdy_alpha = 1-setup$ci;
  resLst$tdy_rema  = tdy;
  resLst$tdy_mxIdx = length(tdy$biomass_by_strata$year);
  resLst$tdy_mxYr  = tdy$biomass_by_strata$year[resLst$tdy_mxIdx];
  resLst$dfrSrvMMB = tdy$biomass_by_strata |>
                                     dplyr::mutate(pred_cv=sqrt(exp(sd_log_pred^2)-1)) |>
                                     dplyr::select(year,
                                                   obs,obs_cv,obs_lci,obs_uci,
                                                   pred,pred_cv,pred_lci,pred_uci) |>
                                     dplyr::mutate(pred=unname(pred),
                                                   pred_cv=unname(pred_cv),
                                                   pred_lci=unname(pred_lci),
                                                   pred_uci=unname(pred_uci));

  if (testing) print(tdy$parameter_estimates);

  #--get one-step-ahead residuals
  osas = rema::get_osa_residuals(mdl);
  resLst$osas = osas;
  if (testing) print(plotOSAs(osas));

  #--do MCMC analysis
  if (doMCMC){
    options(mc.cores = parallel::detectCores());
    mcmc = tmbstan::tmbstan(mdl,algorithm="NUTS",laplace=FALSE,chains=8,warmup=5000,iter=25000,thin=4);
    resLst$mcmc_params = list(algorithm="NUTS",laplace=FALSE,chains=8,warmup=5000,iter=25000,thin=4);
    resLst$mcmc = mcmc;
    #----extract diagnostics
    mon = rstan::monitor(mcmc);
    rHat = max(mon$Rhat);    #--potential scale reduction factor: value < 1.05 indicates convergence
    ESS  = min(mon$Tail_ESS);#--minimum effective sample size (should be > 100)
    resLst$mcmc_rHat = rHat;
    resLst$mcmc_ESS  = ESS;
    if (testing) {
      cat("rHat =",rHat,"\n");
      cat("ESS =",ESS,"\n");
    }
    #----extract posterior densities for derived quantities
    dfrPDs = wtsMCMC::mcmcTMB_RunReport(mdl,mcmc,vars="biomass_pred");
    resLst$dfrPDs = dfrPDs;
    if (testing) {
      p = wtsMCMC::mcmcTMB_PlotParamDiagnostics(dfrPDs,
                                                paste0("biomass_pred[",resLst$tdy_mxIdx,"]"),
                                                label="terminal MMB");
      print(p)
    }
  }#--doMCMC
  return(resLst);
}

#--list of rema models and results----
lstAll$remas = list();
zeros1 = list(assumption = 'NA');
zeros2 = list(assumption = 'small_constant',
              options_small_constant = c(0.01, 1.5)); # values: 1) small constant, 2) assumed CV
zeros3 = list(assumption = 'tweedie',
              options_tweedie=list(zeros_cv=1.5));

opt1 = "0's as NAs";
opt2 = "small constant";
opt3 = "Tweedie";
lstAll$remas[[opt1]] = runREMA(opt1,dfr,zeros=zeros1);
lstAll$remas[[opt2]] = runREMA(opt2,dfr,zeros=zeros2);
lstAll$remas[[opt3]] = runREMA(opt3,dfr,zeros=zeros3,doMCMC=FALSE);#--mcmc very slow

#--compare model estimates----
require(ggplot2)
lstTmp = list();
for (opt in names(lstAll$remas)){
  #--testing: opt = names(lstAll$remas)[1];
  rma = lstAll$remas[[opt]];
  lstTmp[[opt]] = rma$dfrSrvMMB |> dplyr::mutate(option=opt);
}
dfrTmp = dplyr::bind_rows(lstTmp); rm(lstTmp);
compareOpts<-function(dfrTmp,ymax=NA,dfrLast=NULL){
  dfrDat = dfrTmp |> dplyr::distinct(year,obs,obs_lci,obs_uci);
  p1 = ggplot(dfrTmp,aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci,colour=option,fill=option)) +
        geom_errorbar(data=dfrDat,mapping=aes(x=year,y=obs,ymin=obs_lci,ymax=obs_uci),inherit.aes=FALSE) +
        geom_point(data=dfrDat,mapping=aes(x=year,y=obs),inherit.aes=FALSE);
  if (!is.null(dfrLast)){
    p1 = p1 +
          geom_ribbon(data=dfrLast,mapping=aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci),
                      colour=NA,fill="red",alpha=0.5,inherit.aes=FALSE) +
          geom_line(data=dfrLast,mapping=aes(x=year,y=pred),colour="red",inherit.aes=FALSE)
  }
  p1 = p1 +
        geom_ribbon(alpha=0.5,colour=NA) + geom_line() +
        scale_y_continuous(limits=c(0,ymax),oob=scales::squish) +
        scale_fill_viridis_d(aesthetics=c("fill","colour")) +
        labs(y="biomass (t)") +
        theme(legend.position=c(0.99,0.99),
              legend.justification=c(1,1)) +
        wtsPlots::getStdTheme();
  p2 = p1 + scale_y_log10();
  pg = cowplot::plot_grid(p1+theme(axis.title.x=element_blank()),
                          p2+theme(legend.position="none"),
                          ncol=1);
  return(pg);
}
lstAll$compareOpts = compareOpts;
lstAll$dfrPrds = dfrTmp;

compareOpts(lstAll$dfrPrds,ymax=50000,
            dfrLast=lstAll$dfrSrvMMB_last);
compareOpts(lstAll$dfrPrds |> dplyr::filter(year>=2010))
compareOpts(lstAll$dfrPrds |> dplyr::filter(year>=2010),
            dfrLast=lstAll$dfrSrvMMB_last |> dplyr::filter(year>=2010));

#--save list object----
wtsUtilities::saveObj(lstAll,setup$fnSrvMMB);
rm(setup,assmt_yr,dfr,dirInp,dirOut,dirPrj,dirThs,
   inp,lstSrv,ssr_last,testing);
rm(compareOpts,opt,opt1,opt2,opt3,plotOSAs,rma,runREMA,zeros1,zeros2,zeros3,dfrTmp)
if (!rstudio) rm(lstAll); rm(rstudio);
