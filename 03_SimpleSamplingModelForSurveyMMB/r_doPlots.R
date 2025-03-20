#--create figures for RW model output (just to check)
library(ggplot2)
library(knitr)

dirPrj = rstudioapi::getActiveProject();
#--if file is being run from RStudio, set dirOut to path to file for output
rstudio = rstudioapi::isAvailable();#--is file run from RStudio?
if (rstudio){
  dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
  dirInp = dirThs;
  dirOut = dirThs;
  testing = TRUE;
} else {
  testing = FALSE;
}
#--otherwise dirOut should have been already set
if (!exists("dirInp")) stop("RStudio not running and dirInp does not exist");
if (!exists("dirOut")) stop("RStudio not running and dirOut does not exist");

#--read in results from r_doCalcs
rr = wtsUtilities::getObj(file.path(dirInp,"rda_RWModelForSurveyMMB.RData"));

lstAll = list(); #--list for output
lstAll

old_thm = ggplot2::theme_set(cowplot::theme_cowplot(font_size = 10) +
                             cowplot::background_grid() +
                             cowplot::panel_border());
thm = wtsPlots::getStdTheme();

#--create mcmc diagnostic plots for parameters and derived quantities
p1 = wtsMCMC::mcmcTMB_PlotParamDiagnostics(rr$mcmc,"log_PE");
p2 = wtsMCMC::mcmcTMB_PlotParamDiagnostics(rr$mcmc,paste0("log_biomass_pred[",rr$tdy_mxIdx,"]"),label="ln-scale terminal MMB");
p3 = wtsMCMC::mcmcTMB_PlotParamDiagnostics(rr$dfrPDs,paste0(  "biomass_pred[",rr$tdy_mxIdx,"]"),label="terminal MMB");
if (testing){
  print(p1); print(p2); print(p3);
}

#--plot the rema time series results
plts = rema::plot_rema(rr$tdy_rema);
#----observed + predicted
p = plts$biomass_by_strata + 
        geom_ribbon(data=rr$dfrSrvMMB_last,mapping=aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci),
                    colour=NA,fill="red",alpha=0.5,inherit.aes=FALSE) + 
        geom_line(data=rr$dfrSrvMMB_last,mapping=aes(x=year,y=pred),colour="red",inherit.aes=FALSE) + 
        thm;       
if (testing) print(p + ylab("Biomass (t)"));
if (testing) print(p + ylab("Biomass (t)") + scale_y_log10());
p = plts$total_predicted_biomass + 
        geom_ribbon(data=rr$dfrSrvMMB_last,mapping=aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci),
                    colour=NA,fill="red",alpha=0.5,inherit.aes=FALSE) + 
        geom_line(data=rr$dfrSrvMMB_last,mapping=aes(x=year,y=pred),colour="red",inherit.aes=FALSE) + 
        thm; #--predicted only
p1 = p + ylab("Biomass (t)");
p2 = p + ylab("Biomass (t)") + scale_y_log10();
pg = cowplot::plot_grid(p1,p2,ncol=1);

#--plot one-step-ahead residuals diagnostics
p1 = rr$osas$plots$biomass_resids + thm;
p2 = rr$osas$plots$biomass_fitted + thm;
p3 = rr$osas$plots$biomass_hist   + thm;
p4 = rr$osas$plots$biomass_qqplot + thm;
pg1 = cowplot::plot_grid(p1,p2,nrow=1);
pg2 = cowplot::plot_grid(p3,p4,nrow=1);
pg = cowplot::plot_grid(pg1,pg2,ncol=1);
if (testing) print(pg);

ggplot2::theme_set(old_thm);
