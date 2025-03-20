#--create tables (just to check)
require(tables)
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
tdy      = rr$tdy_rema;      #--tidy version of rema results
dfr_last = rr$dfrSrvMMB_last;#--results from last assessment

#--create table for parameter estimate
prs = tdy$parameter_estimates;
kbl = prs |> kableExtra::kbl(booktabs=TRUE,digits=4) |>
             kableExtra::kable_styling(bootstrap_options="bordered",
                                       latex_options="HOLD_position")
kbl

#--create table for time series comparison
curr = as.character(rr$assmt_yr);
tbl = dplyr::bind_rows(
        tdy$biomass_by_strata |> dplyr::mutate(type="observed") |> 
          dplyr::select(type,year,value=obs,lci=obs_lci,uci=obs_uci),
        dfr_last |> dplyr::mutate(type="base") |> 
           dplyr::select(type,year,value=pred,lci=pred_lci,uci=pred_uci),
        tdy$biomass_by_strata |> dplyr::mutate(type=curr) |> 
          dplyr::select(type,year,value=pred,lci=pred_lci,uci=pred_uci)
      ) |> 
      dplyr::mutate(type=factor(type,levels=c("observed","base",curr)))
tblr = tabular(Factor(year)~DropEmpty("--","cell") * Format(digits=2) * 
                             (Heading("estimate") * Factor(type) * value*mean + 
                              Heading("lci")   * Factor(type) * lci*mean + 
                              Heading("uci")   * Factor(type) * uci*mean),
               data=tbl);
colLabels(tblr) = colLabels(tblr)[c(3,2),];
colLabels(tblr)[1,] = c("","value","","","lci","","","uci","");
kbl = tblr |> tables::toKable(format="html",booktabs=TRUE) |> 
              kableExtra::kable_classic(html_font="Cambria") |> 
              kableExtra::column_spec(c(1,4,7),border_right=TRUE); 
kbl
