#--run sdmTMB for PIBKC----
require(ggplot2);
require(sdmTMB)

dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
##--get EBS survey data----
if (FALSE){
  ###--survey years to include----
  yrs = 1975:2024;

  ###--get my PIBKC crab survey results (run r_wts.PIBKC.doCalcs.R)----
  lstWTS = wtsUtilities::getObj(file.path(dirThs,"../02_CrabpackDataComparisons/rda_wts.PIBKC.NMFS.RData"));
  #####--Calculate annual survey abundance and biomass by sex/size group for the Pribilof District
  ######--calculate CPUE by haul and population category (encoded as SEX)
  verbosity=FALSE;
  dfrCPUE.ByXM <-tcsamSurveyData::calcCPUE.ByHaul(lstWTS$dfr.HD,
                                                  lstWTS$dfr.ID,
                                                  bySex=TRUE,
                                                  byMaturity=FALSE,
                                                  byShellCondition=FALSE,
                                                  bySize=FALSE,
                                                  export=FALSE,
                                                  verbosity=verbosity);
  #####--calculate CPUE by station
  dfrCPUE.ByXM <-tcsamSurveyData::calcCPUE.ByStation(lstWTS$dfr.SD,
                                                     dfrCPUE.ByXM,
                                                     export=FALSE,
                                                     verbosity=verbosity);
  #####--calculate biomass/abundance by stratum
  dfrACD.ByXM.ByStrata <-tcsamSurveyData::calcAB.ByStratum(tbl_strata=lstWTS$dfr.SD,
                                                           tbl_cpue=dfrCPUE.ByXM,
                                                           export=FALSE,
                                                           verbosity=verbosity);
  #####--calculate biomass/abundance for the Pribilof District
  dfrACD.ByXM<-tcsamSurveyData::calcAB.EBS(dfrACD.ByXM.ByStrata,
                                            export=FALSE,
                                            verbosity=verbosity) |>
                 dplyr::mutate(STRATUM="Pribilof District");
  dfrACD_MM = dfrACD.ByXM |>
                 dplyr::filter(SEX=="4. mature males"); #--keep mature males only
  ggplot(dfrACD_MM |> dplyr::mutate(p=numNonZeroHauls/numHauls),
         aes(x=totABUNDANCE,y=p)) +
    labs(x="total abundance (millions)",y="proportion of non-zero hauls") +
    geom_point() + wtsPlots::getStdTheme();
  ggplot(dfrACD_MM |> dplyr::mutate(p=numNonZeroHauls/numHauls),
         aes(x=log(totABUNDANCE),y=log(p))) +
    labs(x="total abundance (millions)",y="log(proportion of non-zero hauls)") +
    geom_point() +
    wtsPlots::getStdTheme();
  wtsUtilities::saveObj(dfrACD_MM,file.path(dirThs,"rda_ACD_MM.RData"));

  ###--for CPUE,select only mature males, add environmental covariates,
  ####--convert lon/lat to utm (might want to use Alaska Albers)
  dfrCPUE_MM = dfrCPUE.ByXM |> dplyr::filter(SEX=="4. mature males") |>        #--mature males only
                 dplyr::select(YEAR,GIS_STATION,numIndivs,numCPUE,wgtCPUE) |>  #--keep only
                 dplyr::inner_join(lstWTS$evs.stns |> sf::st_drop_geometry() |>
                                     dplyr::select(YEAR,GIS_STATION,LON,LAT,BOTTOM_DEPTH,BOTTOM_TEMP),
                                   by=dplyr::join_by(YEAR,GIS_STATION)) |>
                 sdmTMB::add_utm_columns(ll_names=c("LON","LAT"),
                                                     utm_names=c("X","Y"),
                                                     units="km");
  wtsUtilities::saveObj(dfrCPUE_MM,file.path(dirThs,"rda_CPUE_MM.RData"));
} else {
  dfrACD_MM  = wtsUtilities::getObj(file.path(dirThs,"rda_ACD_MM.RData"));
  dfrCPUE_MM = wtsUtilities::getObj(file.path(dirThs,"rda_CPUE_MM.RData"));
}
  ##--make meshes----
if (FALSE){
  ###--make inla_mesh using fmesher----
  #dfrXY = tibble::as_tibble(dfrCPUE_MM |> dplyr::filter(YEAR==2023)) |> dplyr::select(X,Y);
  dfrXY = tibble::as_tibble(dfrCPUE_MM) |> dplyr::select(X,Y);
  inla_mesh <- fmesher::fm_mesh_2d_inla(loc = cbind(dfrXY$X, dfrXY$Y), # coordinates
                                        max.edge = c(25, 50), # max triangle edge length; inner and outer meshes
                                        offset = c(5, 25),  # inner and outer border widths
                                        cutoff = 5 # minimum triangle edge length
                                        );
  ggplot() +
    fmesher::geom_fm(data = inla_mesh) +
    labs(subtitle="inla_mesh") +
    geom_point(data=dfrCPUE_MM,mapping=aes(x=X,y=Y)) + wtsPlots::getStdTheme();
  ###--make sdmTMB mesh based on inla_mesh----
  mesh1 = sdmTMB::make_mesh(dfrCPUE_MM,c("X","Y"),mesh=inla_mesh);
  ggplot() + inlabru::gg(mesh1$mesh) +
    labs(subtitle="sdmTMB inla_mesh") +
    geom_point(data=dfrCPUE_MM,mapping=aes(x=X,y=Y)) + wtsPlots::getStdTheme();
  ###--make sdmTMB mesh based on "kmeans" option with 40 knots----
  mesh2 = sdmTMB::make_mesh(dfrCPUE_MM,xy_cols=c("X","Y"),type="kmeans",n_knots=40);
  ggplot() + inlabru::gg(mesh2$mesh) +
    labs(subtitle="kmeans (40 knots)") +
    geom_point(data=dfrCPUE_MM,mapping=aes(x=X,y=Y)) + wtsPlots::getStdTheme();
  ###--make sdmTMB mesh based on cutoff=5----
  mesh3<-sdmTMB::make_mesh(dfrCPUE_MM,c("X","Y"),cutoff=5,type="cutoff");
  ggplot() + inlabru::gg(mesh3$mesh) +
    labs(subtitle="cutoff=5") +
    geom_point(data=dfrCPUE_MM,mapping=aes(x=X,y=Y)) + wtsPlots::getStdTheme();
  ###--compare all meshes----
  ggplot() + inlabru::gg(mesh3$mesh) + inlabru::gg(mesh2$mesh) + inlabru::gg(mesh1$mesh) +
    labs(subtitle="all meshes") +
    geom_point(data=dfrCPUE_MM,mapping=aes(x=X,y=Y)) + wtsPlots::getStdTheme();
  wtsUtilities::saveObj(list(mesh1=mesh1,mesh2=mesh2,mesh3=mesh3),
                        file.path(dirThs,"rda_meshes.RData"));
} else {
  mesh3 = wtsUtilities::getObj(file.path(dirThs,"rda_meshes.RData"))$mesh3;
}
  ##--fit CPUE using mesh 3----
if (FALSE){
  dfr = dfrCPUE_MM |> dplyr::mutate(log_depth=log(BOTTOM_DEPTH),
                                    log_temp = log(BOTTOM_TEMP)) |>
          dplyr::select(year=YEAR,X,Y,wgtCPUE,log_depth,log_temp);
  ###--model m3----
  mdl_m3 <- sdmTMB::sdmTMB(
              wgtCPUE ~ s(log_depth),
              data = dfr,
              mesh = mesh3,
              time="year",
              extra_time = c(2020),
              family = sdmTMB::delta_gamma(type = "poisson-link"),
              spatial = "on",
              spatiotemporal="rw",
              anisotropy = TRUE,
              reml = TRUE,
              offset = NULL,
              do_fit=TRUE,
              do_index=FALSE,
              control=sdmTMB::sdmTMBcontrol(nlminb_loops=4,
                                            newton_loops=4),
              silent = FALSE
            );
  ###--model m3_tw----
  mdl_tw <- sdmTMB::sdmTMB(
              wgtCPUE ~ s(log_depth),
              data = dfr,
              mesh = mesh3,
              time="year",
              extra_time = c(2020),
              family = sdmTMB::tweedie(),
              spatial = "on",
              spatiotemporal="rw",
              anisotropy = TRUE,
              reml = TRUE,
              offset = NULL,
              do_fit=TRUE,
              do_index=FALSE,
              control=sdmTMB::sdmTMBcontrol(nlminb_loops=4,
                                            newton_loops=4),
              silent = FALSE
            );
  mdls = list(mdl_m3=mdl_m3,mdl_tw=mdl_tw);
  wtsUtilities::saveObj(mdls,file.path(dirThs,"rda_mdls_log_depth-st_rw-dg.RData"));
} else {
  mdls = wtsUtilities::getObj(file.path(dirThs,"rda_mdls_log_depth-st_rw-dg.RData"));
}
if (FALSE){
  dfr = dfrCPUE_MM |> dplyr::mutate(log_depth=log(BOTTOM_DEPTH),
                                    log_temp = log(BOTTOM_TEMP)) |>
          dplyr::select(year=YEAR,X,Y,wgtCPUE,log_depth,log_temp);
  ###--model m3_ar----
  mdl_m3_ar = sdmTMB::sdmTMB(
                wgtCPUE ~ s(log_depth),
                data = dfr,
                mesh = mesh3,
                time="year",
                extra_time = c(2020),
                family = sdmTMB::delta_gamma(type = "poisson-link"),
                spatial = "on",
                spatiotemporal="ar1",
                anisotropy = TRUE,
                reml = TRUE,
                offset = NULL,
                do_fit=TRUE,
                do_index=FALSE,
                control=sdmTMB::sdmTMBcontrol(nlminb_loops=4,
                                              newton_loops=4),
                silent = FALSE
              );
  ###--model tw_ar----
  mdl_tw_ar = sdmTMB::sdmTMB(
                wgtCPUE ~ s(log_depth),
                data = dfr,
                mesh = mesh3,
                time="year",
                extra_time = c(2020),
                family = sdmTMB::tweedie(),
                spatial = "on",
                spatiotemporal="ar1",
                anisotropy = TRUE,
                reml = TRUE,
                offset = NULL,
                do_fit=TRUE,
                do_index=FALSE,
                control=sdmTMB::sdmTMBcontrol(nlminb_loops=4,
                                              newton_loops=4),
                silent = FALSE
              );
  mdls_ar = list(mdl_m3_ar=mdl_m3_ar,mdl_tw_ar=mdl_tw_ar);
  wtsUtilities::saveObj(mdls_ar,file.path(dirThs,"rda_mdls_ar_log_depth-st_rw-dg.RData"));
} else {
  mdls_ar = wtsUtilities::getObj(file.path(dirThs,"rda_mdls_ar_log_depth-st_rw-dg.RData"));
}

if (FALSE){
  mdls = c(mdls,mdls_ar);
  wtsUtilities::saveObj(mdls,file.path(dirThs,"rda_mdls-all_log_depth-st_rw-dg.RData"));
} else {
  mdls = wtsUtilities::getObj(file.path(dirThs,"rda_mdls-all_log_depth-st_rw-dg.RData"));
}

  ##--print results----
  for (mdl in names(mdls)){
    cat("###--Results for",mdl,"----\n");
    print(mdls[[mdl]]);
  }

#--compare AICs----
lst = list();
for (mdl in names(mdls)){
  lst[[mdl]] = tibble::tibble(mdl=mdl,AIC=AIC(mdls[[mdl]]));
}
dfrAICs = dplyr::bind_rows(lst);

  ##--check parameters----
  for (mdl in names(mdls)){
    cat("Results for",mdl,"\n");
    sdmTMB::tidy(mdls[[mdl]], conf.int = TRUE)
    sdmTMB::tidy(mdls[[mdl]], "ran_pars", conf.int = TRUE)
  }

  ##--get prediction objects for mesh 3
  lstPrds = list();
  for (mdl in names(mdls)){
    lstPrds[[mdl]] = predict(mdls[[mdl]],
                             newdata=NULL,
                             type="response",
                             se_fit=FALSE,
                             return_tmb_object=TRUE);
  }
  wtsUtilities::saveObj(lstPrds,file.path(dirThs,"rda_prds_log_depth-st_rw-dg.RData"));

  ##--calculate indices and save objects
  lstIdxs = list();
  for (mdl in names(mdls)){
    lstIdxs[[mdl]] = get_index(lstPrds[[mdl]],bias_correct=TRUE);
  }
  wtsUtilities::saveObj(lstIdxs,file.path(dirThs,"rda_idxs_log_depth-st_rw-dg.RData"));

  ##--extract and plot sdmTMB indices----
  lstIdxs1 = list();
  for (mdl in names(mdls)){
    lstIdxs1[[mdl]] = lstIdxs[[mdl]] |>
                        dplyr::select(year,pred=est,pred_lci=lwr,pred_uci=upr,se) |>
                        dplyr::mutate(option=mdl,
                                      pred_cv=se/pred) |>
                        dplyr::select(!se);
  }
  dfrIdxs = dplyr::bind_rows(lstIdxs1);
  ggplot(dfrIdxs,aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci,colour=option,fill=option)) +
    geom_ribbon(alpha=0.3) + geom_point() + geom_line() +
    geom_point(data=dfrACD_MM,mapping=aes(x=YEAR,y=3.4299*totBIOMASS),colour="blue",inherit.aes=FALSE) +
    labs(y="biomass (1,000's t)") +
    wtsPlots::getStdTheme();
  ggplot(dfrIdxs,aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci,colour=option,fill=option)) +
    geom_ribbon(colour=NA,alpha=0.5) + geom_line() +
    geom_vline(xintercept=2020) +
    geom_vline(xintercept=2023:2024,linetype=3) +
    geom_point(data=dfrACD_MM,mapping=aes(x=YEAR,y=3.4299*totBIOMASS),colour="blue",inherit.aes=FALSE) +
    scale_y_log10() +
    labs(y="biomass (1,000's t)") +
    wtsPlots::getStdTheme();

  ##--plot comparisons with REMA models----
  lstREMAs = wtsUtilities::getObj(file.path("../03_SimpleSamplingModelForSurveyMMB/rda_Zeros_REMA_Models.RData"));
  options = c(unique(lstREMAs$dfrPrds$option)[3:5],unique(dfrIdxs$option));
  dfrIdxs1 = lstREMAs$dfrPrds |>
              dplyr::filter(option %in% options) |>
              dplyr::mutate(obs=3.4299*obs/1000,
                           obs_lci=3.4299*obs_lci/1000,
                           obs_uci=3.4299*obs_uci/1000,
                           pred=3.4299*pred/1000,
                           pred_lci=3.4299*pred_lci/1000,
                           pred_uci=3.4299*pred_uci/1000
                           ) |>
              dplyr::bind_rows(dfrIdxs) |>
              dplyr::mutate(option=factor(option,levels=c(options,"sdmTMB")));
  ggplot(dfrIdxs1,aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci,colour=option,fill=option)) +
    geom_ribbon(alpha=0.3,colour=NA) + geom_point() + geom_line() +
    geom_point(data=dfrACD_MM,mapping=aes(x=YEAR,y=3.4299*totBIOMASS),colour="blue",inherit.aes=FALSE) +
    labs(y="predicted biomass (1,000's t)") +
    wtsPlots::getStdTheme();
  ggplot(dfrIdxs1,aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci,colour=option,fill=option)) +
    geom_vline(xintercept=2020) +
    geom_vline(xintercept=2023:2024,linetype=3) +
    geom_ribbon(alpha=0.3,colour=NA) + geom_point() + geom_line() +
    geom_point(data=dfrACD_MM,mapping=aes(x=YEAR,y=3.4299*totBIOMASS),colour="blue",inherit.aes=FALSE) +
    scale_y_log10() +
    labs(y="predicted biomass (1,000's t)") +
    wtsPlots::getStdTheme();
  ggplot(dfrIdxs1 |> dplyr::filter(year>2000),
         aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci,colour=option,fill=option)) +
    geom_vline(xintercept=2020) +
    geom_vline(xintercept=2023:2024,linetype=3) +
    geom_ribbon(alpha=0.3,colour=NA) + geom_point() + geom_line() +
    geom_point(data=dfrACD_MM |> dplyr::filter(YEAR>2000),
               mapping=aes(x=YEAR,y=3.4299*totBIOMASS),colour="blue",inherit.aes=FALSE) +
    scale_y_log10() +
    labs(y="predicted biomass (1,000's t)") +
    wtsPlots::getStdTheme();

  ##--calculate DHARMA residuals----
  set.seed(12345);
  lstRes = list();
  for (mdl in names(mdls)){
    sim = simulate(mdls[[mdl]],nsim=1000,type="mle-mvn");#--type="mle-mvn" recommended
    lstRes[[mdl]] = dharma_residuals(
                              sim,
                              mdls[[mdl]],
                              plot=TRUE,
                              return_DHARMa=TRUE,
                              test_dispersion=TRUE,
                              test_uniformity=TRUE,
                              test_outliers=TRUE
                             );
  }
  wtsUtilities::saveObj(lstRes,file.path(dirThs,"rda_resids_log_depth-st_rw-dg.RData"));
  for (mdl in names(mdls)){
    plot(lstRes[[mdl]]);
  }

