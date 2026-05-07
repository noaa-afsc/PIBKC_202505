#--run sdmTMB for PIBKC----
##--converting CPUE from t/sq. nm. to t/km^2 prior to fitting models
##--BIG DIFFERENCE FROM 2025-05: using station areas to calculate biomass index
##----from sdmTMB predictions of CPUE in t/km^2 (previously, "area" was 1, the default in sdmTMB::get_index)
require(ggplot2);
require(sdmTMB)

dirThs = dirname(rstudioapi::getSourceEditorContext()$path);
##--get EBS survey data----
if (FALSE){
  ###--survey years to include----
  yrs = 1975:2024;

  ###--get my PIBKC crab survey results (run r_wts.PIBKC.doCalcs.R)----
  lstWTS = wtsUtilities::getObj(file.path(dirThs,"../02_CrabpackDataComparisons/rda_wts.PIBKC.NMFS.RData"));
  #####--Calculate annual survey abundance and biomass by sex/size group for the Pribilof District
  ######--calculate CPUE by haul and population category (encoded as SEX)
  ######--CPUE units are number/sq.nmi for numCPUE and t/sq.nmi. for wgtCPUE
  verbosity=FALSE;
  dfrCPUE.ByXM <-tcsamSurveyData::calcCPUE.ByHaul(lstWTS$dfr.HD,
                                                  lstWTS$dfr.ID,
                                                  bySex=TRUE,
                                                  byMaturity=FALSE,
                                                  byShellCondition=FALSE,
                                                  bySize=FALSE,
                                                  export=FALSE,
                                                  verbosity=verbosity); #--wgtCPUE in t/sq. nmi.
  #####--calculate CPUE by station
  dfrCPUE.ByXM <-tcsamSurveyData::calcCPUE.ByStation(lstWTS$dfr.SD,
                                                     dfrCPUE.ByXM,
                                                     export=FALSE,
                                                     verbosity=verbosity); #--wgtCPUE in t/sq. nmi.
  #####--calculate biomass/abundance by stratum
  #####--units are millions for abundance, 1,000's t for biomass
  dfrACD.ByXM.ByStrata <-tcsamSurveyData::calcAB.ByStratum(tbl_strata=lstWTS$dfr.SD,
                                                           tbl_cpue=dfrCPUE.ByXM,
                                                           export=FALSE,
                                                           verbosity=verbosity);
  #####--calculate biomass/abundance for the Pribilof District
  ######--totAbundance in millions
  ######--totBiomass in 1,000's t
  dfrACD.ByXM<-tcsamSurveyData::calcAB.EBS(dfrACD.ByXM.ByStrata,
                                            export=FALSE,
                                            verbosity=verbosity) |>
                 dplyr::mutate(STRATUM="Pribilof District");
  dfrACD_MM = dfrACD.ByXM |>
                 dplyr::filter(SEX=="4. mature males"); #--keep mature males only
  ggplot(dfrACD_MM,aes(x=YEAR,y=totBIOMASS)) + geom_line() + geom_point();
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
  ####--convert wgtCPUE from t/sq.nmi. to t/km^2
  sqnmi2sqkm = 3.4299; # x/sqnmi / sqnmi2sqkm -> x/sqkm; A[sqnmi]*sqnmi2sqkm -> A[sqkm];
  dfrCPUE_MM = dfrCPUE.ByXM |> dplyr::filter(SEX=="4. mature males") |>                #--mature males only
                 dplyr::select(YEAR,STRATUM,GIS_STATION,numIndivs,numCPUE,wgtCPUE) |>  #--keep only (wgtCPUE in t/sq. nmi)
                 dplyr::inner_join(lstWTS$dfr.SD |> dplyr::select(YEAR,STRATUM,GIS_STATION,
                                                                 STRATUM_AREA,STATION_AREA,STRATUM_AREA_BYSTATION),
                                   by=c("YEAR","STRATUM","GIS_STATION")) |>
                 dplyr::inner_join(lstWTS$evs.stns |> sf::st_drop_geometry() |>
                                     dplyr::select(YEAR,GIS_STATION,LON,LAT,BOTTOM_DEPTH,BOTTOM_TEMP),
                                   by=dplyr::join_by(YEAR,GIS_STATION)) |>
                 sdmTMB::add_utm_columns(ll_names=c("LON","LAT"),
                                         utm_names=c("X","Y"),
                                         units="km") |>
                 dplyr::mutate(numCPUE=numCPUE/sqnmi2sqkm,           #--convert to numbers/km^2
                               wgtCPUE=wgtCPUE/sqnmi2sqkm,           #--convert to t/km^2
                               STRATUM_AREA=sqnmi2sqkm*STRATUM_AREA, #--convert to km^2
                               STATION_AREA=sqnmi2sqkm*STATION_AREA, #--covert to km^2
                               STRATUM_AREA_BYSTATION=sqnmi2sqkm*STRATUM_AREA_BYSTATION #--convert to km^2
                              );
  wtsUtilities::saveObj(dfrCPUE_MM,file.path(dirThs,"rda_CPUE_MM.RData"));
} else {
  dfrACD_MM  = wtsUtilities::getObj(file.path(dirThs,"rda_ACD_MM.RData"));  #--biomass in 1,000's t
  dfrCPUE_MM = wtsUtilities::getObj(file.path(dirThs,"rda_CPUE_MM.RData")); #--wgtCPUE in t/km^2.
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
  #--NOT DOING THIS ONE
  dfr = dfrCPUE_MM |> dplyr::mutate(log_depth=log(BOTTOM_DEPTH),
                                    log_temp = log(BOTTOM_TEMP)) |>
          dplyr::select(year=YEAR,X,Y,wgtCPUE,log_depth,log_temp);    #--wgtCPUE in t/km^2
  ###--model m3----
  mdl_m3 <- sdmTMB::sdmTMB(
              wgtCPUE ~ s(log_depth),  #--in t/km^2.
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
              wgtCPUE ~ s(log_depth),  #--in t/km^2
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
  #--do these (or just mdl_tw_ar)----
  dfr = dfrCPUE_MM |> dplyr::mutate(log_depth=log(BOTTOM_DEPTH),
                                    log_temp = log(BOTTOM_TEMP)) |>
          dplyr::select(year=YEAR,X,Y,wgtCPUE,log_depth,log_temp);  #--in t/km^2.
  ###--model m3_ar----
  mdl_m3_ar = sdmTMB::sdmTMB(
                wgtCPUE ~ s(log_depth),  #--in t/sq.nmi.
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
                wgtCPUE ~ s(log_depth),  #--in t/km^2.
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

##--test sanity----
  for (mdl in names(mdls)){
    cat("###--Sanity for",mdl,"----\n");
    print(sanity(mdls[[mdl]]));
  }
##--sanity(mdl_m3_ar): Non-linear minimizer did not converge: do not trust this model!
##--sanity(mdl_tw_ar): OK!
mdls = mdls["mdl_tw_ar"]; #--proceed with only this one

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
if (FALSE){
  lstPrds = list();
  for (mdl in names(mdls)){
    lstPrds[[mdl]] = predict(mdls[[mdl]],
                             newdata=NULL,
                             type="response",
                             se_fit=FALSE,
                             return_tmb_object=TRUE);
  }
  wtsUtilities::saveObj(lstPrds,file.path(dirThs,"rda_prds_log_depth-st_rw-dg.RData"));
} else {
  lstPrds = wtsUtilities::getObj(file.path(dirThs,"rda_prds_log_depth-st_rw-dg.RData"));
}

  ##--check predictions vs input data
  for (mdl in names(mdls)){
    dfr = lstPrds[[mdl]]$data |> dplyr::select(wgtCPUE,est);
    p = ggplot(dfr,aes(x=wgtCPUE,y=est)) +
          geom_point() + geom_abline(slope=1) +
          labs(x="observed",y="predicted") +
          wtsPlots::getStdTheme();
    print(p)
  }

  ##--calculate indices and save objects
  lstIdxs = list();
  for (mdl in names(mdls)){
    lstIdxs[[mdl]] = get_index(lstPrds[[mdl]],
                               bias_correct=TRUE,
                               area=dfrCPUE_MM$STATION_AREA); #<--BIG CHANGE FROM 2025-05!! (STATION_AREA in km^2)
  }
  wtsUtilities::saveObj(lstIdxs,file.path(dirThs,"rda_idxs_log_depth-st_rw-dg.RData"));

  ##--extract and plot sdmTMB indices----
  lstIdxs1 = list();
  for (mdl in names(mdls)){
    lstIdxs1[[mdl]] = lstIdxs[[mdl]] |>
                        dplyr::select(year,
                                      pred=est,      #--predicted biomass in t now that area was correctly specified (!!)
                                      pred_lci=lwr,
                                      pred_uci=upr,
                                      se) |>         #--for log-scale estimate
                        dplyr::mutate(option=mdl,
                                      pred_cv=se/pred) |>
                        dplyr::select(!se);
  }
  dfrIdxs = dplyr::bind_rows(lstIdxs1);
  p = ggplot(dfrIdxs,aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci,colour=option,fill=option)) +
      geom_ribbon(alpha=0.3) + geom_point() + geom_line() +
      geom_point(data=dfrACD_MM,
                 mapping=aes(x=YEAR,y=totBIOMASS*1000), #--convert from 1,000's t to t now (!!)
                 colour="blue",
                 inherit.aes=FALSE) +
      labs(y="biomass (t)") +
      wtsPlots::getStdTheme();
  print(p);
  print(p+scale_y_log10());

  ##--plot comparisons with REMA models----
  lstREMAs = wtsUtilities::getObj(file.path(dirThs,"../04_SamplingBasedZerosSubstitutions/rda_Zeros_REMA_Models.RData"));
  options = c(unique(lstREMAs$dfrPrds$option)[3:5],unique(dfrIdxs$option));
  dfrIdxs1 = lstREMAs$dfrPrds |>
                dplyr::filter(option %in% options) |>
                dplyr::mutate(obs=obs,               #--original scale is t
                             obs_lci=obs_lci,
                             obs_uci=obs_uci,
                             pred=pred,              #--original scale is t
                             pred_lci=pred_lci,
                             pred_uci=pred_uci
                             ) |>
              dplyr::bind_rows(dfrIdxs) |>
              dplyr::mutate(option=factor(option,levels=c(options,"sdmTMB")));
  p = ggplot(dfrIdxs1,aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci,colour=option,fill=option)) +
        geom_ribbon(alpha=0.3,colour=NA) + geom_point() + geom_line() +
        geom_point(data=dfrACD_MM,
                   mapping=aes(x=YEAR,y=totBIOMASS*1000), #--scale to t from 1,000's t
                   colour="blue",
                   inherit.aes=FALSE) +
        labs(y="predicted biomass (t)") +
        wtsPlots::getStdTheme();
  print(p);
  print(p+scale_y_log10());
  p = ggplot(dfrIdxs1 |> dplyr::filter(year>2000),
            aes(x=year,y=pred,ymin=pred_lci,ymax=pred_uci,colour=option,fill=option)) +
        geom_vline(xintercept=2020) +
        geom_vline(xintercept=2023:2024,linetype=3) +
        geom_ribbon(alpha=0.3,colour=NA) + geom_point() + geom_line() +
        geom_point(data=dfrACD_MM |> dplyr::filter(YEAR>2000),
                   mapping=aes(x=YEAR,y=totBIOMASS),
                   colour="blue",inherit.aes=FALSE) +
        scale_y_log10() +
        labs(y="predicted biomass (1,000's t)") +
        wtsPlots::getStdTheme();
  print(p);

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

