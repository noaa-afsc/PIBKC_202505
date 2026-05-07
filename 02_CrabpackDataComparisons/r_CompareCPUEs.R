#--compare CPUEs from crabpack and tcsamSurveyData for PIBKC
dirThs = dirname(rstudioapi::getSourceEditorContext()$path)

#--get my PIBKC crab survey results (run r_wts.PIBKC.doCalcs.R)----
lstWTS = wtsUtilities::getObj(file.path(dirThs,"rda_wts.PIBKC.NMFS.RData"));
dfrWS = lstWTS$dfrCPUE.ByX |>
          dplyr::select(YEAR,
                        STATION_ID=GIS_STATION,
                        SEX,
                        COUNT=numIndivs,
                        CPUE=numCPUE,
                        CPUE_MT=wgtCPUE)


#--get crabpack specimen data
lstCP = wtsUtilities::getObj(file.path(dirThs,"rda_crabpack_specimen_data.PIBKC.RData"));

#--calculate crabpack cpue
dfrCP = crabpack::calc_cpue(lstCP,
                          species="BKC",
                          region="EBS",
                          district="PRIB",
                          years=1975:2024,
                          sex="all_categories", #  c("male","female"),
                          crab_category="all_categories") |>
          dplyr::filter(CATEGORY!="legal_male") |>
          # dplyr::filter(CATEGORY=="mature_male") |>
          dplyr::rename(SEX=SEX_TEXT) |>
          dplyr::group_by(YEAR,STATION_ID,SEX) |>
          dplyr::summarize(COUNT=sum(COUNT),
                           CPUE=sum(CPUE),            #--in numbers/sq.nmi.
                           CPUE_MT=sum(CPUE_MT)) |>   #--in t/sq.nmi.
          dplyr::ungroup() |>
          dplyr::mutate(SEX=paste0(SEX,"s"));

dfr = dplyr::bind_rows(dfrWS |> dplyr::mutate(TYPE="WS",.before=1),
                       dfrCP |> dplyr::mutate(TYPE="CP",.before=1)) |>
        dplyr::arrange(YEAR,STATION_ID,SEX,TYPE);
##--CPUE and CPUE_MT values agree between crabpack and tcsamSurveyData
###--CPUE in numbers/sq.nmi.
###--CPUE_MT in t/sq.nmi.

#--compare ABs
  ##--extract TCSAM02 biomass/abundance time series for PIBKC
  dfrW_ABs = lstWTS$dfrACD.ByX.PD |>
                dplyr::select(YEAR,
                              SEX,
                              ABUNDANCE=totABUNDANCE,   #--in ??
                              BIOMASS_MT=totBIOMASS) |> #--in ??
                dplyr::filter(SEX %in% c("all males","all females")) |>
                dplyr::group_by(YEAR,SEX) |>
                dplyr::summarize(ABUNDANCE=sum(ABUNDANCE),
                                 BIOMASS_MT=sum(BIOMASS_MT)) |>
                dplyr::ungroup();
  ##--get crabpack biomass/abundance time series----
  dfrC_ABs = crabpack::calc_bioabund(lstCP,
                                    species="BKC",
                                    region="EBS",
                                    district="PRIB",
                                    years=1975:2024,
                                    sex="all_categories", #  c("male","female"),
                                    crab_category="all_categories",
                                    spatial_level="district") |>
          dplyr::filter(CATEGORY!="legal_male") |>
          dplyr::select(YEAR,
                        SEX=SEX_TEXT,
                        ABUNDANCE,           #--in numbers (?)
                        BIOMASS_MT) |>       #--in t (?)
          dplyr::group_by(YEAR,SEX) |>
          dplyr::summarize(ABUNDANCE=sum(ABUNDANCE),
                           BIOMASS_MT=sum(BIOMASS_MT)) |>
          dplyr::ungroup() |>
          dplyr::mutate(SEX=paste0(SEX,"s"));

dfr = dplyr::bind_rows(dfrW_ABs |> dplyr::mutate(TYPE="WS",.before=1),
                       dfrC_ABs |> dplyr::mutate(TYPE="CP",.before=1)) |>
        dplyr::arrange(YEAR,SEX,TYPE);
##--abundance units differ: crabpack in numbers, tcsamSurveyData in millions
###--accounting for scale, values are identical except for known exceptions
##--biomass units differ: crabpack in t, tcsamSurveyData in 1,000's t
###--accounting for scale, values are identical except for known exceptions


