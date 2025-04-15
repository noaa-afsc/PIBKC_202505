#--plot maps of mature male biomass PIBKC from the NMFS EBS survey
require(ggplot2);

#--get data----
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
dfrCPUE_MM = wtsUtilities::getObj(file.path(dirThs,"rda_CPUE_MM.RData"));

#--plot maps----
gisPath = "~/Work/StockAssessments-Crab/Data/GIS/Shapefiles";
#----create EBS basemap
crsWGS84<-wtsGIS::get_crs("WGS84");
bbox = wtsGIS::getBBox(list(bottomleft=list(lon=-175,lat=55),
                            topright  =list(lon=-166,lat=59)))
basemap = wtsGIS::gg_CreateBasemapLayers(bbox=bbox);
#--get survey grid layers
grid = (tcsamSurveyData::gisGetSurveyGridLayers()$grid) |>
         wtsGIS::cropFeaturesToBBox(bbox);

# #----create PI Habitat Conservation Area layers
# hca = wtsGIS::readShapefile(file.path(gisPath,
#                                       "ConservationAreas/pribilof_hca.shp"),
#                             crs=crsWGS84);#--no need to crop features

#--create sf dataset for PIBKC MMB
sfMMB = dfrCPUE_MM |> dplyr::left_join(grid,by=c("GIS_STATION"="STATION_ID")) |>
          sf::st_as_sf();
#--EBS survey grid
p = ggplot() + basemap +
      geom_sf(data=sfMMB |> dplyr::filter((YEAR %in% 1982:1985)|(YEAR %in% 2018:2022),wgtCPUE>0.01),
              mapping=aes(fill=wgtCPUE),colour=NA) +
      geom_sf(data=grid,fill=NA) +
      scale_fill_viridis_c(option="magma") +
      facet_wrap(~YEAR,nrow=2) +
      labs(fill="CPUE \n(1,000's t/sq. nmi)") +
      wtsPlots::getStdTheme() +
      theme(legend.direction="horizontal",
            legend.position="inside",
            plot.title=element_blank(),
            axis.title.x=element_blank(),
            legend.position.inside=c(0.01,0.01),
            legend.justification.inside=c(0.0,0.0));
print(p);
ggplot2::ggsave(file.path(dirThs,"fig_Map-CPUE.pdf"),device="pdf",plot=p,width=9,height=5.0,units="in")
ggplot2::ggsave(file.path(dirThs,"fig_Map-CPUE.png"),device="png",plot=p,width=9,height=5.0,units="in")

