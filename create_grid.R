library(raster)

# env_folder = "R:/IMAS/Antarctic_Seafloor/environmental_data/500mResolution"
# env_list = list("depth"="_bathy500m_shelf_gebco2020_depth.gri", "slope"="_bathy500m_shelf_gebco2020_slope.gri", "npp"="_NPP_SummerAverage_shelf.gri", "currents"="_waom2k_seafloorcurrents_500mInterpolation_shelf.gri", "temp"="_waom2k_seafloortemperature_500mInterpolation_shelf.gri")
# extent = extent(-1.05e6, -1e6, 1e6, 1.05e6)
# fname_out = "C:/Users/cgros/data/20201207_coralnetfull/count/SXY.csv"

create_env_grid <- function(env_folder, env_list, extent, fname_out=NULL) {
  out_df = NULL
  for (env_name in names(env_list)) {
    print(env_name)
    env_file <- paste(c(env_folder, 
                        paste0("Circumpolar_EnvData", env_list[env_name])), 
                      collapse="/")
    layer <- raster(env_file, ext=extent)
    # Convert to dataframe
    layer_df <- as.data.frame(layer, xy=TRUE, na.rm=TRUE)
    # Add col name
    colnames(layer_df) <- c("Longitude", "Latitude", env_name)
    
    if (is.null(out_df)) {
      out_df = layer_df
    }
    else {
      out_df = merge(out_df, layer_df, by=c("Longitude", "Latitude"))
    }
  }
  
  if (!is.null(fname_out)) {
    print(paste0("Saving env grid in: ", fname_out))
    write.csv(out_df, fname_out)
  }
  
  out_df
}