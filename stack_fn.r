# clean working environment
rm(list=ls())

# set the wprking directory. provide path of GlacierSMBM_ICIMOD_2020_day_4_5
setwd('/home/nirav/Desktop/ICIMOD_R_Training_2020/GlacierSMBM_ICIMOD_2020_day_4_5/')

# define output path/relative path as the output folder is inside working directory
outpath='./processed_input_utm45N/'
newout = './stack/test/'
# load required library
library(raster)

# load srtm 90m dem
srtm_dem=raster("./srtm_dem/New Dem/SRTM_Kali.tif")

# visualize srtm dem
plot(srtm_dem)

# load polygon boundary of RS glacier
RS=shapefile("./Rikha_Samba_polygon/Kali_Gandaki_Glacier.shp")

#visualize RS polygon, overlay on srtm dem
plot(RS,add=T)

# crop srtm dem and overlay RS glacier boundry
plot(crop(srtm_dem,RS))
plot(RS,add=T)

# increase extent of RS
extent(RS)
inc_extent=extent(83.17, 84.35 , 28.46 ,29.34)

# crop srtm dem using new RS glacier extent and make plot
srtm.RS=crop(srtm_dem,inc_extent)
plot(srtm.RS)
plot(RS,add=T)

# change projection of RS glacier boundary
utm44.RS=spTransform(RS,CRS('+init=epsg:32644'))
projection(utm44.RS)

# change projection of srtm RS and resample to 30m. We want to run model with 30m resolution
utm44.30m.srtm.RS=projectRaster(srtm.RS,res=30,method="bilinear",crs="+init=epsg:32644")
projection(utm44.30m.srtm.RS)

# mask RS srtm 30m using RS boundary
srtm.30m.RS.masked=mask(utm44.30m.srtm.RS,utm44.RS)
plot(srtm.30m.RS.masked)
plot(utm44.RS,add=T)

# export 30m srtm elevation of RS glacier
writeRaster(srtm.30m.RS.masked,filename = paste0(outpath,'RS_srtm_elv_30m.tif'),format="GTiff",overwrite=T)

# create standard glacier mask from this dem
glacier_mask.30m=srtm.30m.RS.masked
glacier_mask.30m[glacier_mask.30m>0]<-1 # replace values greater than 1 0 by 1 # format require for SMBM
plot(glacier_mask.30m)

# export glacier mask
writeRaster(glacier_mask.30m,paste0(outpath,'RS_glacier_mask.30m_2010.tif'),format="GTiff",overwrite=T)

#import era5 land topography
orography=raster("./era5_land_orography/geo_1279l4_0.1x0.1.grib2_v4_unpack.nc")

# get detail of era5 and orography
orography
library(ncdf4)
nc_open("./era5_land_orography/geo_1279l4_0.1x0.1.grib2_v4_unpack.nc")

# make a plot and overlay RS glacier
plot(orography)
plot(RS,add=T)

# crop era5 orography and vusualize
plot(crop(orography,extent(RS)))
plot(RS,add=T)

#crop using extended extent, we devided by 9.8to convert geopotential height to elevation
era5.land.oro.RS=crop(orography,inc_extent)/9.8
plot(era5.land.oro.RS)
plot(RS,add=T)

# reproject era5 orography to utm45n and 30 m using to make this consistent with srtm 30m
era5.land.30m=disaggregate(era5.land.oro.RS,10)

# reproject era5 oRS elevation to match extent
era5.land.30m=projectRaster(era5.land.30m,glacier_mask.30m)

# mask era5 land resampled elevation using glacier boundary
era5.land.30m.masked=mask(era5.land.30m,glacier_mask.30m)
plot(era5.land.30m.masked);plot(utm44.RS,add=T)

# export RS glacier era5 30m elevation
writeRaster(era5.land.30m.masked,filename = paste0(outpath,'RS_era5_land_elv_30m.tif'),format="GTiff",overwrite=T)

# create MODIS path
modis_path='./improved_daily_modis_snow_2010/'

# create list of all MODIS file in 2010
modis_list=list.files(modis_path,pattern="2010",full.names=T)

# print files in modis list
modis_list

# make a loop and create snow and ice mask
daily_snow_stack=stack()
daily_ice_stack=stack()

is_first <- TRUE
for  (i in 1:length(modis_list)){
  i_modis=raster(modis_list[i])
  modis_resample=projectRaster(i_modis,glacier_mask.30m,method="ngb")
  modis_mask=mask(modis_resample,glacier_mask.30m)
  
  snow_mask=modis_mask
  snow_mask[snow_mask==200]<-1
  snow_mask[snow_mask==242]<-1
  snow_mask[snow_mask==252]<-1
  snow_mask[snow_mask!=1]<-0
  
  # ice_mask=modis_mask
  # ice_mask[ice_mask==250]<-1
  # ice_mask[ice_mask==238]<-1
  # ice_mask[ice_mask==239]<-1
  # ice_mask[ice_mask!=1]<-0
  
  if (is_first == TRUE){
    writeRaster(snow_mask,paste0(outpath,'RS_daily_snow_mask_2010_30m.tif'),format='GTiff',overwrite=T)
  #  writeRaster(ice_mask,paste0(outpath,'RS_daily_ice_mask_2010_30m.tif'),format='GTiff',overwrite=T)
    print(i)
    is_first  <- FALSE
    next 
    }
  temp_snow_stack = stack()
#  temp_ice_stack = stack()
  # old_snow_stacked_tiff = )
  # old_ice_stacked_tiff = raster()
  
 
  temp_snow_stack = stack(temp_snow_stack,stack("./processed_input_utm45N/RS_daily_snow_mask_2010_30m.tif"))
# print("Type of Temp_Snow_Stack", typeof(temp_snow_stack))
#  temp_ice_stack = stack(temp_ice_stack,stack("./processed_input_utm45N/RS_daily_ice_mask_2010_30m.tif"))

  temp_snow_stack = stack(temp_snow_stack, snow_mask)
#  temp_ice_stack = stack(temp_ice_stack, ice_mask)

  # export output modis mask stack
  writeRaster(temp_snow_stack, paste0(outpath,'RS_daily_snow_mask_2010_30m.tif'),format='GTiff',overwrite=T)
 # print("Number of Layers in Snow_Mask", nlayers(raster(./stack/test/snow_mask.tif)))
#  writeRaster(temp_ice_stack, paste0(outpath,'RS_daily_ice_mask_2010_30m.tif'),format='GTiff',overwrite=T)
  
  print(i)
  gc()

}
  
  
# for (i in 1:length(modis_list)){
#   #i=1
#   gc()
  
#   i_modis=raster(modis_list[i])
#   modis_resample=projectRaster(i_modis,glacier_mask.30m,method="ngb")
#   modis_mask=mask(modis_resample,glacier_mask.30m)
  
#   snow_mask=modis_mask
#   snow_mask[snow_mask==200]<-1
#   snow_mask[snow_mask==242]<-1
#   snow_mask[snow_mask==252]<-1
#   snow_mask[snow_mask!=1]<-0
#   daily_snow_stack=stack(daily_snow_stack,snow_mask)
  
#   ice_mask=modis_mask
#   ice_mask[ice_mask==250]<-1
#   ice_mask[ice_mask==238]<-1
#   ice_mask[ice_mask==239]<-1
#   ice_mask[ice_mask!=1]<-0
#   daily_ice_stack=stack(daily_ice_stack,ice_mask)
  
#   print(i)
# }

# export output modis mask stack
# writeRaster(daily_ice_stack,paste0(outpath,'RS_daily_ice_mask_2010_30m.tif'),format='GTiff',overwrite=T)
# writeRaster(daily_snow_stack,paste0(outpath,'RS_daily_snow_mask_2010_30m.tif'),format='GTiff',overwrite=T)

# prepare debris mask/ rikhs samba being debris free lacier
RS_debris_mask=glacier_mask.30m
RS_debris_mask[RS_debris_mask==1]<-0

# mask using glacier mask
RS_debris_mask=mask(RS_debris_mask,glacier_mask.30m)

# export debris mask
writeRaster(RS_debris_mask,paste0(outpath,'RS_debris_mask_30m.tif'),format='GTiff',overwrite=T)

# prepare debris thickness map
RS_debris_thickness=RS_debris_mask

# mask using glacier mask
RS_debris_thickness=mask(RS_debris_thickness,glacier_mask.30m)

# export debris mask
writeRaster(RS_debris_thickness,paste0(outpath,'RS_debris_thickness_mask_30m.tif'),format='GTiff',overwrite=T)


