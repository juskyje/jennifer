#########################################################################
# Joy-Giovanni Matabishi, github: https://github.com/Joy-Giovanni
# Project: calculate NDVI from GEE using R
# study area: South Tyrol
# July 2021
#########################################################################

# first install rgee, set python environment, connect to google account
# check installation link: https://r-spatial.github.io/rgee/#installation
# script based on: https://csaybar.github.io/blog/2020/06/15/rgee_02_io/

#************************************************************#
#         1. NDVI from sentinel 2                            #
#____________________________________________________________#
## 0. INSTALL PACKAGES
install.packages("rgeeExtra")
install.packages("reticulate")
install.packages("gdalcubes")
install.packages("stars")
install.packages("cptcity")
install.packages("tmap")
install.packages("gifski")

## 1. IMPORT PACKAGES
library(rgeeExtra)
library(reticulate)
py_install("pandas") #to be used just the first time
ee_install_set_pyenv(py_env = "YOUR_ENV") #to be used just the first time
rgee::ee_install_upgrade() #to be used just the first time
reticulate::py_install('earthengine-api==0.1.262')  #to be used just the first time
library(gdalcubes)
library(stars)
library(cptcity)
library(tmap)
library(gifski)

# 2. INITIALIZE GEE & GOOGLE DRIVE
ee_install() #to be used just the first time
ee_check() #to be used just the first time to see if you have successfully downloaded all the necessary dependencies
ee_Initialize()
ee_Initialize(drive = TRUE)

##set roi (= region of interest)
sud_ty <- ee$Geometry$Point(c(11.346776, 46.494666))$buffer(1000)

# 3. FOLLOWING SEARCH FUNCTION NOT AVAILALBLE FOR THIS VERSION
#   ee_search_dataset() %>%
#   ee_search_title("sentinel") %>%
#   ee_search_title("MSI") %>%
#   ee_search_dataset_display()                 

# 4. CONNECT TO IMAGE COLLECTION
s2 <- ee$ImageCollection("COPERNICUS/S2_SR")

# 5. FILTER OUT POOR QUALITY PIXELS WITH QUALITY CONTROL BITS BAND

getQABits <- function(image, qa) {
  # Convert decimal (character) to decimal (little endian)
  qa <- sum(2^(which(rev(unlist(strsplit(as.character(qa), "")) == 1))-1))
  # Return a single band image of the extracted QA bits, giving the qa value.
  image$bitwiseAnd(qa)$lt(1)
}

s2_clean <- function(img) {
  # Estimate the NDVI from B8 and B4
  ndvi <- img$normalizedDifference(c("B8", "B4"))
  
  # Extract quality band
  ndvi_qa <- img$select("QA60")
  
  # Select pixels to mask
  quality_mask <- getQABits(ndvi_qa, "110000000000")
  
  # Mask pixels with value zero.
  ndvi$updateMask(quality_mask)$copyProperties(
    img, 
    c('system:id', 'system:time_start','system:time_end')
  )
}


# 6. PREPARE FOR DOWNLOAD
s2_sud_ty <- s2$
  filterBounds(sud_ty)$
  filter(ee$Filter$lte("CLOUDY_PIXEL_PERCENTAGE", 20))$
  filter(ee$Filter$date("2017-01-01", as.character(Sys.Date())))$
  filter(ee$Filter$calendarRange(6, field = "month"))$
  map(s2_clean)

# 7. GET DATES & IDs OF SELECTED IMAGES
nimages <- s2_sud_ty$size()$getInfo()
ic_date <- ee_get_date_ic(s2_sud_ty)

# 8. INTERACTIVE DISPLAY OF RESULTS
Map$centerObject(sud_ty,zoom = 8)
s2_img_list <- list() 
for (index in seq_len(nimages)) {
  py_index <- index - 1
  s2_img <- ee$Image(s2_sud_ty$toList(1, py_index)$get(0))
  s2_img_list[[index]] <- Map$addLayer(
    eeObject = s2_img,
    visParams = list(min = -0.1, max = 0.8, palette = cptcity::cpt("grass_ndvi", 10)),
    name = ic_date$id[index]
  )
}
Reduce('+', s2_img_list)

# 9. DOWNLOAD IMAGES TO PC THROUGH GOOGLE DRIVE
s2_ic_local <- ee_imagecollection_to_local(
  ic = s2_sud_ty,
  scale = 10,
  region = sud_ty,
  via = 'drive'
)

# 10. IMPORT DOWNLOADED IMAGES
rastlist <-
  list.files(
    path = "C:/Users/matjo/OneDrive/Desktop/EURAC/EURAC",
    pattern = '.tif$',
    all.files = TRUE,
    full.names = FALSE
  )

# 11. CREATE RASTER DATA CUBE WITH DIMENSIONS (X,Y,NDVI)
s2_stars <- rastlist %>% 
  read_stars %>% 
  merge %>% 
  st_set_dimensions(names = c("x", "y", "NDVI")) %>% 
  `names<-`("NDVI")

# 12. DEFINE TITLE OF EACH IMAGE
s2_stars %>% 
  st_get_dimension_values(3) %>% 
  substr(
    start = 2,
    stop = 9
  ) %>% 
  as.Date(format="%Y%m%d") %>% 
  as.character() %>% 
  sprintf("South Tyrol, Italy: %s", .) ->  
  s2_new_names

# 13. PUT EVERYTHING IN A TMAP OBJECT
m1 <- tm_shape(s2_stars) +
  tm_raster(
    palette = cpt("grass_ndvi", 20),
    n = 20, 
    title = "NDVI",
    style = "fisher") +
  tmap_style(style = "natural") +
  tm_facets(nrow = 1, ncol = 1) +
  tm_layout(
    frame.lwd = 2,
    panel.label.bg.color = NA,
    attr.outside = TRUE,
    panel.show = FALSE,
    legend.title.size = 1,
    legend.title.fontface = 2,
    legend.text.size = 0.7,
    legend.frame = FALSE,
    legend.outside = TRUE,
    legend.position = c(0.20, 0.15),
    legend.bg.color = "white",
    legend.bg.alpha = 1,
    main.title = sprintf( s2_new_names),
    main.title.size = 1.2,
    main.title.fontface = 2
  )+
  tm_credits(
    text = "Source: Sentinel-2 MSI: MultiSpectral Instrument, Level-2A",
    size = 1,
    just = "right"
  ) 

# grDevices::dev.size("px")
tmap_animation(tm = m1, width = 699*3,height = 555*3,delay = 100)
