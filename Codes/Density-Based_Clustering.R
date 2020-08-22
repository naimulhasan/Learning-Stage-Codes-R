# Downloading the file in the Data Scientist Workbench
download.file("https://ibm.box.com/shared/static/th5dsj5txtw052dt2k9i5hekkjhydda2.csv","G:/ML0151/WeatherStations.csv")

# Reading the csv file
WeatherStations <- read.csv("G:/ML0151/WeatherStations.csv", sep =',')
shah<-read.csv("G:/ML0151/Road_Survey.csv", sep =',')
sh<-shah[,c(1:2)]

# What does the data look like?
head(WeatherStations)
head(sh)
# Let's check general information about the data!
str(WeatherStations)

# Stn_Name	Station Name
# Lat	Latitude (North+, degrees)
# Long	Longitude (West - , degrees)
# Prov	Province
# Tm	Mean Temperature (°C)
# DwTm	Days without Valid Mean Temperature
# D	Mean Temperature difference from Normal (1981-2010) (°C)
# Tx	Highest Monthly Maximum Temperature (°C)
# DwTx	Days without Valid Maximum Temperature
# Tn	Lowest Monthly Minimum Temperature (°C)
# DwTn	Days without Valid Minimum Temperature
# S	Snowfall (cm)
# DwS	Days without Valid Snowfall
# S%N	Percent of Normal (1981-2010) Snowfall
# P	Total Precipitation (mm)
# DwP	Days without Valid Precipitation
# P%N	Percent of Normal (1981-2010) Precipitation
# S_G	Snow on the ground at the end of the month (cm)
# Pd	Number of days with Precipitation 1.0 mm or more
# BS	Bright Sunshine (hours)
# DwBS	Days without Valid Bright Sunshine
# BS%	Percent of Normal (1981-2010) Bright Sunshine
# HDD	Degree Days below 18 °C
# CDD	Degree Days above 18 °C
# Stn_No	Climate station identifier (first 3 digits indicate drainage basin, last 4 characters are for sorting alphabetically).


# Creating the main subset: It contains the Station Name, latitude, longitude, highest monthly and lowest monthly temperature
WeatherStations.submain <- subset(WeatherStations, select = c(Stn_Name,Lat,Long,Tx,Tn))
head(WeatherStations.submain)

# Changing column names to cleaner names
colnames(WeatherStations.submain) <- c("Stn_Name","Lat", "Long", "Tmax", "Tmin")
head(WeatherStations.submain)

WeatherStations.submain <- WeatherStations.submain[complete.cases(WeatherStations.submain),]
head(WeatherStations.submain)

library(leaflet)
# Packages used to display the maps in this notebook
library(htmlwidgets)
library(IRdisplay)


# Establishing the limits of our default visualization
lower_lon = -140
upper_lon = -50
lower_lat = 40
upper_lat = 65

lower_lon = 91.847
upper_lon = 91.849
lower_lat = 24.90
upper_lat = 24.92

# Establishing the center of our default visualization
center_lon = (lower_lon + upper_lon)/2
center_lat = (lower_lat + upper_lat)/2
# Setting the default zoom of our default visualization
zoom = 4

subset <- WeatherStations.submain #'subset' stores the subset we will be using for the leaflet map
weather_map <- leaflet(subset) %>% #creating a leaflet map
  setView(center_lon,center_lat, zoom)%>% #setting the default view for our map
  addProviderTiles("OpenStreetMap.Mapnik")%>% #setting the map that we want to use as background
  addCircleMarkers(lng = subset$Long, #the longitude is the longitude of our subset!
                   lat = subset$Lat, #the latitude is the latitude of our subset!
                   popup = subset$Stn_Name, #pop-ups will show the name of station if you click in a data point
                   fillColor = "Black", #colors of the markers will be black
                   fillOpacity = 1, #the shapes will have maximum opacity
                   radius = 4, #radius determine the size of each shape
                   stroke = F) #no stroke will be drawn in each data point
weather_map
saveWidget(weather_map, file="G:/ML0151/weather_map.html", selfcontained = F) #saving the leaflet map in html
display_html(paste("<iframe src=' ", 'weather_map.html', " ' width='100%' height='400'","/>")) #displaying the map !

subset <- sh #'subset' stores the subset we will be using for the leaflet map
road_map <- leaflet(subset) %>% #creating a leaflet map
  setView(center_lon,center_lat, 17)%>% #setting the default view for our map
  addGraticule() %>% 
  addProviderTiles(providers$Esri.WorldImagery)%>%
  addMiniMap(
    tiles = providers$Esri.WorldImagery,
    toggleDisplay = TRUE)%>%#setting the map that we want to use as background
  addCircleMarkers(lng = subset$Longitude, #the longitude is the longitude of our subset!
                   lat = subset$Latitude, #the latitude is the latitude of our subset!
                    #pop-ups will show the name of station if you click in a data point
                   fillColor = "Blue", #colors of the markers will be black
                   fillOpacity = 1, #the shapes will have maximum opacity
                   radius = 4, #radius determine the size of each shape
                   stroke = F) #no stroke will be drawn in each data point
road_map
saveWidget(road_map, file="G:/ML0151/road_map.html", selfcontained = F) #saving the leaflet map in html


# installing the library 'dbscan'
install.packages("dbscan", dependencies = TRUE)
library('dbscan')

# creating a subset containing only the location information
WeatherStations.sub1 <- subset(WeatherStations.submain, select = c(Lat,Long))

# preparing the dataset for dbscan: center and scale.
scaled_WS.sub1 <- scale(WeatherStations.sub1, center = TRUE, scale = TRUE)
head(scaled_WS.sub1)

# assigning clusters for each
clusters_assignments1 <- dbscan(scaled_WS.sub1, eps = 0.138, minPts = 12)
clusters_assignments1

# clusters must be converted to factor before plotting in different colors
clusters_assignments1$cluster <- as.factor(clusters_assignments1$cluster)

# linking the assigned cluster to each station
WeatherStations.sub1$cluster_no <- clusters_assignments1$cluster
head(WeatherStations.sub1)

# function that calculates the centroid of a given cluster_no and dataframe
cluster_centroid <- function(cluster_no, df){
  cluster_df <- df[ which(df$cluster_no==cluster_no), ]
  lat_mean <- mean(cluster_df[["Lat"]])
  long_mean <-  mean(cluster_df[["Long"]])
  return(c(lat_mean,long_mean))
}

# function that calculates all the centroids latitudes of a given dataframe
all_lats_centroids <- function(df){
  all_lats <- numeric ()
  for (cluster in unique (df$cluster_no)){
    all_lats[cluster] <- cluster_centroid (cluster,df) [1]
  }
  return (all_lats)
}
# function that calculates all the centroids longitudes of a given dataframe
all_longs_centroids <- function(df){
  all_longs <- numeric ()
  for (cluster in unique (df$cluster_no)){
    all_longs[cluster] <- cluster_centroid (cluster,df) [2]
  }
  return (all_longs)
}
# storing the latitude of all cluster centroids
lats_centroids1 <- all_lats_centroids (WeatherStations.sub1)
lats_centroids1
# storing the longitude of all cluster centroids
longs_centroids1 <- all_longs_centroids (WeatherStations.sub1)
longs_centroids1

# plotting the graph
subset <- WeatherStations.sub1
color <- colorFactor("Set1", as.factor(subset$cluster_no))
weather_map1 <- leaflet(subset) %>%
  setView(center_lon,center_lat, zoom)%>% 
  addProviderTiles("OpenStreetMap.BlackAndWhite")%>% 
  addCircleMarkers(lng = subset$Long,
                   lat = subset$Lat, 
                   popup = subset$Stn_Name,
                   fillColor = ~color(cluster_no),
                   fillOpacity = 1,
                   radius = 4,
                   stroke = F) %>% 
  addMarkers(lng = longs_centroids1,
             lat = lats_centroids1,
             popup = unique(subset$cluster_no))%>% 
  addLegend("bottomleft",
            pal = color,
            values = ~cluster_no,
            opacity = 1,
            title = "Cluster")

saveWidget(weather_map1, file="G:/ML0151/weather_map1.html", selfcontained = F)

# Preparing the data (Subsetting and Normalizing)

# creating a subset containing only the location
WeatherStations.sub2 <- subset(WeatherStations.submain, select = c(Lat,Long,Tmax,Tmin))

# preparing the dataset for dbscan: center and scale.
scaled_WS.sub2 <- scale(WeatherStations.sub2, center = TRUE, scale = TRUE)
head(scaled_WS.sub2)

# assigning clusters for each
clusters_assignments2 <- dbscan(scaled_WS.sub2, eps = 0.27, minPts = 12)
clusters_assignments2

# clusters must be converted to factor before plotting in different colors
clusters_assignments2$cluster <- as.factor(clusters_assignments2$cluster)
# linking the assigned cluster to each station
WeatherStations.sub2$cluster_no <- clusters_assignments2$cluster
head(WeatherStations.sub2)

# storing the latitude of all cluster centroids
lats_centroids2 <- all_lats_centroids (WeatherStations.sub2)
lats_centroids2

# storing the longitude of all cluster centroids
longs_centroids2 <- all_longs_centroids (WeatherStations.sub2)
longs_centroids2

# plotting the graph
subset <- WeatherStations.sub2
color <- colorFactor("Set1", as.factor(subset$cluster_no))
weather_map2 <- leaflet(subset) %>%
  setView(center_lon,center_lat, zoom)%>% 
  addProviderTiles("OpenStreetMap.BlackAndWhite")%>% 
  addCircleMarkers(lng = subset$Long,
                   lat = subset$Lat, 
                   popup = subset$Stn_Name,
                   fillColor = ~color(cluster_no),
                   fillOpacity = 1,
                   radius = 4,
                   stroke = F) %>%
  addMarkers(lng = longs_centroids2,
             lat = lats_centroids2,
             popup = unique(subset$cluster_no)) %>%
  addLegend("bottomleft",
            pal = color,
            values = ~cluster_no,
            opacity = 1,
            title = "Cluster")

saveWidget(weather_map2, file="G:/ML0151/weather_map2_1.html", selfcontained = F)

install.packages('d3heatmap')
library(d3heatmap)
d3heatmap(mtcars, scale="column", colors="Blues")

install.packages('networkD3')

library(networkD3)
View(MisLinks)
View(MisNodes)
data(MisLinks, MisNodes)
fn<-forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
             Target = "target", Value = "value", NodeID = "name",
             Group = "group", opacity = 0.4)
saveWidget(fn, file="G:/BD0211/network3D.html", selfcontained = F)
?forceNetwork()

