#Maps

install.packages("leaflet")
library(leaflet)
map<-leaflet()
map
map<-leaflet() %>% addTiles()
map<-leaflet() %>% addTiles() %>% addMarkers(lng=-73.9851,lat=40.7589)
map<-leaflet() %>% addTiles() %>% 
  addMarkers(lng=-73.9851,lat=40.7589,popup='Times Square')
map<-leaflet() %>% addProviderTiles("Stamen.Watercolor") %>%
        addMarkers(lng=2.2945,lat=48.8584,
                   popup="Eiffel Tower")
#Stamen.TonerHybrid
#https://leaflet-extras.github.io/leaflet-providers/preview/
head(quakes)
map<-leaflet(quakes) %>% addTiles() %>%
      addCircleMarkers(lng=quakes$long,lat=quakes$lat)
map<-leaflet(quakes) %>% addTiles() %>%
  addMarkers(lng=quakes$long,lat=quakes$lat)
map<-leaflet(quakes) %>% addTiles() %>%
  addCircles(lng=quakes$long,lat=quakes$lat)
install.packages("htmlwidgets")
library(htmlwidgets)
install.packages("IRdisplay")
library(IRdisplay)

map<-leaflet(quakes) %>% addTiles() %>%
  addMarkers(clusterOptions=markerClusterOptions())
saveWidget(map,file="map7.html",selfcontained=F)
display_html(paste("<iframe src=' ",'map7.html',"'width='100%' height='300'","/>"))

map<-leaflet(quakes) %>% addTiles() %>%
  addMarkers(lng=86.92,lat=27.99,
  popup ="Mount Everest") %>%
  addRectangles(86.9,27.95,87,28.05)
