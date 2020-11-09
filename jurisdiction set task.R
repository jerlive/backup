library(raster)
library(sp)
library(sf)
library(ptinpoly)
library(spatialEco)
library(rgeos)
library(maptools)
library(tidyverse)



"
created:jvithaya@ | 29/10/2020
last modified : 5/11/2020

"

main_table <- data.frame(Area= numeric(0),Package_Count= numeric(0))
getRaw <- function(sc){
    href="C:/Users/jvithaya/Downloads/juris_dependency/bigshapefile.shp"
    a  =shapefile(href)  
    atrops <- read.csv("C:/Users/jvithaya/Desktop/raw_data/atrops.csv")
    vc<-unique(atrops[atrops$Station_Code==sc,]$max_postal_code)
    a=a[a$pincode %in% vc,]
    data <- read.csv("C:/Users/jvithaya/Downloads/bigWeek39.csv")
    data=data[data$station_code %in% sc,]
    val=sc
    result_table <- data.frame(Area_sq_km= numeric(0),Package_Count= numeric(0),Number_of_Subpolygons_disjoint_beyond_500m=numeric(0))
    
    b = aggregate(a, dissolve = TRUE) 
    
    #DT <- data.table(longitude=data$rabbit_long,latitude=data$rabbit_lat)
    #DT_sf = st_as_sf(DT, coords = c("longitude", "latitude"),crs = 4326, agr = "constant")
    
    c <- disaggregate(b) 
    extrac <- buffer(c, width = 0.0045045045045045, dissolve = TRUE)
    extrac <- disaggregate(extrac)
    plot(c)
    disjoint=length(extrac)
    
    f = st_as_sf(c, coords = c("longitude", "latitude"),crs = 4326, agr = "constant")
    #mcat("\n COUNT OF SUBPOLYGONS : ",length(c),"\n\n")
    
    sum_area=round(area(a)/1000000,digits=3)
    #cat(" DETAILS OF SUBPOLYGONS: \n\n")
    for(j in 1:length(c)){
          res_area=as.numeric(round(st_area(f[j,])/1000000,digits=2))
          
      
          query_x=as.vector(data[data$station_code==val,]$rabbit_lat)
          query_y=as.vector(data[data$station_code==val,]$rabbit_long)
      
      
      
      
         extractCoords <- function(sp.df)
         {
            results <- list()
           for(i in 1:length(sp.df@polygons[[1]]@Polygons))
           {
             results[[i]] <- sp.df@polygons[[1]]@Polygons[[i]]@coords
           }
            results <- Reduce(rbind, results)
            results
          }
          vertices<-extractCoords(c[j,])
      
      
          vertices_x <- as.vector(vertices[,1])
          vertices_y <- as.vector(vertices[,2])
      
      
         result=point.in.polygon(query_x, query_y, vertices_y, vertices_x)
         result_table[nrow(result_table)+1, ] <- c(res_area,sum(result==1),disjoint)
  
         
    }
    
    result_table$Area_distribution_percent <- round((result_table$Area_sq_km / sum(result_table$Area_sq_km))*100,digits=2)
    result_table$Package_distribution_percent <- round((result_table$Package_Count / sum(result_table$Package_Count))*100,digits=3)
    result_table$Area_greater_than_eightypercent <- ifelse(result_table$Area_distribution_percent>=80, "YES", "NO")
    result_table$historically_no_packages <- ifelse(result_table$Package_Count==0, "YES", "NO")
    if(j<=3)
    {result_table=cbind(three_or_less_subpolygons = "YES", result_table)}
    else
    {result_table=cbind(three_or_less_subpolygons = "NO", result_table)}
    if(disjoint<=3)
    {result_table=cbind(three_or_less_subpolygons_beyond_500_m = "YES", result_table)}
    else
    {result_table=cbind(three_or_less_subpolygons_beyond_500_m = "NO", result_table)}
    
    result_table=cbind(Number_of_Subpolygons = j, result_table)
    #result_table=cbind(Number_of_Subpolygons_disjoint_beyond_500m = length(extrac), result_table)
    result_table=cbind(Subpolygons=paste(val,rownames(result_table),sep="_"),result_table)
    result_table=cbind(Station_Code = val, result_table)
    result_table
}


sc_list=  c('AMDE','BLRF','DELF','HYDQ','VTZF')
for (entry in sc_list)
{ main_table=rbind(main_table,getRaw(entry))}
cat("\n SUBPOLYGON LEVEL DATA \n\n")
print(main_table)
#write.csv(main_table,"polygonlevel-rawdata.csv")


" EXCEPTION LEVEL WITH TOLERANCE "

s.sf2 <- shapefile("C:/Users/jvithaya/Downloads/juris_dependency/bigshapefile.shp")
data <- read.csv("C:/Users/jvithaya/Downloads/bigWeek39.csv")
stationcode=unique(data$station_code)
atrops <- read.csv("C:/Users/jvithaya/Desktop/raw_data/atrops.csv")

result_table <- data.frame(Station_Code=character(0),Total= numeric(0), Inside_MMI= numeric(0),Outside_MMI= numeric(0), MMI_Exception_percent= numeric(0),Inside_MMI_500=numeric(0),Outside_MMI_500=numeric(0),MMI_500_Exception_percent=numeric(0),NDL_Unresolved_Count=numeric(0),NDL_Percent_Unresolved=numeric(0))
for (val in stationcode){
  vc<-unique(atrops[atrops$Station_Code==val,]$max_postal_code)
  datasc=data[data$station_code==val,]
  
  regexcondition="(._[0-9]{6}$)"
  is_resolvedlist=ifelse(datasc$manifest_route_code=="NO_SORT_CODE",0,ifelse(str_detect(datasc$manifest_route_code,regexcondition),1,0))
  #0 is unresolved and 1 is resolved
  unresolved_count=(sum(is_resolvedlist==0))
  
  scoperesult=ifelse(datasc$scope>27,0,datasc$scope/datasc$scope)
  
  
  t2=s.sf2[s.sf2$pincode %in% vc,]
  
  
  p <- spTransform(t2,CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  p <- gUnaryUnion(p)
  p<-remove.holes(p)
  
  
  b <- gBuffer(t2, width=0.0045045045045045, quadsegs=8,capStyle="ROUND",joinStyle="ROUND", mitreLimit=3)
  b <- gUnaryUnion(b)
  b <-remove.holes(b)
  
  e <- erase(b, p)                         
  plot(e, col="red")    
  
  
  query_x=as.vector(data[data$station_code==val,]$rabbit_lat)
  query_y=as.vector(data[data$station_code==val,]$rabbit_long)
  
  
  
  
  extractCoords <- function(sp.df)
  {
    results <- list()
    for(i in 1:length(sp.df@polygons[[1]]@Polygons))
    {
      results[[i]] <- sp.df@polygons[[1]]@Polygons[[i]]@coords
    }
    results <- Reduce(rbind, results)
    results
  }
  vertices<-extractCoords(p)
  vertices_tolerance<-extractCoords(b)
  
  
  vertices_x <- as.vector(vertices[,1])
  vertices_y <- as.vector(vertices[,2])
  
  
  result=point.in.polygon(query_x, query_y, vertices_y, vertices_x)
  result=ifelse(scoperesult==0,0,result)
  vertices_x<-as.vector(vertices_tolerance[,1])
  vertices_y<-as.vector(vertices_tolerance[,2])
  result_tolerance<-point.in.polygon(query_x, query_y, vertices_y, vertices_x)
  result_tolerance=ifelse(scoperesult==0,0,result_tolerance)
  
  inside<-sum(result==1)
  outside<-sum(result==0)
  out_tol<-sum(result_tolerance==0)
  in_tol=sum(result_tolerance==1)
  
  tot=outside+inside
  result_table[nrow(result_table)+1, ] <- c(val,tot,inside,outside,tot*100,in_tol,out_tol,round((out_tol/in_tol)*100,digits = 3),unresolved_count,round(((unresolved_count)/tot)*100,digits = 3))
}

 
sclevel=data.frame(matrix("", ncol = 0, nrow = length(sc_list)))
sclevel=cbind(sclevel,Station_Code=unique(main_table$Station_Code))
sclevel=cbind(sclevel,Number_of_Subpolygons=aggregate(Number_of_Subpolygons ~ Station_Code, main_table, max)$Number_of_Subpolygons)
sclevel=cbind(sclevel,Number_of_Subpolygons_disjoint_beyond_500m=aggregate(Number_of_Subpolygons_disjoint_beyond_500m ~ Station_Code, main_table, max)$Number_of_Subpolygons_disjoint_beyond_500m)
sclevel=cbind(sclevel,three_or_less_subpolygons=aggregate(three_or_less_subpolygons ~ Station_Code, main_table, max)$three_or_less_subpolygons)
sclevel=cbind(sclevel,three_or_less_subpolygons_beyond_500_m=aggregate(three_or_less_subpolygons_beyond_500_m ~ Station_Code, main_table, max)$three_or_less_subpolygons_beyond_500_m)
sclevel=cbind(sclevel,Majority_area_polygon_more_than_80 = aggregate(Area_greater_than_eightypercent ~ Station_Code, main_table, max)$Area_greater_than_eightypercent)
sclevel=cbind(sclevel,Area_of_biggest_subpolygon_sq_km = aggregate(Area_sq_km ~ Station_Code, main_table, max)$Area_sq_km)
sclevel=cbind(sclevel,historically_no_package_polygon_present=aggregate(historically_no_packages ~ Station_Code, main_table, max)$historically_no_packages)
sclevel=merge(sclevel,result_table,by="Station_Code")

cat("\n STATION LEVEL DATA \n\n")
print(sclevel)
write.csv(main_table,"polygonlevel.csv")
write.csv(sclevel,"stationlevel.csv")

