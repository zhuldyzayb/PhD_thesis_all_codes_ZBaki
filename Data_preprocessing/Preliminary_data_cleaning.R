require(spatstat)
require(spatstat.geom)
require(spatstat.explore)
require(rgdal)
require(lubridate)
require(gdata)
require(sf)
require(sp)
require(MASS)



##################################### GAS FIELD ##############################################

gasfield  = gasfield=readOGR(file.choose())  ## Retrieved in April
unique(gasfield$FIELD_NAME)
groningen = subset(gasfield, gasfield$FIELD_NAME=='Groningen')
g=as(groningen, "SpatialPolygons") 
g=st_as_sf(g, crs=CRS("+proj=utm +zone=31N +units=m")) 
W=as.owin(g)


##################################### WELL LOCATIONS ##############################################
wells = read.csv2(file.choose())  ## well_locations.csv 
wells=wells[, c(1:4)]
wells$Type = rep("", nrow(wells))

production_wells = c("AMR", "BIR","EKL","EKR1", "EKR2", "FRB", "KPD", "LRM", "MWD", "NBR",'NWS', 'OVS', 
          'OWG', 'PAU', 'POS', 'SAP', 'SCB', 'SDB', 'SLO', 'SPI1', 'SPI2', 'SZW1', 'SZW2',
          'TJM', 'UTB', 'ZND', 'TUS', 'ZDP', 'ZVN')

observational_wells = c('BIR', 'BOL', 'BRH', 'DZL', 'FRM', 'HGL', 'HND', 'HRS', 'KHM', 'MDN', 'ODP', 'OLD', 'ROT',
        'SDM', 'SMR', 'SPH', 'SWO', 'TBR', 'UHM', 'UHZ', 'ZBR', 'ZRP', 'ZWD')

injection_wells = c('BRW')

abandoned_wells = c('BRW', 'EKR1', 'HRS', 'PAU', 'SAP', 'SZW2', 'SLO', 'UHM', 'WBL', 'ZND', 'ZVN')
## WBL not present in the coordinates data

## Check for wells statuses
intersect(production_wells, observational_wells)  ## BIR has prod and observation
intersect(production_wells, abandoned_wells) ## EKR1, PAU, SAP, SLO, ZND, ZVN abandoned production wells
intersect(observational_wells, abandoned_wells)  ## HRS UHM abandoned observation wells
intersect(injection_wells, abandoned_wells)  ## BRW abandoned injection well



for(i in 1:nrow(wells)){
  if(length(intersect(production_wells, wells$Well.code[i]))!=0){ wells$Type[i] = 'P'}
  else if(length(intersect(observational_wells, wells$Well.code[i]))!=0){ wells$Type[i] = 'O'}
}
wells$Type[wells$Well.code=='BRW'] = 'I'

for(i in 1:nrow(wells)){
  if(length(intersect(abandoned_wells, wells$Well.code[i]))!=0){ wells$Type[i] = paste(wells$Type[i], "/A", sep='')}
}
wells$Type[wells$Well.code=='BIR'] = 'P/O'
wells$Ind = substr(wells$Type, 1,1)
colnames(wells) = c('Name', 'Code', 'X', 'Y', 'Type', 'Ind')

########  Convert coordinates to UTM system ########
upd <- st_as_sf(wells[,c(3:4)], coords = c("X", "Y"), crs = st_crs(28992))
# Transform to UTM 31
utm31_sf <- st_transform(upd, crs = st_crs(32631))
# Get the transformed coordinates
utm31_coords <- st_coordinates(utm31_sf)

wells$UTM1 = utm31_coords[,1]
wells$UTM2 = utm31_coords[,2]

##################################### PRODUCTION FIGURES ##############################################

gas21 = read.csv2(file.choose()) ## Groningen gas production.csv -> from production figures upto 2021
gas21 = na.omit(gas21)
gas21 = gas21[, -1]
gas21$Month = as.Date(gas21$Month, format = '%d/%m/%Y')
colnames(gas21) = c('Well', 'Date', 'Days', 'Gasprod')

gas23 = read.csv2(file.choose()) ## production figures between 2022-2023
gas23 = na.omit(gas23)
head(gas23)
gas23$Date = as.Date(gas23$Date, format = '%d/%m/%Y')
gas23 = gas23[,-1]
colnames(gas23) = c('Name', 'Well', 'Date', 'Gasprod')

gy = seq(year(min(gas21$Date)), year(max(gas21$Date)),1)
anngas = gy * 0  ## annual production values upto 2021
for(i in 1:length(gy)){
  temp = gas21[year(gas21$Date)==gy[i],]
  anngas[i] = sum(temp$Gasprod)/(10^9) # convert to Nbcm
}

############### Create total monthly and annual production values ################
tm = seq(min(gas21$Date), max(gas23$Date), by='months')
ty = seq(min(gas21$Date), max(gas23$Date), by='year')
gasM = rep(0, length(tm))
gasY = rep(0, length(ty))
## Monthly
for(i in 1:length(tm)){
  if(year(tm[i])<2022){
    temp = gas21[gas21$Date==tm[i],]
    gasM[i] = sum(temp$Gasprod)/(10^9)
  }
  else {
    temp = gas23[gas23$Date==tm[i],]
    gasM[i] = sum(temp$Gasprod)/(10^9)
  }
}
## Annual
for(i in 1:length(ty)){
  if(year(ty[i])<2022){
    temp = gas21[year(gas21$Date)==year(ty[i]),]
    gasY[i] = sum(temp$Gasprod)/(10^9)
  }
  else {
    temp = gas23[year(gas23$Date)==year(ty[i]),]
    gasY[i] = sum(temp$Gasprod)/(10^9)
  }
}
## Cumulative annual
gasCum = rep(0, length(gasY))
for(i in 1:length(gasCum)){
  gasCum[i] = sum(gasY[1:i])
}


##################################### EARTHQUAKES ##############################################
eqs = read.csv(file.choose()) ## all_induces.csv
# convert dates into correct format
eqs$Date = paste(eqs$YYMMDD%%100, "-", (eqs$YYMMDD%%10000 - eqs$YYMMDD%%100)/100, "-", eqs$YYMMDD%/%10000, sep='')
eqs$Date = as.Date(eqs$Date, format = "%d-%m-%Y")
# convert Amersford RD coordinates into UTM-31
cord.dec = SpatialPoints(cbind(eqs$LON,eqs$LAT),proj4string=CRS("+proj=longlat"))
cord.utm = spTransform(cord.dec, CRS("+proj=utm +zone=31 +datum=WGS84 +units=m"))
eqs$UTM1 = cord.utm@coords[,1]
eqs$UTM2 = cord.utm@coords[,2]
## filter earthquakes only within the gas field
ind = inside.owin(eqs$UTM1, eqs$UTM2, W)
eqs = eqs[ind,]
## filter events that happened starting from 1995 (reliable data)
eqs = eqs[year(eqs$Date)>=1995,]
## filter with magnitude threshold M>=1.5
eqs = eqs[eqs$MAG>=1.5,]
## year upto 2021 for model fitting purposes (not compulsory)
eqs21 = eqs[year(eqs$Date)<=2021,]
ann = rep(0, year(max(eqs$Date)) - year(min(eqs$Date))+1)
years = seq(min(eqs$Date), max(eqs$Date), by = 'year')
years = year(years)
## annual counts of earthquakes
ann = years * 0
for(i in 1:length(ann)){
  temp = eqs[year(eqs$Date)==years[i],]
  ann[i] = nrow(temp)
}


##################################### PRESSURE ##############################################
pd = read.csv(file.choose()) ##  Observed pressure data cleaned
pd$DATE = as.Date(pd$DATE, format = '%m/%d/%Y')
yy = seq(year(min(pd$DATE)), year(max(pd$DATE)), 1)
## count available measurements per year
cc = yy*0  
for(i in 1:length(yy)){
  tt = pd[year(pd$DATE)==yy[i],]
  cc[i] = nrow(tt)
}
## Excel format swith from "," to "." for decimal sign
replace_symbol <- function(data, old_symbol, new_symbol) {
  data <- as.data.frame(lapply(data, function(column) {
    if (is.character(column)) {
      gsub(old_symbol, new_symbol, column)
    } else {
      column
    }
  }))
  return(data)
}
## Replace and save a s numeric values
df_replaced <- replace_symbol(pd, ",", ".")
pd = df_replaced
## Convert to a datframe useful for modelling
dates = rep(pd$DATE, ncol(pd)-1)
press = NA
for(i in 2:ncol(pd)){
  if(i==2){press = as.matrix(pd[,i])}
  else { press = rbind(press, as.matrix(pd[,i]))}
}
cols = colnames(pd)[-1]
cols = rep(cols, nrow(pd))

presdf = as.data.frame(cols)
presdf$date = dates 
presdf$press = press
head(presdf)

colnames(presdf) = c('Well', 'Date', 'Press')
ind=(presdf$Press==presdf$Press[1])
presdf = presdf[!ind,]
colnames(presdf) = c('Well', 'Date', 'Press')
presdf$Press = as.numeric(presdf$Press)

## Assign Eemskanaal13 observations into Haarkstede block (another option is to regard them completely) 
for(i in 1:nrow(presdf)){
  if(presdf$Well[i]=='E13'){ presdf$Well[i] = 'HRS'}
}
















