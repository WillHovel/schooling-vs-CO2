#Comparative Kinematics Mega Project
#Import code for midlines data
#By E. Goerig and T. Castro-Santos

#September 20, 2017. Modified many times after...last 01-25-2019.

#The purpose of this code is to compile the files produced by the Matlab curve-fitting output (CurvMapper)into a single dataframe
#Using Tidy structure to facilitate processing and analysis of data.

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)


# Setting up your directory, depending on your computer environment
#setwd("I:/Analysis")
#setwd("C:/Users/valen/Desktop/")
setwd("C:/Users/valen/Desktop/ZEBEGG/swimming/Videos/zebegg analyses")

rm(list=ls()) # kill everything first

#Multisp<-read_excel("I:/3. multispecies spreadsheet/Multi-species spreadsheet.xlsx")
Multisp<-read_excel("C:/Users/valen/Desktop/ZEBEGG/swimming/Videos/zebegg analyses/mega zebegg spreadsheet.xlsx")

Multisp$clip_id <- as.integer(Multisp$clip_id)

Multisp$fn<-paste(Multisp$trimmed_filename,".mp4_CURVES.xls",sep="")

#First we need a blank dataframe to hold all the data from the Excel files
frames<-data.frame(matrix(ncol=14))
names(frames)<-c("clip_id",1:13)
df<-data.frame(matrix(ncol=28))
names(df)<-c("clip_id","location","x01","y01","x02","y02","x03","y03","x04","y04","x05","y05","x06","y06","x07","y07","x08","y08","x09","y09","x10","y10","x11","y11","x12","y12","x13","y13")

# accomodate files with <> 10 midlines
pad2<-data.frame(matrix(ncol=2))
pad4<-data.frame(matrix(ncol=4))
pad6<-data.frame(matrix(ncol=6))

# create a variable for location along the body (1-200)
location<-c(1:200)

setwd("C:/Users/valen/Desktop/ZEBEGG/swimming/Videos/zebegg curvemapper")
#setwd("J:/Mega kinematic all videos/2. curvemapper midlines xls")

for (i in as.vector(Multisp$clip_id))

  {
  #i<-9
  fn<-Multisp$fn[Multisp$clip_id==i]
  file<-read_excel(fn,col_names=FALSE)
  nc<-ncol(file)
  ifelse(nc==20,file<-cbind(file,pad6),file<-file)
  ifelse(nc==22,file<-cbind(file,pad4),file<-file)
  ifelse(nc==24,file<-cbind(file,pad2),file<-file)
  
  framesi<-file[1,c(1,3,5,7,9,11,13,15,17,19,21,23,25)]
  framesi$clip_id<-i
  
  names(framesi)<-c(1:13,"clip_id")
  frames<-rbind(frames,framesi)
  
  datai<-file[2:201,]
  datai$clip_id<-i
  datai$location<-location
  names(datai)<-c("x01","y01","x02","y02","x03","y03","x04","y04","x05","y05","x06","y06","x07","y07","x08","y08","x09","y09","x10","y10","x11","y11","x12","y12","x13","y13","clip_id","location")
  df<-rbind(df,datai)
}

# Arrange all the x data into a single column
xdata1<-select(df,clip_id,location,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13)
xdata<-gather(xdata1,key=midline,value=x,x01,x02,x03,x04,x05,x06,x07,x08,x09,x10,x11,x12,x13)
xdata<-arrange(xdata,clip_id,midline,location)
xdata$midline<-substr(xdata$midline,2,3)
xdata<-select(xdata,clip_id,midline,location,x)

# Arrange all the y data into a single column
ydata1<-select(df,clip_id,location,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12, y13)
ydata<-gather(ydata1,key=midline,value=y,y01,y02,y03,y04,y05,y06,y07,y08,y09,y10,y11,y12, y13)
ydata<-arrange(ydata,midline,clip_id)
ydata$midline<-substr(ydata$midline,2,3)
ydata<-select(ydata,clip_id,midline,location,y)

midlines<-inner_join(xdata,ydata)

#Delete midlines where x and y = NA. Those were introduced above to accomodate for variability in the number of digitized midlines.
row.has.na <- apply(midlines, 1, function(x){any(is.na(x))})
sum(row.has.na)
midlines<- midlines[!row.has.na,]

midlines<-group_by(midlines,clip_id)
midlines$midline <- as.numeric(midlines$midline)
nb_midlines<-summarize(midlines,nb=max(midline))
unique(midlines$clip_id)
hist(nb_midlines$nb)

frames<-rename(frames,'01'='1','02'='2','03'='3','04'='4','05'='5','06'='6','07'='7','08'='8','09'='9')

frames<-gather(frames,key=midline,value=frame,'01','02','03','04','05','06','07','08','09','10','11','12','13')
frames<-arrange(frames,clip_id,midline)

frames<-na.omit(frames)

frames$midline<-as.numeric(frames$midline)

midlines<-inner_join(midlines,frames)

midlines<-select(midlines,clip_id,midline,location,frame,x,y)

# Save data
#write.csv(midlines, file = "I:/Analysis/Midlines.csv")
write.csv(midlines, file = "C:/Users/valen/Desktop/ZEBEGG/swimming/Videos/zebegg analyses/MidlinesZebEgg.csv")




#################################################################################

# End Data compilation
# from here we will begin with translation and rotation exercise with the code 'Transformer'

#################################################################################











