# Code for rotating and translating midlines in prep for kinematics analysis
#Dependent on Compiler, which brings data into a single tidy matrix
#This code adjusts for rotation of camera angle and centers the average position of the nose at (0,0)
# T. Castro-Santos & E. Goerig
# September 23, 2017
####################################################


library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)


# Set directory 
#setwd("J:/Mega kinematic all videos/Analysis")
setwd("C:/Users/valen/Desktop/ZEBEGG/swimming/Videos/zebegg analyses")

rm(list=ls())

midlines<-read.csv("MidlinesZebEgg.csv")
midlines<-select(midlines,2:7)#Check against raw file to make sure structure is the same

#Multisp<-read_excel("J:/Mega kinematic all videos/3. multispecies spreadsheet/Multi-species spreadsheet.xlsx")

Multisp<-read_excel("C:/Users/valen/Desktop/ZEBEGG/swimming/Videos/zebegg analyses/mega zebegg spreadsheet.xlsx")
Multisp$clip_id <- as.integer(Multisp$clip_id)

######################################################################################################################################

#### STANDARDIZATION ####
#First step is to standardize everything to body lengths
#Adjust length values. Those with no known TL, or digitized with a TL of 1 = 1. These are stored in a 2nd variable for TL: TL2
# These adjusted values will be used mainly for this initial standardization.Otherwise, it is hard to plot the midlines as they are on a different scale than the other ones.
TL<-select(Multisp,clip_id,TL,species)
TL$TL<-as.numeric(TL$TL)
TL$TL<-ifelse(is.na(TL$TL),1,TL$TL)
TL$TL2 <- ifelse(TL$clip_id %in% c(8:11,16:18,30,49,62,101:121,144:146,155:157,160:187),1,TL$TL)
#Add TL values to the Midlines dataset.
midlines<-inner_join(midlines,TL,by="clip_id")
midlines$TL<-as.numeric(midlines$TL2)

#Adjust midlines scale for clip 11, which appears to have been off by a factor of 100 (typo when digitizing?)
midlines$x<-ifelse(midlines$clip_id==11,midlines$x/100,midlines$x)
midlines$y<-ifelse(midlines$clip_id==11,midlines$y/100,midlines$y)

#Standardize to BL using TL2
midlines$xBL<-midlines$x/midlines$TL2#xBL and yBL are x and y, converted to body lengths
midlines$yBL<-midlines$y/midlines$TL2

###############################################################################################################################
#### ROTATION ####


#TCS effort
#I have now convinced myself that we don't want to translate the x and y coordinates to the origin until after we have rotated the image.  The good news is that we have the 'location' variable, so we can still select just the parts of the fish that are least vulnerable to extreme values/outliers
#For now I will try doing this twice, once with the full length of the fish, and once based on 20:150 (10-75% of body length).  By leaving the axis of motion in place (i.e. no translation) we should actually improve the accuracy of the centering.  At least that is what I am telling myself.


#Clip ID 116 (Tuna) is reversed, so let's re-reverse it
midlines$xBL<-with(midlines,ifelse(clip_id==116,-xBL+0.5,xBL))


midlines$clip_id<-as.numeric(midlines$clip_id)

chop<-midlines[(midlines$location>=10&midlines$location<=190),]#This is to allow us to select a subset of the total midline for fitting regressions...i.e. we might want to fit only to the lead 90% of the fish's body length
clipid<-unique(chop$clip_id)

models<-data.frame(matrix(ncol=4))

names(models)=c("clip_id","a","b","R2")#parameters for intercept (a),slope (b), and r-squared

for (i in as.vector(clipid))
{
  a<-subset(chop,clip_id==i)
  lmi<-lm(yBL~xBL, data=a)
  sum.lmi<-summary(lmi)
  a <- sum.lmi$coef[1,1]   #this gives you the intercept
  b <- sum.lmi$coef[2,1]    #this gives you the the slope
  R2 <- sum.lmi$r.squared  #this gives you the R2
  modeli<-c(i,a,b,R2)
  models<-rbind(models,modeli)
  
}

#All slopes suggest very little rotation actually, but let's take a look at the data and make sure things look right.

# for (i in as.vector(clipid))
# {
#  plot<-paste("Plots/",i,".png", sep="")#this tells it to place the plots in the "Plots" subdirectory (in the working directory)
# png(plot,width=7.5, height=5, units="in",res=300)#we'll be saving plots as .png files
# title=paste("Clip ID - ",i,sep="")#title for each plot
# a<-subset(chop,clip_id==i)#we will do just one plot ta a time
# p<-''#clear the ggplot object
# p<-ggplot(data=a,aes(x=xBL,y=yBL)) + 
#  geom_path(aes(group=midline, color=midline))+
#  xlim(-0.5,1.5)+ylim(-1, 1)+
#  ggtitle(title)+
#  theme (panel.background = 
#             element_rect(fill = "transparent", colour = NA), 
#           plot.background = 
#     element_rect(fill = "transparent", colour = NA),
#             axis.line = element_line(color = 'black'))
#    #set up plots with clear backgrounds
#    #note we use geom_path instead of geom_line. This is because we have >y per x in some cases and geom-path forces it to maintain the order in which the     data were entered.
#   print(p)
#   dev.off()
#  }
# 

#TCS: OK it is clear from these plots that we have some issues still.
#Some, like Clip 67 look like they have some kind of systematic error...I don't really believe the fish were swimming like this...body always bent in one direction...

models<-na.omit(models)#get rid of na's
models$alpha<-atan(models$b)#calculate the angle

models$theta<-2*pi-models$alpha#This adjusts theta to ensure we rotate in the proper direction...negative theta for clockwise rotation


#The linear formulas for adjusting x and y are 
#x'= x cos(theta) - y sin(theta)
#y' = x sin(theta) + y cos(theta)
#note though that our origin varies for each regression model.
#The location of the origin is (0, a), where a is the intercept term of the model
#This means that in order to rotate about this, we need to adjust the y coordinate to account for the y-intercept
# thus the formulas above must be modified:
# x'= x*cos(theta) - (y-a)*sin(theta)
# y' = x*sin(theta) + (y-a)*cos(theta) + a

#again, we could probably do this efficiently with a matrix operation,
#but a for-loop should work nicely

for (i in as.vector(clipid))
{
  a<-models$a[models$clip_id==i]
  theta<-models$theta[models$clip_id==i]
  midlines$xrot[midlines$clip_id==i]<-
    midlines$xBL[midlines$clip_id==i]*cos(theta)-
    (midlines$yBL[midlines$clip_id==i]-a)*sin(theta)
  midlines$yrot[midlines$clip_id==i]<-
    midlines$xBL[midlines$clip_id==i]*sin(theta) + 
    (midlines$yBL[midlines$clip_id==i]-a)*cos(theta)+a
}


 for (i in as.vector(clipid))
 {
   plot<-paste("PlotsRot/",i,".png", sep="")#this tells it to place the plots in the "Plots" subdirectory (in the working directory)
   png(plot,width=7.5, height=5, units="in",res=300)
   title=paste("Clip ID - ",i,sep="")
   a<-subset(midlines,clip_id==i)
   p<-''
  p<-ggplot(data=a,aes(x=xrot,y=yrot)) + 
     geom_path(aes(group=midline, color=midline))+
     xlim(-0.5,1.5)+ylim(-1, 1)+
     ggtitle(title)+
     theme (panel.background = element_rect(fill = "transparent", colour = NA), 
            plot.background = element_rect(fill = "transparent", colour = NA),
            axis.line = element_line(color = 'black'))
   print(p)
   dev.off()
 }



###################################################################################################
#### TRANSLATION ####

#  Translate images so the nose centers on 0,0


# note this might not work well for images that failed to rotate neatly (Amphioxus)
# Also note that fish with asymmetrical strides will not rotate properly
# In those cases we might need to either select a different stride or else consider digitizing more than one stride per fish 

midlines<-group_by(midlines,clip_id,midline)
minx<-(summarize(midlines,minx=min(xrot)))#I am sticking minx in a temporary dataframe because we don't really need to retain this value
minx<-ungroup(minx)
midlines<-ungroup(midlines)
midlines<-inner_join(midlines,minx,by=c("clip_id","midline"))#looks like we might have to join them
midlines<-ungroup(midlines)
midlines$X<-midlines$xrot-midlines$minx #The final X, shifted so the nose is at (0,y)
midlines$Y<-midlines$yrot#y has already been adjusted with the rotation procedure above.  Note that we do NOT set all noses at 0,0, but allow oscillation about y.  This is important for the analyses to come
rm(minx)


for (i in as.vector(clipid))
{
  plot<-paste("PlotsAdj/",i,".png",sep="")#this tells it to place the plots in the "Plots" subdirectory (in the working directory)
 # plot1<-paste("PlotsAdj/",i,".emf",sep="")
 # plot2<-paste("PlotsAdj/",i,".pdf",sep="")
  png (plot,width=7.5, height=5, units="in",res=300)# save in .png
  title=paste("Clip ID - ",i,sep="")
  a<-subset(midlines,clip_id==i)
  p<-''
  p<-ggplot(data=a,aes(x=X,y=Y)) +
    geom_path(aes(group=midline, color=midline))+
    xlim(-0.5,1.5)+ylim(-1, 1)+
    ggtitle(title)+
    theme (panel.background = element_rect(fill = "transparent", colour = NA),
           plot.background = element_rect(fill = "transparent", colour = NA),
           axis.line = element_line(color = 'black'))
  print(p)
 # ggsave(plot1,width=7.5, height=5,dpi=600)# save in .emf
 # ggsave(plot2,width=7.5, height=5,dpi=600)# save in .pdf
  dev.off()
}

########################################################################

# save data with a different name 'just in case'
write.csv(midlines, file = "C:/Users/valen/Desktop/ZEBEGG/swimming/Videos/zebegg analyses/MidlinesZebEggAdj.csv")


# END

