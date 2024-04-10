library(ggplot2)
library(ggtree)
library(ggimage)
library(maps)
library(maptools)
library(scatterpie)
library(tidyverse)

library(rgdal)
library(data.table)


###
# load function
###

# convert "Nine-dash line" shp data
## notes: works in the same way as shp2finaldf()
convertl9 <- function(shp_data,centerlong, centerlat, projection){
  l9 <- fortify(shp_data)
  #shift <- 360 - centerlong
  #l9$long.new <- l9$long + shift
  #l9$long.new <- ifelse(l9$long.new > 180, l9$long.new-360, l9$long.new)
  
  # project coordinates
  PROJ <- paste("+proj=",projection," +lon_0=",as.character(centerlong)," +lat_0=",as.character(centerlat)," +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep='')
  l9[, c("X","Y")] <- data.table(project(cbind(l9$long, l9$lat), proj=PROJ))
  
  return(l9)
}


###########
# load map
###########
# presetting
centerlong = 80
centerlat = 30
projection = 'longlat'

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         panel.border = element_blank(),
                         #legend.background = element_rect(fill=rgb(red = 242, green = 242, blue = 242, max = 255)),#设施图例背景色
                         #legend.key = element_rect(colour = rgb(red = 242, green = 242, blue = 242, max = 255),
                         #fill = rgb(red = 242, green = 242, blue = 242, max = 255)),#设置图例填充色
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         plot.title = element_text(size=10)))

# 1. load basic maps
mymap <-rgdal::readOGR("./map/ESRI世界地图中文版.shp") #读取地图
mymap <- map_data("world")
mymapd <- fortify(mymap)  #打散地图为数据框，方便ggplot读取
p <- ggplot(data = mymapd) +
  # allow to add variant legend and population label
  # guides(fill=FALSE)+ 
  theme_opts+
  geom_polygon(aes(x = long, y = lat, group = group,fill=group), color="grey64", fill = "grey95")
print(p)


# 2. load Nanhai 9 dash line
southChinaSea <- readOGR(dsn = "./resource/中国标准世界国家区划/南海九段线.shp")
southChinaSea <- convertl9(southChinaSea, centerlong, 0, projection)


###########
# load population information & allele frequency
###########
pop_info = read_delim("./reference/KGP.pop.location.tsv", delim = "\t")
allele_freq = read_delim("./input/2_109513601.freq.txt", delim = "\t")
names(allele_freq)
# at least four columns, 
# [1] "pop"     "long"    "lat"     "ALTfreq"


# calculate ref allele freq
allele_freq = allele_freq %>%
  mutate(REFfreq = 1 - ALTfreq)
allele_freq.long = allele_freq %>%
  pivot_longer(cols = c("REFfreq", "ALTfreq"),
               names_to = "FrequencyType",
               values_to = "Value")


###########
# plot
###########
# create basic world map; allowing duplicates
prs = p  
# add population label
## notes: please modify to fit your plot
popList = c("CHB", "CHS", "CDX", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI")
# set pie size
allele_freq$radius = 4  # adjust the pie size
# set pie color
altRefColors <- c("ALTfreq" = "#006d77", "REFfreq" = "#83c5be")  # set pie composition colors

# plot
p_world = prs + 
  # 1) add allele frequency composition: ALF, REF
  geom_scatterpie(data=allele_freq,
                  aes(x=long, y=lat, r=radius),
                  cols=c('ALTfreq','REFfreq'), 
                  color='white', alpha=0.8, size=0.25, sorted_by_radius=F) +
  # 2) set pie color
  scale_fill_manual(values=altRefColors, 
                    name="Variant ID",
                    breaks=c('ALTfreq','REFfreq'), 
                    labels=c('Derived allele frequency','Ancestral allele frequency')) + 
  # 3) add population label
  # geom_label(data = subset(allele_freq, pop %in% popList),
  #            aes(x = long, y = lat, label = pop),
  #            # Adjust this to position your labels
  #            nudge_y = 2,                         
  #            # Adjust text size
  #            size = 3,                            
  #            # Background color of the label
  #            # fill = "lightblue",           
  #            # Adjust padding around text
  #            label.padding = unit(0.2, "lines"),  
  #            # Adjust border size
  #            label.size = 0.25) + 
  geom_text(data = subset(allele_freq, pop %in% popList),
            aes(x = long, y = lat, label = pop),
            nudge_y = 1,  # Adjusts the position of the text above its actual point for clarity
            color = "black",  # Text color
            size = 3,  # Text size, adjust as necessary for visibility
            hjust = 0.5,  # Horizontal adjustment to center text
            vjust = -0.5) +  # Vertical adjustment to move text up or down
  # 4) set world map boundaries
  coord_cartesian(ylim = c(-14,70),xlim = c(-111, 145)) + guides(color=FALSE) + 
  # 5) 
  theme(legend.position = c(0, 1), 
        legend.justification = c(0, 1)) + 
  # 6) ensures that one unit on the x-axis is the same length as one unit on the y-axis
  theme(panel.grid = element_blank(), panel.background = element_blank()) +
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
  # 7) Nine line
  geom_path(data = southChinaSea,aes(x = X, y = Y, group = group), colour = "gray65",size = 0.2)

p_world
ggsave("./WorldMap-allele_freq.pdf", p_world, width = 7.9, height = 3, dpi = 600)
