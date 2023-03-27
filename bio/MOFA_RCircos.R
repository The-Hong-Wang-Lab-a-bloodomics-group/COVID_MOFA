#library(RCircos)
library(statnet)
library(circlize)
library(MOFA2)
# factor 11 c(1:3,5:9,11:12)
# factor 3 c(1:5,7:11)
data_circos <- data.frame()
for (i in 1:4){
  weights <- get_weights(MOFAobject.trained, 
                         views = names(MOFAobject.trained@data)[i], 
                         factors = 15,
                         as.data.frame = TRUE 
  )
  weights <- weights[order(abs(weights$value),decreasing = T),]
  weights <- weights[1:10,]
  weights$x <- c(-2,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6)
  weights$feature <- gsub(paste0("_",names(MOFAobject.trained@data)[i]),"",weights$feature)
  data_circos <- rbind(data_circos,weights)
  }

# initialize first track -------------------------------------------------
WGCNA::sizeGrWindow(3, 3)
circos.clear()
circos.initialize(data_circos$view, xlim = c(-2.2,1.8))
circos.track(data_circos$view, y = data_circos$value,
             track.height = mm_h(8),
             panel.fun = function(x, y) {
             })
a <- c(names(MOFAobject.trained@data)[1:3],"plasma metabolome")
for(i in 1:4){
  circos.update(sector.index = names(MOFAobject.trained@data)[i], 
                track.index = 1)
  facing <- c("bending.outside",
              "bending.outside",
              "bending.inside",
              "bending.inside")
  circos.text(CELL_META$xcenter, 
              CELL_META$cell.ylim[2] + mm_y(2), # mm_y(5)
              a[i],
              facing = facing[i],
              cex = 0.6)
}
for(i in 1:4){
  weights <- subset(data_circos,data_circos$view == names(MOFAobject.trained@data)[i])
  for(z in 1:10){
    if(weights$value[z] > 0){
      weights$col[z] <- "#ED1E1A"
    }
    if(weights$value[z] < 0){
      weights$col[z] <- "#3F66E1"
    }
  }
  circos.update(sector.index = names(MOFAobject.trained@data)[i],
                track.index = 1)
  circos.barplot(value = weights$value,pos = weights$x,bar_width = 0.3,
                 col = weights$col)
  circos.lines(c(-2.2,2), c(0,0), col = "black")
  circos.yaxis(labels.cex = 0.3,sector.index =names(MOFAobject.trained@data)[i])

}
circos.yaxis(labels.cex = 0.3,sector.index =names(MOFAobject.trained@data)[1])

## initialize second track ---------------------------------------------------------------
circos.track(data_circos$view, y = data_circos$value,
             track.height = mm_h(15),
             panel.fun = function(x, y) {
             })

for(i in 1:4){
  weights <- subset(data_circos,
                    data_circos$view == names(MOFAobject.trained@data)[i])
  weights <- weights[order(abs(weights$value),
                           decreasing = T),]
  weights <- weights[1:10,]
  weights$x <- c(-1.6,-1.4,-1.1,-0.8,-0.4,0,0.4,0.8,1.2,1.5)
  weights$feature <- gsub(paste0("_", names(MOFAobject.trained@data)[i]),
                          "",
                          weights$feature)
  for(z in 1:10){
    if(weights$value[z] > 0){
      weights$col[z] <- "#ED1E1A"
    }
    if(weights$value[z] < 0){
      weights$col[z] <- "#3F66E1"
    }
  }
  circos.update(sector.index = names(MOFAobject.trained@data)[i], 
                track.index = 2)
  circos.text(c(seq(CELL_META$cell.xlim[1]+mm_x(4),
                    CELL_META$cell.xlim[2]-mm_x(4),
                    length=10)),
              c(rep(seq(CELL_META$cell.ylim[2]-mm_y(2),
                        CELL_META$cell.ylim[1]+mm_y(2),
                        length=5),2)), 
              labels = weights$feature,
              col = weights$col,
              font = 2,
              niceFacing = T,
              cex = 0.2+(0.3/(max(abs(data_circos$value))-
                                min(abs(data_circos$value))))*(abs(weights$value)-min(abs(data_circos$value))))
}





