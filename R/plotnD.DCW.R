plotnD.DCW <- function(xvals, yvals, zvals, subclonal.fraction_x, subclonal.fraction_y,
          pngFile, samplename_x, samplename_y, max.plotted.value = NA,
          cluster.locations = NULL, plot_mutations = F)
{
  print("plotnD.DCW running")
  
  colours = colorRampPalette(c("white", "red"))
  if (!is.na(max.plotted.value)) {
    range = list(c(-0.1, max.plotted.value), c(-0.1, max.plotted.value))
    zvals = data.matrix(zvals[xvals <= max.plotted.value,
                              yvals <= max.plotted.value])
    xvals = data.matrix(xvals[xvals <= max.plotted.value])
    yvals = data.matrix(yvals[yvals <= max.plotted.value])
  }
  else {
    #DCW na.rm added 100417
    range = list(c(floor(min(subclonal.fraction_x,na.rm=T) * 10) -
                     1, ceiling(max(subclonal.fraction_x,na.rm=T) * 10) + 1)/10,
                 c(floor(min(subclonal.fraction_y,na.rm=T) * 10) - 1, ceiling(max(subclonal.fraction_y,na.rm=T) *
                                                                        10) + 1)/10)
  }
  plot.data = cbind(subclonal.fraction_x, subclonal.fraction_y)
  #DCW 080417 remove points with zero no.chrs.bearing mut
  plot.data = plot.data[!is.na(rowSums(plot.data)),]
  
  if (!is.na(max.plotted.value)) {
    plot.data = plot.data[plot.data[, 1] <= max.plotted.value &
                            plot.data[, 2] <= max.plotted.value, ]
  }
  if (!is.null(cluster.locations) & !plot_mutations) {
    panel_function = function(...) {
      panel.levelplot(...)
      panel.abline(h = 0:floor(max(plot.data[, 2])))
      panel.abline(v = 0:floor(max(plot.data[, 1])))
      lpoints(cluster.locations, pch = ".", cex = 6, col = "black")
    }
  }
  else if (plot_mutations) {
    panel_function = function(...) {
      panel.levelplot(...)
      panel.abline(h = 0:floor(max(plot.data[, 2])))
      panel.abline(v = 0:floor(max(plot.data[, 1])))
      lpoints(plot.data, pch = ".", cex = 4, col = "black")
    }
  }
  else {
    panel_function = function(...) {
      panel.levelplot(...)
      panel.abline(h = 0:floor(max(plot.data[, 2])))
      panel.abline(v = 0:floor(max(plot.data[, 1])))
    }
  }
  png(filename = pngFile, width = 1500, height = 1000)
  image.wid = 500 * (range[[1]][2] - range[[1]][1])
  image.ht = 500 * (range[[2]][2] - range[[2]][1])
  fig = levelplot(zvals, row.values = xvals, column.values = yvals,
                  xlim = range[[1]], ylim = range[[2]], xlab = list(label = samplename_x,
                                                                    cex = 2), ylab = list(label = samplename_y, cex = 2),
                  scales = list(x = list(cex = 1.5), y = list(cex = 1.5)),
                  col.regions = colours, colorkey = F, panel = panel_function)
  print(fig)
  dev.off()
}
