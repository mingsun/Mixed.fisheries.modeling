plotbiomass=function (sim, species=dimnames(sim@params@initial_n)$sp[!is.na(sim@params@A)], time_range , y_ticks = 6, ylim = c(NA, 
                                                                                NA), total = FALSE, background = TRUE, highlight = NULL, 
                      ...) 
{
  if (missing(time_range)) {
    time_range <- c(as.numeric(dimnames(sim@n)[[1]][1]),as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]))
  }
  start_time<-range(as.numeric(time_range))[1]
  end_time<-range(as.numeric(time_range))[2]
b <- getBiomass(sim, ...)
 
  b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) & 
           (as.numeric(dimnames(b)[[1]]) <= end_time), , drop = FALSE]
  b_total = rowSums(b)
  if (total) {
    b <- cbind(b, Total = b_total)
    species <- c("Total", species)
  }
  names(dimnames(b)) <- c("time", "Species")
  bm <- reshape2::melt(b)
  bm$Species <- as.factor(bm$Species)
  min_value <- 1e-20
  bm <- bm[bm$value >= min_value & (is.na(ylim[1]) | bm$value >= 
                                      ylim[1]) & (is.na(ylim[2]) | bm$value <= ylim[1]), 
  ]
  spec_bm <- bm[bm$Species %in% species, ]
  x_label <- "Year"
  y_label <- "Biomass [g]"
  
  p <- ggplot(spec_bm, aes(x = time, y = value)) + scale_y_continuous(trans = "log10", 
                                                                        breaks = log_breaks(n = y_ticks), labels = prettyNum, 
                                                                        name = y_label) + scale_x_continuous(name = x_label) + 
    scale_colour_manual(values = sim@params@linecolour) + 
    scale_linetype_manual(values = sim@params@linetype)
  if (background) {
    back_sp <- dimnames(sim@n)$sp[is.na(sim@params@A)]
    back_bm <- bm[bm$Species %in% back_sp, ]
    if (nrow(back_bm) > 0) {
      p <- p + geom_line(aes(group = Species), data = back_bm, 
                         colour = sim@params@linecolour["Background"], 
                         linetype = sim@params@linetype["Background"])
    }
  }
  linesize <- rep(0.8, length(sim@params@linetype))
  names(linesize) <- names(sim@params@linetype)
  linesize[highlight] <- 1.6
  p <- p + scale_size_manual(values = linesize) + geom_line(aes(colour = Species, 
                                                                linetype = Species, size = Species))
  return(p)
}
