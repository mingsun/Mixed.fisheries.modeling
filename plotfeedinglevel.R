plotfeedinglevel=function (object, species = NULL, time_range, highlight = NULL, 
                           all.sizes = FALSE, include_critical = FALSE, ...) 
{
  if (is(object, "MizerSim")) {
    if (missing(time_range)) {
      time_range <- max(as.numeric(dimnames(object@n)$time))
    }
    params <- object@params
    feed <- getFeedingLevel(object, time_range = time_range, 
                            drop = FALSE)
  }
  else {
    params <- validParams(object)
    feed <- getFeedingLevel(params, drop = FALSE)
  }
  if (length(dim(feed)) == 3) {
    feed <- apply(feed, c(2, 3), mean)
  }
  if (is.null(species)) {
    species <- dimnames(params@initial_n)$sp[!is.na(params@A)]
  }
  sel_sp <- as.character(dimnames(feed)$sp) %in% species
  feed <- feed[sel_sp, , drop = FALSE]
  plot_dat <- data.frame(value = c(feed), Species = factor(dimnames(feed)$sp, 
                                                           levels = dimnames(feed)$sp), w = rep(params@w, each = length(species)))
  if (!all.sizes) {
    for (sp in species) {
      plot_dat$value[plot_dat$Species == sp & (plot_dat$w < 
                                                 params@species_params[sp, "w_min"] | plot_dat$w > 
                                                 params@species_params[sp, "w_inf"])] <- NA
    }
    plot_dat <- plot_dat[complete.cases(plot_dat), ]
  }
  if (include_critical) {
    feed_crit <- getCriticalFeedingLevel(params)[sel_sp, 
                                                 , drop = FALSE]
    plot_dat_crit <- data.frame(value = c(feed_crit), Species = factor(dimnames(feed)$sp, 
                                                                       levels = dimnames(feed)$sp), w = rep(params@w, each = length(species)))
    if (!all.sizes) {
      for (sp in species) {
        plot_dat_crit$value[plot_dat_crit$Species == 
                              sp & (plot_dat_crit$w < params@species_params[sp, 
                                                                            "w_min"] | plot_dat_crit$w > params@species_params[sp, 
                                                                                                                               "w_inf"])] <- NA
      }
      plot_dat_crit <- plot_dat_crit[complete.cases(plot_dat_crit), 
      ]
    }
    p <- ggplot() + geom_line(aes(x = w, y = value, colour = Species, 
                                  linetype = Species, size = Species, alpha = "actual"), 
                              data = plot_dat) + geom_line(aes(x = w, y = value, 
                                                               colour = Species, linetype = Species, alpha = "critical"), 
                                                           data = plot_dat_crit) + scale_discrete_manual("alpha", 
                                                                                                         name = "Feeding Level", values = c(actual = 1, critical = 0.5))
  }
  else {
    p <- ggplot() + geom_line(aes(x = w, y = value, colour = Species, 
                                  linetype = Species, size = Species), data = plot_dat)
  }
  linesize <- rep(0.8, length(params@linetype))
  names(linesize) <- names(params@linetype)
  linesize[highlight] <- 1.6
  p <- p + scale_x_continuous(name = "Size [g]", trans = "log10") + 
    scale_y_continuous(name = "Feeding Level", limits = c(0, 
                                                          1)) + scale_colour_manual(values = params@linecolour) + 
    scale_linetype_manual(values = params@linetype) + scale_size_manual(values = linesize)
  return(p)
}