plotm2=function (object, species = NULL, time_range, highlight = NULL, 
                 ...) 
{
  if (is(object, "MizerSim")) {
    if (missing(time_range)) {
      time_range <- max(as.numeric(dimnames(object@n)$time))
    }
    params <- object@params
  }
  else {
    params <- validParams(object)
  }
  pred_mort <- getPredMort(object, time_range = time_range, 
                           drop = FALSE)
  if (length(dim(pred_mort)) == 3) {
    pred_mort <- apply(pred_mort, c(2, 3), mean)
  }
  if (is.null(species)) {
    species <- dimnames(params@initial_n)$sp[!is.na(params@A)]
  }
  species_levels <- c(as.character(params@species_params$species), 
                      "Background", "Resource", "Total")
  pred_mort <- pred_mort[as.character(dimnames(pred_mort)[[1]]) %in% 
                           species, , drop = FALSE]
  plot_dat <- data.frame(value = c(pred_mort), Species = factor(dimnames(pred_mort)[[1]], 
                                                                levels = species_levels), w = rep(params@w, each = length(species)))
  p <- ggplot(plot_dat) + geom_line(aes(x = w, y = value, colour = Species, 
                                        linetype = Species, size = Species))
  linesize <- rep(0.8, length(params@linetype))
  names(linesize) <- names(params@linetype)
  linesize[highlight] <- 1.6
  p <- p + scale_x_continuous(name = "Size [g]", trans = "log10") + 
    scale_y_continuous(name = "Predation mortality [1/year]", 
                       limits = c(0, max(plot_dat$value))) + scale_colour_manual(values = params@linecolour) + 
    scale_linetype_manual(values = params@linetype) + scale_size_manual(values = linesize)
  return(p)
}