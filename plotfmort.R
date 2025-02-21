
plotfmort=function (object, species = NULL, time_range, highlight = NULL, 
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
 getfmort= function (object, effort, time_range, drop = TRUE) 
  {
    if (is(object, "MizerParams")) {
      params <- validParams(object)
      if (missing(effort)) {
        effort <- params@initial_effort
      }
      n <- params@initial_n
      n_pp <- params@initial_n_pp
      n_other <- params@initial_n_other
      no_gears <- dim(params@catchability)[[1]]
      f <- get(params@rates_funcs$FMort)
      if (length(dim(effort)) == 2) {
        f_mort <- t(apply(effort, 1, function(x) f(params, 
                                                   n = n, n_pp = n_pp, n_other = n_other, effort = x, 
                                                   t = t, e_growth = getEGrowth(params, n = n, n_pp = n_pp, 
                                                                                n_other = n_other, t = t), pred_mort = getPredMort(params, 
                                                                                                                                   n = n, n_pp = n_pp, n_other = n_other, time_range = t))))
        dim(f_mort) <- c(dim(effort)[[1]], dim(params@initial_n))
        dimnames(f_mort) <- c(list(time = dimnames(effort)[[1]]), 
                              dimnames(params@initial_n))
        return(f_mort)
      }
      else if (length(effort) == 1) {
        fmort <- f(params, n = n, n_pp = n_pp, n_other = n_other, 
                   effort = rep(effort, no_gears), t = t, e_growth = getEGrowth(params, 
                                                                                n = n, n_pp = n_pp, n_other = n_other, t = t), 
                   pred_mort = getPredMort(params, n = n, n_pp = n_pp, 
                                           n_other = n_other, time_range = t))
        dimnames(fmort) <- dimnames(params@metab)
        return(fmort)
      }
      else if (length(effort) == no_gears) {
        fmort <- f(params, n = n, n_pp = n_pp, n_other = n_other, 
                   effort = effort, t = t, e_growth = getEGrowth(params, 
                                                                 n = n, n_pp = n_pp, n_other = n_other, t = t), 
                   pred_mort = getPredMort(params, n = n, n_pp = n_pp, 
                                           n_other = n_other, time_range = t))
        dimnames(fmort) <- dimnames(params@metab)
        return(fmort)
      }
      else {
        stop("Invalid effort argument")
      }
    }
    else {
      sim <- object
      if (missing(time_range)) {
        time_range <- dimnames(sim@effort)$time
      }
      time_elements <- get_time_elements(sim, time_range)
      f_mort <- getFMort(sim@params, sim@effort)
      return(f_mort[time_elements, , , drop = drop])
    }
  }
  f <- getfmort(object, time_range = time_range, drop = FALSE)
  if (length(dim(f)) == 3) {
    f <- apply(f, c(2, 3), mean)
  }
  if (is.null(species)) {
    species <- dimnames(params@initial_n)$sp[!is.na(params@A)]
  }
  species_levels <- c(as.character(params@species_params$species), 
                      "Background", "Resource", "Total")
  f <- f[as.character(dimnames(f)[[1]]) %in% species, , drop = FALSE]
  plot_dat <- data.frame(value = c(f), Species = factor(dimnames(f)[[1]], 
                                                        levels = species_levels), w = rep(params@w, each = length(species)))
  linesize <- rep(0.8, length(params@linetype))
  names(linesize) <- names(params@linetype)
  linesize[highlight] <- 1.6
  p <- ggplot(plot_dat) + geom_line(aes(x = w, y = value, colour = Species, 
                                        linetype = Species, size = Species))
  p <- p + scale_x_continuous(name = "Size [g]", trans = "log10") + 
    scale_y_continuous(name = "Fishing mortality [1/Year]", 
                       limits = c(0, max(plot_dat$value))) + scale_colour_manual(values = params@linecolour) + 
    scale_linetype_manual(values = params@linetype) + scale_size_manual(values = linesize)
  return(p)
}
