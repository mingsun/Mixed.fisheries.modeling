plotspectra<- function (sim, species = NULL, time_range = max(as.numeric(dimnames(sim@n)$time)), 
                        min_w = min(sim@params@w)/100, ylim = c(NA, NA), power = 1,
                        biomass = TRUE, print_it = TRUE, total = FALSE, plankton = TRUE, 
                        background = TRUE,highlight = NULL, ...) 
{library(ggplot2)
  if (missing(power)) {
    power <- as.numeric(biomass)
  }
  time_range <- range(as.numeric(time_range))
  sim_times <- as.numeric(dimnames(sim@effort)[[1]])
  sim_time_range <- range(sim_times)
  if ((time_range[1] < sim_time_range[1]) | (time_range[2] > 
                                             sim_time_range[2])) 
    stop("Time range is outside the time range of the model")
  time_elements <- (sim_times >= time_range[1]) & (sim_times <= 
                                                     time_range[2])
  names(time_elements) <- dimnames(sim@effort)$time
  n <- apply(sim@n[time_elements, , , drop = FALSE], 
             c(2, 3), mean)
  n_pp <- apply(sim@n_pp[time_elements, , drop = FALSE], 
                2, mean)
  params<-sim@params
  if (total) {
    fish_idx <- (length(params@w_full) - length(params@w) + 
                   1):length(params@w_full)
    total_n <- n_pp
    total_n[fish_idx] <- total_n[fish_idx] + colSums(n)
    total_n <- total_n * params@w_full^power
  }
  if (is.null(species)) {
    species <- params@species_params$species[!is.na(params@A)]
  }
  if (power %in% c(0, 1, 2)) {
    y_label = c("Number density [1/g]", "Biomass density", 
                "Biomass density [g]")[power + 1]
  }
  else {
    y_label = paste0("Number density * w^", power)
  }
  n <- sweep(n, 2, params@w^power, "*")
  spec_n <- n[as.character(dimnames(n)[[1]]) %in% species, 
              , drop = FALSE]
  plot_dat <- data.frame(value = c(spec_n), Species = as.factor(dimnames(spec_n)[[1]]), 
                         w = rep(params@w, each = dim(spec_n)[[1]]))
  if (plankton) {
    max_w <- min(params@species_params$w_mat)
    if (is.na(max_w)) {
      max_w <- Inf
    }
    plankton_sel <- params@w_full >= min_w & params@w_full < 
      max_w
    w_plankton <- params@w_full[plankton_sel]
    plank_n <- n_pp[plankton_sel] * w_plankton^power
    plot_dat <- rbind(plot_dat, data.frame(value = c(plank_n), 
                                           Species = "Plankton", w = w_plankton))
  }
  if (total) {
    plot_dat <- rbind(plot_dat, data.frame(value = c(total_n), 
                                           Species = "Total", w = params@w_full))
  }
  plot_dat <- plot_dat[(plot_dat$value > 0) & (plot_dat$w >= 
                                                 min_w), ]
  if (!is.na(ylim[1])) {
    plot_dat <- plot_dat[plot_dat$value < ylim[1], ]
  }
  if (is.na(ylim[2])) {
    ylim[2] <- 1e-20
  }
  plot_dat <- plot_dat[plot_dat$value > ylim[2], ]
  linesize <- rep(0.8, length(params@linetype))
  names(linesize) <- names(params@linetype)
  linesize[highlight] <- 1.6
  p <- ggplot(plot_dat, aes(x = w, y = value)) + scale_x_continuous(name = "Size [g]", 
                                                                    trans = "log10", breaks =scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_continuous(trans = "log10", 
                       breaks =scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)),
                       name = y_label)+ 
    scale_colour_manual(values = params@linecolour) + 
    scale_linetype_manual(values = params@linetype)+ scale_size_manual(values = linesize)
  if (background) {
    back_n <- n[is.na(params@A), , drop = FALSE]
    plot_back <- data.frame(value = c(back_n), Species = as.factor(dimnames(back_n)[[1]]), 
                            w = rep(params@w, each = dim(back_n)[[1]]))
    plot_back <- plot_back[(plot_back$value > 0) & (plot_back$w >= 
                                                      min_w), ]
    if (!is.na(ylim[1])) {
      plot_back <- plot_back[plot_back$value < ylim[1], 
      ]
    }
    plot_back <- plot_back[plot_back$value > ylim[2], ]
    p <- p + geom_line(aes(group = Species, size = Species), colour =NA, 
                       data = plot_back)
  }
  if ((length(species) + plankton + total) >51) {
    p <- p + geom_line(aes(group = Species, size = Species))
  }
  else {
    p <- p + geom_line(aes(colour = Species, linetype = Species, size = Species))
  }
  if (print_it) {
    print(p)}
  return(p)
}
