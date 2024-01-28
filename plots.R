library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(data.table)

library(ggthemes)

calculate_epi_params <- function(params) {
  c(
    N = params[["N"]],
    R0 = params[["beta"]] / params[["gamma"]],
    latent = 1 / params[["sigma"]],
    infectious = 1/ params[["gamma"]]
  )
}

strain_info <- function(params) {
  return("")
}

plot_SEIR <- function(
    dist,
    style = c("combined", "exploded", "incidence")
) {
  theme_set(theme_minimal())
  
  data <- dist$overview
  
  if (style == "combined") {
    # melt the data to be long, for ggplot2
    long_data <- pivot_longer(data, cols = c("S", "E", "I", "R"))
    
    result <- ggplot(data=long_data, aes(x=time, y=value, colour=name)) + 
      geom_line(linewidth=0.8) +
      event_markers(dist$events) +

      
      # manually set legend
      scale_colour_manual(
        values = c("seagreen3", "black", "red", "deepskyblue3"),
        breaks=c("S", "E", "I", "R"),
        name = "Type",
        labels = c("susceptible", "exposed", "infected", "removed")
      ) +

      # general theming
      theme(legend.title=element_blank(), legend.position = c(0.88, 0.5)) +
      labs(x = "day (t)", y = "quantity")
    
    return(result)
  }
  
  else if (style == "exploded") {
    p1 <- ggplot(data, aes(x=time, y=S)) + geom_line(linewidth = 0.8, colour="seagreen3") +
      labs(x = "day (t)", y = "", title = "susceptible") +
      event_markers(dist$events) +
      theme(plot.title = element_text(size=10))
    
    p2 <- ggplot(data, aes(x=time, y=R)) + geom_line(linewidth = 0.8, colour="deepskyblue3") +
      labs(x = "day (t)", y = "", title = "removed") +
      event_markers(dist$events) +
      theme(plot.title = element_text(size=10))
    
    p3 <- ggplot(data, aes(x=time, y=E)) + geom_line(linewidth = 0.8, colour="black") +
      labs(x = "day (t)", y = "", title = "exposed") +
      event_markers(dist$events) +
      theme(plot.title = element_text(size=10))
    
    p4 <- ggplot(data, aes(x=time, y=I)) + geom_line(linewidth = 0.8, colour="red") +
      labs(x = "day (t)", y = "", title = "infectious") +
      event_markers(dist$events) +
      theme(plot.title = element_text(size=10))
    
    grid.arrange(p1, p2, p3, p4)
  }
  
  else if (style == "incidence") {
    result <- ggplot() +
      geom_line(data=data, aes(x=time, y=incidence), linewidth = 0.8) +
      event_markers(dist$events) +
      labs(x = "day (t)", y = "incidence")
    
    return(result)
  }
}

plot_SEIR_incidences <- function(default, additional, parameters) {
  epi_params <- calculate_epi_params(parameters)
  
  title <- sprintf("SEIR model\nwith N = %s, R_0 = %#.2f, lat = %#.1f, inf = %#.1f\n(1 default, %d additional)",
                   epi_params[["N"]],
                   epi_params[["R0"]],
                   epi_params[["latent"]],
                   epi_params[["infectious"]],
                   length(additional)
  )
  
  default$overview$.id = 0
  
  # combine the multiple df's into one single one with a .id col
  data <- rbindlist(lapply(additional, \(x) x$overview), idcol=TRUE)
  data <- rbind(data, default$overview)
  
  result <- ggplot() + 
    geom_line(data[.id != 0], mapping = aes(x=time, y=incidence, group=.id), linewidth = 0.5, alpha = 0.1) +
    geom_line(data[.id == 0], mapping = aes(x=time, y=incidence), linewidth = 1.3, colour="seagreen2") +
    labs(x = "day (t)", y = "incidence", title = title)
  
  result
}

plot_multiple_SEIR <- function(
    dist,
    style = c("combined", "exploded", "incidence")
) {
  
  data <- dist$overview
  
  epi_params = calculate_epi_params(data$params)
  
  title <- paste(
    c(
      "SEIR model with",
      lapply(
        names(dist$strains),
        function(x) {
          paste(x, ": ", dist$strains[[x]]$params(0)$label, sep="")
        }
      )
    ),
    collapse = "\n"
  )
  
  if (style == "combined") {
    # melt the data to be long, for ggplot2
    long_data <- pivot_longer(data, cols = c("S", "E", "I", "R"))
    
    result <- ggplot(long_data, aes(x=time, y=value, colour=name, linetype=strain)) + 
      event_markers(dist$events) +
      geom_line(linewidth=0.8) +
      
      # manually set legend
      scale_colour_manual(
        values = c("seagreen3", "black", "red", "deepskyblue3"),
        breaks=c("S", "E", "I", "R"),
        name = "Type",
        labels = c("susceptible", "exposed", "infectious", "removed")
      ) +
      
      # general theming
      theme(legend.title=element_blank(), legend.position = "right") +
      labs(x = "day (t)", y = "quantity", title = title)
    
    return(result)
  }
  
  else if (style == "exploded") {
    plot_options <- list(
      list(
        cols = c("S.AB", "S.A", "S.B"),
        label = "susceptible"
      ),
      list(
        cols = c("R.AB", "R.A", "R.B"),
        label = "recovered"
      ),
      list(
        cols = c("E.A1", "E.A2", "E.B1", "E.B2"),
        label = "exposed"
      ),
      list(
        cols = c("I.A1", "I.A2", "I.B1", "I.B2"),
        label = "infectious"
      )
    )
    
    plots <- list()
    
    for (i in seq_along(plot_options)) {
      p <- plot_options[[i]]
      filtered_data <- dist$long_data[dist$long_data$name %in% p$cols, ]
      
      plots[[i]] <- 
        ggplot(filtered_data, aes(x = time, y = value, colour = name, linetype=strain)) + 
        event_markers(dist$events) +
        geom_line(linewidth = 0.8) +
        labs(x = "day (t)", y = "", title = p$label) +
        theme(plot.title = element_text(size=12),
              legend.title=element_blank(),
              legend.position = "bottom"
        )
    }
    
    grid.arrange(grobs = plots,
                 top = textGrob(title, gp=gpar(fontsize=12), x = 0, hjust = 0))
  }
  
  else if (style == "incidence") {
    result <- ggplot() +
      event_markers(dist$events) +
      geom_line(data=data, aes(x=time, y=i, colour=strain), linewidth = 0.8) +
      labs(x = "day (t)", y = "incidence", title = title)
    
    return(result)
  }
}

event_markers <- function(events) {
  geom_vline(data=events, aes(xintercept=time), alpha=0.3, linetype="dashed")
}
