library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(data.table)

calculate_epi_params <- function(params) {
  c(
    N = params[["N"]],
    R0 = params[["beta"]] / params[["gamma"]],
    latent = 1 / params[["sigma"]],
    infectious = 1/ params[["gamma"]]
  )
}

plot_SEIR <- function(
    data,
    params,
    style = c("combined", "exploded", "incidence")
) {
  epi_params = calculate_epi_params(params)
  
  title <- sprintf("SEIR model\nwith N = %s, R_0 = %#.2f, lat = %#.1f, inf = %#.1f",
                   epi_params[["N"]],
                   epi_params[["R0"]],
                   epi_params[["latent"]],
                   epi_params[["infectious"]]
  )
  
  data <- data$overview
  
  if (style == "combined") {
    # melt the data to be long, for ggplot2
    long_data <- pivot_longer(data, cols = c("S", "E", "I", "R"))
    
    result <- ggplot(long_data, aes(x=time, y=value, colour=name)) + 
      geom_line(linewidth=0.8) +
      
      # manually set legend
      scale_colour_manual(
        values = c("seagreen3", "black", "red", "deepskyblue3"),
        breaks=c("S", "E", "I", "R"),
        name = "Type",
        labels = c("susceptible", "exposed", "infected", "removed")
      ) +

      # general theming
      theme(legend.title=element_blank(), legend.position = c(0.88, 0.5)) +
      labs(x = "day (t)", y = "quantity", title = title)
    
    return(result)
  }
  
  else if (style == "exploded") {
    p1 <- ggplot(data, aes(x=time, y=S)) + geom_line(linewidth = 0.8, colour="seagreen3") +
      labs(x = "day (t)", y = "", title = "susceptible") +
      theme(plot.title = element_text(size=10))
    
    p2 <- ggplot(data, aes(x=time, y=R)) + geom_line(linewidth = 0.8, colour="deepskyblue3") +
      labs(x = "day (t)", y = "", title = "removed") +
      theme(plot.title = element_text(size=10))
    
    p3 <- ggplot(data, aes(x=time, y=E)) + geom_line(linewidth = 0.8, colour="black") +
      labs(x = "day (t)", y = "", title = "exposed") +
      theme(plot.title = element_text(size=10))
    
    p4 <- ggplot(data, aes(x=time, y=I)) + geom_line(linewidth = 0.8, colour="red") +
      labs(x = "day (t)", y = "", title = "infected") +
      theme(plot.title = element_text(size=10))
    
    grid.arrange(p1, p2, p3, p4, top = textGrob(title, gp=gpar(fontsize=16), x = 0.09, hjust = 0))
  }
  
  else if (style == "incidence") {
    result <- ggplot() +
      geom_line(data=data, aes(x=time, y=incidence), linewidth = 0.8) +
      labs(x = "day (t)", y = "incidence", title = title)
    
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
