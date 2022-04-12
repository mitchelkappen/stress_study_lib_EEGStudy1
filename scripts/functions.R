# Plot functions for:
#   3.0_Audio_SelfReport_analysis_visualisation.R
#   3.1_IBI_analysis_visualisation.R

# One general plotting function
audio_pretty_plot <-
  function(emmean_dataframe, title){
    ggplot(emmean_dataframe, aes(x=condition, y=emmean, colour = condition)) +
      geom_point(size = 5) + 
      geom_line(aes(group = 1),size = 1, colour = "black", linetype = "dotted")+
      geom_errorbar(width=.25, size = 1.4, aes(ymin=emmean-SE, ymax=emmean+SE))+
      labs(y = title, x = "Feedback condition")+
      scale_colour_manual(values=cbPalette)+
      plot_theme_apa()
  }
# One general theme to clean up code
plot_theme_apa <-
  function(...){
    theme(
      # legend.key.size=unit(1.3, 'cm'),
      # legend.text=element_text(size=13),
      legend.position = "none",
      plot.title = element_text(size=rel(2)),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.grid.major.y = element_line( size=.1, color="#dedede" ),
      axis.text.x=element_text(size=rel(2)),
      axis.title.y=element_text(size=rel(1.5)),
      axis.title.x = element_text(size=rel(1.5)))
  }