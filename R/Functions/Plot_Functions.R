# Functions for plotting results.

PlotBar <- function(df, missLevel, traitType, title, upperThresh) {
  
  # Function for plotting a bar plot depicting error rates for different methods at a certain missingness level.
  
  # df = dataframe containing error rate data
  # traitType = either "numeric" or "categorical".
  # missLevel = missingness level (corresponds to row to subset)
  # title = to include in title of plot
  # upperThresh = upper range of y-axis
  
  # Subset row according to missingness level.
  dfSub <- df[missLevel, ]
  # Remove SE columns for now.
  SE <- grep("_SE", colnames(dfSub))
  dfError <- dfSub[, -SE]
  # Remove missingness level col.
  dfError$missingness_level <- NULL
  # Create df of SE values.
  dfSE <- dfSub[, SE]
  
  # Transpose dfs.
  dfError <- t(dfError)
  dfSE <- t(dfSE)
  # Bind dfs.
  dfPlot <- cbind(dfError, dfSE)
  # Rename cols.
  colnames(dfPlot) <- c("error", "SE")
  
  # Let's create a grouping variable for adding colour to the plots later on. First, get the row names.
  methods <- rownames(dfPlot)
  # Identify elements with method names and replace elements with grouping var.
  methods[grep("KNN", methods)] <- "KNN"
  methods[grep("MICE", methods)] <- "MICE"
  methods[grep("RF", methods)] <- "RF"
  # Convert dfPlot to df.
  dfPlot <- as.data.frame(dfPlot)
  
  # Create a column to hold the group information.
  dfPlot$group <- methods
  # Create a method column for the x-axis.
  dfPlot$method <- rownames(dfPlot)
  
  # Create title of plot.
  plotTitle <- paste(title, missLevel, "BarPlot", ".tiff", sep = "")
  # First, let's order the df according to group so colours are consistent for methods for different traits.
  dfPlot <- dfPlot[order(dfPlot$group), ]
  
  if(traitType == "numeric") {
    
    # Change MeanMode to Mean.
    dfPlot[dfPlot == "MeanMode"] <- "Mean"
    # Convert group to factor variable to ensure correct plotting order. 
    dfPlot$group <- factor(dfPlot$group, levels = c("Mean", "KNN", "RF", "MICE"))
    # Create a vector to hold colours for grouping var.
    cols <- c("#F0E442", "#56B4E9", "#009E73", "#CC79A7")
    # Create a vector to hold colours for SE bars (need one for each error rate).
    SEcolours <- c(rep("#F0E442", 4), "#56B4E9", rep("#009E73", 4), rep("#CC79A7", 4))
    # Fix order so trait-only methods are first.
    methodOrd <- c("Mean", "KNN", "KNN_COI", "KNN_CMOS", "KNN_RAG1", "RF", "RF_COI", "RF_CMOS", "RF_RAG1", "MICE", "MICE_COI", "MICE_CMOS", "MICE_RAG1")
    # Fix names as gg plots alphabetically within group.
    dfPlot <- dfPlot[match(methodOrd, dfPlot$method),]
    # Fix order according to factor.
    dfPlot$method <- factor(dfPlot$method, levels = dfPlot$method[order(dfPlot$group)])
    dfPlot$method  # notice the changed order of factor levels
    
    # Create a vector to hold colours for patterns.
    pattCols <- c("#F0E442", "red", "white", "blue")
    
    theme_set(theme_classic())
    
    tiff(plotTitle, units="in", width=12, height=14, res=600)
    
    plot <- ggplot(data = dfPlot, aes(method, y = error)) +
      geom_col_pattern(
        aes(pattern = group, pattern_fill = group, fill = group),
        colour = "black",
        pattern_angle = 45,
        pattern_density = 0.1,
        pattern_spacing = 0.01,
        position = position_dodge2(preserve = 'single')) +
      #geom_bar(stat = "identity") +
      scale_pattern_manual(values = c("none", "stripe", "circle", "crosshatch"), 
                           guide = guide_legend(title = "Legend", override.aes = list(fill = cols))) +
      scale_colour_manual(values = pattCols) +
      scale_fill_manual(values = cols) + ylab("") + xlab("") + 
      geom_errorbar(aes(ymin = error - SE, ymax = error + SE, width = .2)) +
      ylim(0, upperThresh) +
      labs(title = plotTitle,
           subtitle = "",
           caption = "") +
      theme(axis.text.y = element_text(vjust = 0.6, size = 32, face = "bold"), axis.text.x = element_text(vjust = 0.6, size = 28, angle = 65, face = "bold"))
    
    print(plot)
    
  } else if(traitType == "categorical") {
    
    # Change MeanMode to Mode.
    dfPlot[dfPlot == "MeanMode"] <- "Mode"
    # Convert to factor so we can plot the error rates according to colour.
    dfPlot$method <- factor(dfPlot$method, levels = dfPlot$method[order(dfPlot$group)])
    # Convert group to factor variable to ensure correct plotting order. 
    dfPlot$group <- factor(dfPlot$group, levels = c("Mode", "KNN", "RF", "MICE"))
    # Create a vector to hold colours for grouping var.
    cols <- c("#F0E442", "#56B4E9", "#009E73", "#CC79A7")
    # Create a vector to hold colours for SE bars (need one for each error rate).
    SEcolours <- c(rep("#F0E442", 4), "#56B4E9", rep("#009E73", 4), rep("#CC79A7", 4))
    # Fix order so trait-only methods are first.
    methodOrd <- c("Mode", "KNN", "KNN_COI", "KNN_CMOS", "KNN_RAG1", "RF", "RF_COI", "RF_CMOS", "RF_RAG1", "MICE", "MICE_COI", "MICE_CMOS", "MICE_RAG1")
    # Order dfPlot by methodOrd.
    dfPlot <- dfPlot[match(methodOrd, dfPlot$method),]
    # Fix order according to factor.
    dfPlot$method <- factor(dfPlot$method, levels = dfPlot$method[order(dfPlot$group)])
    
    # Create a vector to hold colours for patterns.
    pattCols <- c("#F0E442", "red", "white", "blue")
    
    theme_set(theme_classic())
    
    tiff(plotTitle, units="in", width=12, height=14, res=600)
    
    plot <- ggplot(data = dfPlot, aes(method, y = error)) +
      geom_col_pattern(
        aes(pattern = group, pattern_fill = group, fill = group),
        colour = "black",
        pattern_angle = 45,
        pattern_density = 0.1,
        pattern_spacing = 0.01,
        position = position_dodge2(preserve = 'single')) +
      #geom_bar(stat = "identity") +
      scale_pattern_manual(values = c("none", "stripe", "circle", "crosshatch"), 
                           guide = guide_legend(title = "Legend", override.aes = list(fill = cols))) +
      scale_colour_manual(values = pattCols) +
      scale_fill_manual(values = cols) + ylab("") + xlab("") + 
      geom_errorbar(aes(ymin = error - SE, ymax = error + SE, width = .2)) +
      ylim(0, upperThresh) +
      labs(title = plotTitle,
           subtitle = "",
           caption = "") +
      theme(axis.text.y = element_text(vjust = 0.6, size = 32, face = "bold"), axis.text.x = element_text(vjust = 0.6, size = 28, angle = 65, face = "bold"))
    
    print(plot)
    
  }
  
  dev.off()
  
  return(plot)
  
}

PlotLine <- function(df, traitType, plotType, title, minError, maxError) {
  
  # Function for plotting error rates (Squamata specific version).
  
  # df = dataframe with error rates.
  # traitType = either "numeric" or "categorical".
  # plotType = either "KNN", "RF", "MICE" or "between". Between refers to error rates between methods (i.e. KNN vs. RF vs. MICE).
  # title = title to give to plot.
  # minError = lower range of y-axis scale for within method comparisons (to facilitate comparisons within trait)
  # maxError = upper range of y-axis scale for within method comparisons (to facilitate comparisons within trait)
  
  # If data type is numeric...
  if(traitType == "numeric") {
    
    # If plotType is "KNN"...
    if(plotType == "KNN") {
      
      # Plot line graph.
      plot <- plot_ly(df, x = ~missingness_level, y = ~KNN, error_y = ~list(array = KNN_SE, color = '#F0E442', opacity = 0.75), name = 'KNN', type = 'scatter', mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#F0E442')), line = list(color = c('#F0E442'), width = 9)) %>%
        add_trace(y = ~COI_KNN_PVR, name = 'KNN_COI', error_y = ~list(array = COI_KNN_PVR_SE, color = '#023FA5', opacity = 0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#023FA5')), line = list(color = c('#023FA5'), width = 9, dash = "dash")) %>%
        add_trace(y = ~RAG1_KNN_PVR, name = 'KNN_RAG1', error_y = ~list(array = RAG1_KNN_PVR_SE, color = '#0070A5', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#0070A5')), line = list(color = c('#0070A5'), width = 9, dash = "dot")) %>%
        add_trace(y = ~CMOS_KNN_PVR, name = 'KNN_CMOS', error_y = ~list(array = CMOS_KNN_PVR_SE, color = 'skyblue', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('skyblue')), line = list(color = c('skyblue'), width = 9, dash = "dashdot")) %>%
        layout(title = paste(title, "Method: KNN"), xaxis = list(title = "", ticks = 'outside', tickfont = list(size = 24), showgrid = T, titlefont = list(size = 31))) %>%
        layout(yaxis = list(titlefont = list(size = 31), range = c(minError, maxError), tickfont = list(size = 24), showgrid = T, title = "")) %>%
        layout(plot_bgcolor='white') %>%
        layout(showlegend = TRUE, legend = list(font = list(size = 24)))
      
      export(plot, file = paste(title, "Method_KNN", ".png", sep = ""))
      
    } else if(plotType == "RF") {
      
      # Plot line graph.
      plot <- plot_ly(df, x = ~missingness_level, y = ~RF, error_y = ~list(array = RF_SE, color = '#CC79A7', opacity=0.75), name = 'RF', type = 'scatter', mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#CC79A7')), line = list(color = c('#CC79A7'), width = 9)) %>%
        add_trace(y = ~COI_RF_PVR, name = 'RF_COI', error_y = ~list(array = COI_RF_PVR_SE, color = '#008A0F', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#008A0F')), line = list(color = c('#008A0F'), width = 9, dash = "dash")) %>%
        add_trace(y = ~RAG1_RF_PVR, name = 'RF_RAG1', error_y = ~list(array = RAG1_RF_PVR_SE, color = '#33B761', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#33B761')), line = list(color = c('#33B761'), width = 9, dash = "dot")) %>%
        add_trace(y = ~CMOS_RF_PVR, name = 'RF_CMOS', error_y = ~list(array = CMOS_RF_PVR_SE, color = '#A8DAB4', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#A8DAB4')), line = list(color = c('#A8DAB4'), width = 9, dash = "dashdot")) %>%
        layout(title = paste(title, "Method: RF"), xaxis = list(title = "", ticks = 'outside', tickfont = list(size = 24), showgrid = T, titlefont = list(size = 31))) %>%
        layout(yaxis = list(titlefont = list(size = 31), range = c(minError, maxError), tickfont = list(size = 24), showgrid = T, title = "")) %>%
        layout(plot_bgcolor='white') %>%
        layout(showlegend = TRUE, legend = list(font = list(size = 24)))
      
      export(plot, file = paste(title, "Method_RF", ".png", sep = ""))
      
    } else if(plotType == "MICE") {
      
      # Plot line graph.
      plot <- plot_ly(df, x = ~missingness_level, y = ~MICE, error_y = ~list(array = MICE_SE, color = '#009E73', opacity=0.75), name = 'MICE', type = 'scatter', mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#009E73')), line = list(color = c('#009E73'), width = 9)) %>%
        add_trace(y = ~COI_MICE_PVR, name = 'MICE_COI', error_y = ~list(array = COI_MICE_PVR_SE, color = "#8E063B", opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c("#8E063B")), line = list(color = c("#8E063B"), width = 9, dash = "dash")) %>%
        add_trace(y = ~RAG1_MICE_PVR, name = 'MICE_RAG1', error_y = ~list(array = RAG1_MICE_PVR_SE, color = '#F79CD4', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#F79CD4')), line = list(color = c('#F79CD4'), width = 9, dash = "dot")) %>%
        add_trace(y = ~CMOS_MICE_PVR, name = 'MICE_CMOS', error_y = ~list(array = CMOS_MICE_PVR_SE, color = '#F4B8C0', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#F4B8C0')), line = list(color = c('#F4B8C0'), width = 9, dash = "dashdot")) %>%
        layout(title = paste(title, "Method: MICE"), xaxis = list(title = "", ticks = 'outside', tickfont = list(size = 24), showgrid = T, titlefont = list(size = 31))) %>%
        layout(yaxis = list(titlefont = list(size = 31), range = c(minError, maxError), tickfont = list(size = 24), showgrid = T, title = "")) %>%
        layout(plot_bgcolor='white') %>%
        layout(showlegend = TRUE, legend = list(font = list(size = 24)))
      
      export(plot, file = paste(title, "Method_MICE", ".png", sep = ""))
      
    } else if(plotType == "between") {
      
      plot <- plot_ly(df, x = ~missingness_level, y = ~MeanMode, error_y = ~list(array = MeanMode_SE, color = '#F0E442', opacity=0.75), name = 'Mean', type = 'scatter', mode = 'lines+markers', symbol = 'none', marker = list(size = 16, color = c('#F0E442')), line = list(color = c('#F0E442'), width = 9)) %>%
        add_trace(y = ~KNN, name = 'KNN', error_y = ~list(array = KNN_SE, color = '#56B4E9', opacity=0.75), mode = 'lines+markers', symbol = 'none', marker = list(size = 30, color = c('#56B4E9')), line = list(color = c('#56B4E9'), width = 9)) %>%
        add_trace(y = ~RF, name = 'RF', error_y = ~list(array = RF_SE, color = '#009E73', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 20, color = c('#009E73')), line = list(color = c('#009E73'), width = 9)) %>%
        add_trace(y = ~MICE, name = 'MICE_PMM', error_y = ~list(array = MICE_SE, color = '#CC79A7', opacity=0.75), mode = 'lines+markers', symbol = 'square', marker = list(size = 18, color = c('#CC79A7')), line = list(color = c('#CC79A7'), width = 9)) %>%
        layout(title = paste(title, "Method comparison"), xaxis = list(title = "", ticks = 'outside', tickfont = list(size = 24), showgrid = T, titlefont = list(size = 31))) %>%
        layout(yaxis = list(titlefont = list(size = 31), tickfont = list(size = 24), showgrid = T, title = "")) %>%
        layout(plot_bgcolor='white') %>%
        layout(showlegend = TRUE, legend = list(font = list(size = 24)), itemsizing='trace')
      
      export(plot, file = paste(title, "Method comparison", ".png", sep = ""))
      
    }
    
  } else if(traitType == "categorical") {
    
    if(plotType == "KNN") {
      
      # Plot line graph.
      plot <- plot_ly(df, x = ~missingness_level, y = ~KNN, error_y = ~list(array = KNN_SE, color = '#F0E442', opacity=0.75), name = 'KNN', type = 'scatter', mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#F0E442')), line = list(color = c('#F0E442'), width = 9)) %>%
        add_trace(y = ~COI_KNN_PVR, name = 'KNN_COI', error_y = ~list(array = COI_KNN_PVR_SE, color = '#023FA5', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#023FA5')), line = list(color = c('#023FA5'), width = 9, dash = "dash")) %>%
        add_trace(y = ~RAG1_KNN_PVR, name = 'KNN_RAG1', error_y = ~list(array = RAG1_KNN_PVR_SE, color = '#0070A5', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#0070A5')), line = list(color = c('#0070A5'), width = 9, dash = "dot")) %>%
        add_trace(y = ~CMOS_KNN_PVR, name = 'KNN_CMOS', error_y = ~list(array = CMOS_KNN_PVR_SE, color = 'skyblue', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('skyblue')), line = list(color = c('skyblue'), width = 9, dash = "dashdot")) %>%
        layout(title = paste(title, "Method: KNN"), xaxis = list(title = "", ticks = 'outside', tickfont = list(size = 24), showgrid = T, titlefont = list(size = 31))) %>%
        layout(yaxis = list(titlefont = list(size = 31), range = c(minError, maxError), tickfont = list(size = 24), showgrid = T, title = "")) %>%
        layout(plot_bgcolor='white') %>%
        layout(showlegend = TRUE, legend = list(font = list(size = 24)))
      
      export(plot, file = paste(title, "Method_KNN", ".png", sep = ""))
      
    } else if(plotType == "RF") {
      
      # Plot line graph.
      plot <- plot_ly(df, x = ~missingness_level, y = ~RF, error_y = ~list(array = RF_SE, color = '#CC79A7', opacity=0.75), name = 'RF', type = 'scatter', mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#CC79A7')), line = list(color = c('#CC79A7'), width = 9)) %>%
        add_trace(y = ~COI_RF_PVR, name = 'RF_COI', error_y = ~list(array = COI_RF_PVR_SE, color = '#008A0F', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#008A0F')), line = list(color = c('#008A0F'), width = 9, dash = "dash")) %>%
        add_trace(y = ~RAG1_RF_PVR, name = 'RF_RAG1', error_y = ~list(array = RAG1_RF_PVR_SE, color = '#33B761', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#33B761')), line = list(color = c('#33B761'), width = 9, dash = "dot")) %>%
        add_trace(y = ~CMOS_RF_PVR, name = 'RF_CMOS', error_y = ~list(array = CMOS_RF_PVR_SE, color = '#A8DAB4', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#A8DAB4')), line = list(color = c('#A8DAB4'), width = 9, dash = "dashdot")) %>%
        layout(title = paste(title, "Method: RF"), xaxis = list(title = "", ticks = 'outside', tickfont = list(size = 24), showgrid = T, titlefont = list(size = 31))) %>%
        layout(yaxis = list(titlefont = list(size = 31), range = c(minError, maxError), tickfont = list(size = 24), showgrid = T, title = "")) %>%
        layout(plot_bgcolor='white') %>%
        layout(showlegend = TRUE, legend = list(font = list(size = 24)))
      
      export(plot, file = paste(title, "Method_RF", ".png", sep = ""))
      
    } else if(plotType == "MICE") {
      
      # Plot line graph.
      plot <- plot_ly(df, x = ~missingness_level, y = ~MICE, error_y = ~list(array = MICE_SE, color = '#009E73', opacity=0.75), name = 'MICE', type = 'scatter', mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#009E73')), line = list(color = c('#009E73'), width = 9)) %>%
        add_trace(y = ~COI_MICE_PVR, name = 'MICE_COI', error_y = ~list(array = COI_MICE_PVR_SE, color = "#8E063B", opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c("#8E063B")), line = list(color = c("#8E063B"), width = 9, dash = "dash")) %>%
        add_trace(y = ~RAG1_MICE_PVR, name = 'MICE_COI', error_y = ~list(array = RAG1_MICE_PVR_SE, color = '#F79CD4', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#F79CD4')), line = list(color = c('#F79CD4'), width = 9, dash = "dot")) %>%
        add_trace(y = ~CMOS_MICE_PVR, name = 'MICE_COI', error_y = ~list(array = CMOS_MICE_PVR_SE, color = '#F4B8C0', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 16, color = c('#F4B8C0')), line = list(color = c('#F4B8C0'), width = 9, dash = "dashdot")) %>%
        layout(title = paste(title, "Method: MICE"), xaxis = list(title = "", ticks = 'outside', tickfont = list(size = 24), showgrid = T, titlefont = list(size = 31))) %>%
        layout(yaxis = list(titlefont = list(size = 31), range = c(minError, maxError), tickfont = list(size = 24), showgrid = T, title = "")) %>%
        layout(plot_bgcolor='white') %>%
        layout(showlegend = TRUE, legend = list(font = list(size = 24)))
      
      export(plot, file = paste(title, "Method_MICE", ".png", sep = ""))
      
      
    } else if(plotType == "between") {
      
      plot <- plot_ly(df, x = ~missingness_level, y = ~MeanMode, error_y = ~list(array = MeanMode_SE, color = '#F0E442', opacity=0.75), name = 'Mode', type = 'scatter', mode = 'lines+markers', symbol = 'none', marker = list(size = 10, color = c('#F0E442')), line = list(color = c('#F0E442'), width = 9)) %>%
        add_trace(y = ~KNN, name = 'KNN', error_y = ~list(array = KNN_SE, color = '#56B4E9', opacity=0.75), mode = 'lines+markers', symbol = 'none', marker = list(size = 30, color = c('#56B4E9')), line = list(color = c('#56B4E9'), width = 9)) %>%
        add_trace(y = ~RF, name = 'RF', error_y = ~list(array = RF_SE, color = '#009E73', opacity=0.75), mode = 'lines+markers', symbol = 'circle', marker = list(size = 20, color = c('#009E73')), line = list(color = c('#009E73'), width = 9)) %>%
        add_trace(y = ~MICE, name = 'MICE_LR', error_y = ~list(array = MICE_SE, color = '#CC79A7', opacity=0.75), mode = 'lines+markers', symbol = 'square', marker = list(size = 18, color = c('#CC79A7')), line = list(color = c('#CC79A7'), width = 9)) %>%
        layout(title = paste(title, "Method comparison"), xaxis = list(title = "", ticks = 'outside', tickfont = list(size = 24), showgrid = T, titlefont = list(size = 31))) %>%
        layout(yaxis = list(titlefont = list(size = 31), tickfont = list(size = 24), showgrid = T, title = "")) %>%
        layout(plot_bgcolor='white') %>%
        layout(showlegend = TRUE, legend = list(font = list(size = 24)), itemsizing='trace')
      
      export(plot, file = paste(title, "Method comparison", ".png", sep = ""))
      
    }
    
  }
  
}

PlotPhySig <- function(df, title) {
  
  # Function for visualizing phylogenetic signal data (barplot).
  # df = dataframe containing phylogenetic signal measures. Must also contain a "gene" and "trait" column.
  
  # Create a vector to hold colours for grouping var (gene).
  cols <- c("#023FA5", "#7D87B9", "#E2E2E2")
  
  # Set the theme.
  theme_set(theme_light())
  
  # Create title.
  plotTitle <- paste(title, ".tiff", sep = "")
  # Write to .tiff file.
  tiff(plotTitle, units="in", width=14, height=8, res=600)
  
  # Barplot.
  plot <- ggplot(data = df, aes(x = trait, y = phylogenetic_signal, fill = gene)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = cols) + ylab("") + xlab("") + 
    theme(axis.text.y = element_text(vjust = 0.6, size = 22, face = "bold"), axis.text.x = element_text(vjust = 0.6, size = 20, face = "bold")) + coord_flip()
  
  # Print to tiff.
  print(plot)
  
  # Turn off dev.
  dev.off()
  
}

RatioError <- function(df, method) {
  
  # Function for dividing baseline error rate by all other method types.
  # df  = dataframe containing error and empty ratio column.
  # method = method type. One of "KNN", "RF", or "MICE".
  
  # Which row is the baseline method?
  baseline <- df[which(df$method == method), "error"]
  # Identify other rows for this method.
  rowIndex <- setdiff(grep(method, df$method), which(df$method == method))
  # Divide these rows by baseline error rate and enter into ratio column.
  df[rowIndex, "ratio"] <- baseline/(df[rowIndex, "error"])
  # Return modified df.
  return(df)
  
}
