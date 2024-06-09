library(shiny)
library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(reshape2)
library(coin)
library(plotly)
library(bslib)
library(shinyjs)


AGGREGATE_COLS <-c( "E11.5_XY_1","E11.5_XY_2","E11.5_XY_3","E11.5_XX_1","E11.5_XX_2","E11.5_XX_3",
                    "E12.5_XY_1","E12.5_XY_2","E12.5_XY_3","E12.5_XX_1","E12.5_XX_2","E12.5_XX_3",
                    "E13.5_XY_1","E13.5_XY_2","E13.5_XY_3","E13.5_XX_1","E13.5_XX_2","E13.5_XX_3",
                    "E15.5_XY_1","E15.5_XY_2","E15.5_XY_3","E15.5_XX_1","E15.5_XX_2","E15.5_XX_3")
COLS_TO_REMOVE <- c(8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 21, 22)

# Function to remove unwanted columns and standardize column names
clean_matrix <- function(tpm_matrix) {
  # Remove the columns from the dataframe
  tpm_matrix <- tpm_matrix[, -COLS_TO_REMOVE]
  
  # Modify column names
  cols <- colnames(tpm_matrix)
  cols_to_modify <- cols[!grepl('gene_id', cols)]
  new_cols <- sub('.*E1', 'E1', cols_to_modify)
  colnames(tpm_matrix)[!grepl('gene_id', cols)] <- new_cols
  cols <- colnames(tpm_matrix)
  pattern <- "(_[XY]{2}).*?(_[123])"
  new_cols <- sub(pattern, "\\1\\2", cols)
  colnames(tpm_matrix) <- new_cols
  
  return(tpm_matrix)
}

# Function to melt data to long format for ggplot2
melt_data <- function(tpm_matrix) {
  group_means_melted <- melt(tpm_matrix, id.vars = "gene_id", 
                             variable.name = "Sample", value.name = "Expression")
  group_means_melted$Stage <- sub("_.*", "", group_means_melted$Sample)
  group_means_melted$Cell_Type <- ifelse(grepl("XY", group_means_melted$Sample), "Sertoli", "Granulosa")
  
  return(group_means_melted)
}

TPM_FILE <- "230606_ALL_rsem.merged.gene_tpm.csv"

tpm_matrix <- read.table(TPM_FILE, sep = ',', header = TRUE)
tpm_matrix <- clean_matrix(tpm_matrix)

filtered_df <- tpm_matrix[!grepl("Rik$", tpm_matrix$gene_id), ]
ord_mtx = filtered_df[order(filtered_df$gene_id),]

# Calculate mean expression values for each group
group_means <- aggregate(tpm_matrix[,AGGREGATE_COLS],
                         by = list(tpm_matrix$gene_id), mean)

# Melt data to long format for ggplot2
group_means_melted <- melt_data(tpm_matrix)

# Define UI ----
ui <- fillPage(
  titlePanel("RNA Seq Results"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Show gene specific expression (in TPM)."),
      
      # Allow multiple gene selections
      selectizeInput("gene_name",
                     label = "Choose gene name(s):",
                     choices = ord_mtx$gene_id,
                     selected = c("Amh"),
                     multiple = TRUE),
      width = 2
    ),
    
    mainPanel(
      htmlOutput('gridScatter'),
      width = 10  # Set width to 12 to occupy entire page
    )
  ),
  
  tags$head(
    tags$style(HTML("
      #gridScatter {
        height: 90vh;
        overflow-y: auto;
      }
    "))
  )
)

# Define server logic ----
server <- function(input, output) {
  sort_data <- function(data, genes_names){
    # Filter data for selected genes
    data <- group_means_melted[group_means_melted$gene_id %in% genes_names, ]
    data$Day <- substring(data$Stage, 2)
    data <- data %>%
      group_by(Day, Cell_Type, gene_id) %>%
      mutate(Expression_mean = mean(Expression, na.rm = TRUE))
    return(data)
  }
  get_color_pallete <- function(genes){
    num_colors <- length(unique(genes))
    # Generate colors from "GnBu" and "PuRd"
    colors_even <- darken(brewer.pal(num_colors, "Blues"), 0.2)
    colors_odd <- darken(brewer.pal(num_colors, "RdPu"), 0.2)
    combined_colors <- list()
    for (i in 1:num_colors) {
      combined_colors[[2*i - 1]] <- colors_odd[i]
      combined_colors[[2*i]] <- colors_even[i]
    }
    return(combined_colors)
  }
  get_pval <- function(data){
    p_df <- data.frame(
      gene_id = character(),
      p_value = numeric(),
      Day = numeric(),
      stringsAsFactors = FALSE # Optional: Avoid automatic conversion to factors
    )
    # # show difference between cell type at each time point
    for (day in unique(data$Day)) {
      p_values <- data %>%
        filter(Day == day) %>%
        group_by(gene_id) %>%
        summarize(p_value = t.test(Expression ~ Cell_Type)$p.value)  %>%
        mutate(Day = day)
      p_df <- rbind(p_df, p_values)
    }
    p_df <- p_df[order(p_df$Day), ]
    alpha_levels <- c(0.05, 0.01)
    p_df$label <- ifelse(p_df$p_value >= alpha_levels[1], "-ns-",
                         ifelse(p_df$p_value >= alpha_levels[2], "-*-", "-**-"))
    return(p_df)
  }
  colors <- function(){
    return(c("#E7BCDE", "#86B6F6"))
  }
  output$gridScatter <- renderUI({
    data <- sort_data(data, input$gene_name)
    # Combine the Cell_Type and gene_id columns to create a new column for grouping
    data$group <- interaction(data$gene_id, data$Cell_Type)
    data <- data[order(data$Day), ]
    unique_days <- as.numeric(unique(data$Day))
    # make genes names italic
    subtitles <- paste0("<i>", input$gene_name, "</i>")
    # Reorder data based on sorted_gene_id otherwise the plots will be alphabetical and the titles not...
    sorted_gene_id <- unique(data$gene_id)[order(match(unique(data$gene_id), input$gene_name))]
    data <- data[order(match(data$gene_id, sorted_gene_id)), ]
    plots <- c()
    plots <- lapply(unique(data$gene_id), function(gene) {
      gene_data <- data[data$gene_id == gene,]
      p_values <- get_pval(gene_data)
      gene_data <- gene_data %>%
        left_join(p_values, by = c("Day", "gene_id")) %>%
        select(-Day, -gene_id)      
      gene_data$maxExp <- 1.01 * max(gene_data$Expression)
      scatter_plot <- plot_ly(x = ~gene_data$Day, y = ~gene_data$Expression,
                              type = 'scatter', mode = 'markers', color = ~gene_data$Cell_Type, colors = colors(), 
                              legendgroup=~gene_data$Cell_Type, 
                              showlegend = ifelse(gene == unique(data$gene_id)[1], TRUE, FALSE)) 
      scatter_plot <- scatter_plot %>% 
        add_trace(y = ~gene_data$Expression_mean,
                  mode = 'lines', name = ~gene_data$Cell_Type, color = ~gene_data$Cell_Type,
                  showlegend = FALSE) %>% 
        layout(xaxis = list(title = list(text ='Embryonic day')), 
               yaxis = list(title = list(text ='TPM')), 
               height = 800) %>%
        add_trace(y = ~gene_data$maxExp, type = 'scatter',
                  mode = 'text', text = ~gene_data$label, textposition = 'middle',
                  textfont = list(color = '#000000', size = 10), showlegend = FALSE)
    })
    num_plots <- length(plots)
    num_cols <- 2  # Number of columns in the subplot grid
    
    # Calculate the number of rows dynamically
    num_rows <- ceiling(num_plots / num_cols)
    fig <- subplot(plots, nrows = num_rows, titleY= TRUE, titleX= TRUE, shareX = TRUE)
    #set y coordinate of titles
    title_y <- rev(seq(0, 1, length.out = num_rows + 1)[-1])
    annotations <- lapply(seq_along(plots), function(i) {
      list(x = ifelse(length(plots) == 1, 0.5, ifelse(i %% 2 == 0, 0.75, 0.25)),
           y = ifelse(i <= 2, 1, title_y[ceiling(i/2)]-0.075),
           text = subtitles[i],
           xref = "paper",
           yref = "paper",
           xanchor = "center",
           yanchor = "bottom",
           font = list(size = 22),
           showarrow = FALSE)})
    config(fig, toImageButtonOptions = list(format= 'svg', filename= 'plot',
                                            height= 800, width= 1200,
                                            scale= 1 )) %>%
      layout(annotations = annotations)
  })
  output$gridBox <- renderUI({
    data <- sort_data(data, input$gene_name)
    data <- data[order(data$Day), ]
    subtitles <- paste0("<i>", input$gene_name, "</i>")
    sorted_gene_id <- unique(data$gene_id)[order(match(unique(data$gene_id), input$gene_name))]
    data <- data[order(match(data$gene_id, sorted_gene_id)), ]
    plots <- c()
    plots <- lapply(unique(data$gene_id), function(gene) {
      gene_data <- data[data$gene_id == gene,]
      p_values <- get_pval(gene_data)
      gene_data <- gene_data %>%
        left_join(p_values, by = c("Day", "gene_id")) %>%
        select(-Day, -gene_id)      
      gene_data$maxExp <- 1.01 * max(gene_data$Expression)
      box <- plot_ly(gene_data, x = ~Day, y = ~Expression, color = ~Cell_Type, colors = colors(), 
                     type = 'box', legendgroup = ~Cell_Type, facet_col = ~gene_id,
                     showlegend = ifelse(gene == unique(data$gene_id)[1], TRUE, FALSE)) %>%
      add_trace(y = ~gene_data$maxExp, type = 'scatter',
                               mode = 'text', text = ~gene_data$label, textposition = 'middle',
                               textfont = list(color = '#000000', size = 10), showlegend = FALSE)
      box <- box %>%
        layout(xaxis = list(title = list(text ='Embryonic day')), 
               yaxis = list(title = list(text ='TPM')))
    })
    num_plots <- length(plots)
    num_cols <- 2  # Number of columns in the subplot grid
    
    # Calculate the number of rows dynamically
    num_rows <- ceiling(num_plots / num_cols)
    fig <- subplot(plots, nrows = num_rows, titleY= TRUE, titleX= TRUE, shareX = TRUE)
    #set y coordinate of titles
    title_y <- rev(seq(0, 1, length.out = num_rows + 1)[-1])
    annotations <- lapply(seq_along(plots), function(i) {
      list(x = ifelse(length(plots) == 1, 0.5, ifelse(i %% 2 == 0, 0.75, 0.25)),
           y = ifelse(i <= 2, 1, title_y[ceiling(i/2)]-0.05),
           text = subtitles[i],
           xref = "paper",
           yref = "paper",
           xanchor = "center",
           yanchor = "bottom",
           showarrow = FALSE)})
    
    config(fig, toImageButtonOptions = list(format= 'svg', filename= 'plot',
                                            height= 800, width= 1200,
                                            scale= 1 )) %>%
      layout(annotations = annotations)
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)

# library(rsconnect)
# deployApp()
