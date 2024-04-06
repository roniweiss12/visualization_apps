library(shiny)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggsignif)
library(coin)
library(RColorBrewer)
library(colorspace)
library(bslib)
library(gridExtra)
library(highcharter)

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

TPM_FILE <- "TPM_allGenes/230606_ALL_rsem.merged.gene_tpm.csv"

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
ui <- fluidPage(
  
  titlePanel("RNA Seq Results"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Show gene specific expression (in TPM)."),
      
      # Allow multiple gene selections
      selectizeInput("gene_name", 
                     label = "Choose gene name(s):",
                     choices = ord_mtx$gene_id,
                     selected = c("Amh"),
                     multiple = TRUE)
    ),
    
    mainPanel(highchartOutput('scatterPlot'), htmlOutput('boxPlot'))
  )
)


# Define server logic ----
server <- function(input, output) {
  # Define export options
  export <- list(
    list(
      text = "PNG",
      onclick = JS("function () {
                   this.exportChart({ type: 'image/png' }); }")
    ),
    list(
      text = "JPEG",
      onclick = JS("function () {
                   this.exportChart({ type: 'image/jpeg' }); }")
    ),
    list(
      text = "SVG",
      onclick = JS("function () {
                   this.exportChart({ type: 'image/svg+xml' }); }")
    ),
    list(
      text = "PDF",
      onclick = JS("function () {
                   this.exportChart({ type: 'application/pdf' }); }")
    )
  )
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
      stringsAsFactors = FALSE
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
    return(p_df)
  }
  output$boxPlot <- renderUI({
    plots <- c()
    for (gene in input$gene_name) {
      data <- sort_data(data, gene)
      
      gene_id_list <- unique(data$gene_id)
      
      myboxplotData <- lapply(gene_id_list, function(gene) {
        subset_data <- subset(data, gene_id == gene)
        data_to_boxplot(subset_data, Expression, Cell_Type, Day) %>%
          mutate(id = gene)  # Add Gene_ID as a column
      })
      myboxplotData <- bind_rows(myboxplotData)
      myboxplotData$name <- paste(myboxplotData$id, myboxplotData$name, sep = " - ")
      myboxplotData$name <- gsub("(.*?) - ", "<i>\\1</i> - ", myboxplotData$name)
      combined_colors <- get_color_pallete(data$gene_id)
      p_values <- get_pval(data)
      alpha_levels <- c(0.05, 0.01)
      p_values$label <- ifelse(p_values$p_value >= alpha_levels[1], "-ns-",
                               ifelse(p_values$p_value >= alpha_levels[2], "-*-", "-**-"))
      data <- data[order(data$Day), ]
      unique_days <- unique(data$Day)
      maxExpression <- max(data$Expression)
      hc <- highchart() %>%
        hc_colors(unlist(combined_colors)) %>%
        hc_xAxis(categories = unique_days, title = list(text = "Day"))%>%
        hc_add_series_list(myboxplotData)%>%
        hc_xAxis(title = list(text = "Embryonic day"))%>%
        hc_yAxis(title = list(text = "Expression (TPM)"))%>%
        hc_title(text= paste("<i>", unique(data$gene_id), "</i>", "gonadal expression")) %>%
        hc_legend(enabled= TRUE)%>%
        hc_plotOptions(area = list(relativeXValue = TRUE))%>%
        hc_exporting(
          enabled = TRUE,
          formAttributes = list(target = "_blank"),
          buttons = list(contextButton = list(
            text = "Export",
            theme = list(fill = "transparent"),
            menuItems = export))) %>%
        hc_chart(backgroundColor = "#ffffff")%>%
        hc_annotations(
          list(labels = list(
              list(point = list(x = 0, y = maxExpression*1.1, xAxis = 0, yAxis = 0), text = as.character(p_values$label[1])),
              list(point = list(x = 1, y = maxExpression*1.1, xAxis = 0, yAxis = 0), text = as.character(p_values$label[2])),
              list(point = list(x = 2, y = maxExpression*1.1, xAxis = 0, yAxis = 0), text = as.character(p_values$label[3])),
              list(point = list(x = 3, y = maxExpression*1.1, xAxis = 0, yAxis = 0), text = as.character(p_values$label[4]))
            ), labelOptions = list(backgroundColor = "white", borderColor = "white")))
      plots <- c(plots, list(hc))
    }
    hw_grid(plots)
  })
  output$scatterPlot <- renderHighchart({
    if (is.null(input$gene_name) || input$gene_name == "") {
      return(tags$p("Please choose genes"))
    }
    data <- sort_data(data, input$gene_name)
    combined_colors <- get_color_pallete(data$gene_id)
    # Combine the Cell_Type and gene_id columns to create a new column for grouping
    data$group <- interaction(data$gene_id, data$Cell_Type)
    data$group <- gsub("\\.", " - ", data$group)
    data <- data[order(data$Day), ]
    unique_days <- as.numeric(unique(data$Day))
    genes <- paste0(unique(data$gene_id), collapse = ', ')
    plot_title <- paste("<i>", genes, "</i>", "gonadal expression")

    # Plot scatter and line
    highchart() %>%
      hc_colors(unlist(combined_colors)) %>%
      hc_add_series(data, "line", 
                    hcaes(x = as.numeric(Day), y = Expression_mean, 
                          group = group),
                    marker = list(enabled = FALSE)) %>%
      hc_add_series(data, "scatter", 
                    hcaes(x = as.numeric(Day), y = Expression, 
                          group = group), 
                    dataLabels = list(enabled = FALSE),
                    showInLegend = FALSE) %>%
      hc_yAxis(title = list(text = "Expression (TPM)")) %>%
      hc_xAxis(tickPositions = unique_days, title = list(text = "Embryonic day")) %>%
      hc_legend(enabled = TRUE) %>%
      hc_tooltip(pointFormat = "Day: {point.x}<br/>Expression: {point.y}<br/>") %>%
      hc_title(text = plot_title) %>%
      hc_chart(backgroundColor = "#ffffff")%>%
      hc_exporting(
        enabled = TRUE,
        formAttributes = list(target = "_blank"),
        buttons = list(contextButton = list(
          text = "Export",
          theme = list(fill = "transparent"),
          menuItems = export
        ))
      ) 
  })
}
# Run the app ----
shinyApp(ui = ui, server = server)
