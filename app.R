# Load required libraries
library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(arrow)
library(DBI)
library(duckdb)

# # Create directories if they don't exist
# if (!dir.exists("data")) {
#     dir.create("data")
# }
# if (!dir.exists("parquet_data")) {
#     dir.create("parquet_data")
# }



#################### UI ####################

ui <- page_sidebar(
    title = "Mass Spectrometry Data Analysis",
    sidebar = sidebar(
        selectInput("file_selection", "Select Dataset:",
                    choices = c("Indoor Dust" = "./data/IndoorDust_SusFut_2025_GroupA_A1_cropped.parquet",
                                "Alkanes C9-C20 + Benzophenone@965 sec" = "./data/alkanestd_C7_C40_bensophenone_1_cropped.parquet"),
                    selected = "./data/alkanestd_C7_C40_bensophenone_1_cropped.parquet"),
        hr(),
        actionButton("add_range", "Add m/z Range", class = "btn-primary"),
        uiOutput("mz_inputs"),
        actionButton("openEIC", "Extract Ion Chromatogram", class = "btn-success")
    ),
    
    # Scrollable container with cards sized to fit screen for first two
    div(
        style = "height: 100vh; overflow-y: auto; padding-right: 10px;",
        # First two cards in a row to share screen space
        layout_column_wrap(
            width = 1/2, # Each takes half the width
            heights_equal = "row",
            card(
                card_header("Total Ion Chromatogram (TIC)"),
                plotlyOutput("ms_plot", height = "80vh", width = "100%"), # Use viewport height
                fill = TRUE
            ),
            card(
                card_header("Mass Spectrum at Selected Retention Time"),
                plotlyOutput("spectrum_plot", height = "80vh", width = "100%"), # Use viewport height
                fill = TRUE
            )
        ),
        # Third card gets full width and extends below
        card(
            card_header("File Information"),
            verbatimTextOutput("file_info"),
            fill = TRUE,
            style = "margin-top: 10px; min-height: 300px;" # Add some spacing and minimum height
        )
    )
)

#################### SERVER ####################

server <- function(input, output, session) {
    
    mz_range_count <- reactiveVal(1)
    
    observeEvent(input$add_range, {
        mz_range_count(mz_range_count() + 1)
    })
    
    output$mz_inputs <- renderUI({
        n <- mz_range_count()
        lapply(1:n, function(i) {
            fluidRow(
                column(6,
                       numericInput(paste0("min_mass_", i), paste("Min m/z (", i, "):"), value = 100 + (i - 1) * 100, min = 0)
                ),
                column(6,
                       numericInput(paste0("max_mass_", i), paste("Max m/z (", i, "):"), value = 200 + (i - 1) * 100, min = 0)
                )
            )
        })
    })
    
    db_conn <- reactiveVal(NULL)
    # Removed default_parquet_file - now using input$file_selection directly
    
    onSessionEnded(function() {
        if (!is.null(db_conn())) {
            dbDisconnect(db_conn())
        }
    })
    
    # Initial connection - now uses the selected file from input
    observe({
        req(input$file_selection)
        req(file.exists(input$file_selection))
        con <- setup_duckdb(input$file_selection)
        db_conn(con)
        
        output$file_info <- renderText({
            req(db_conn())
            get_file_info_duckdb(db_conn())
        })
        
        output$ms_plot <- renderPlotly({
            req(db_conn())
            plot_tic_duckdb(db_conn())
        })
        
        output$eic_plot <- renderPlotly({ NULL })  # Clear EIC initially
    })
    
    # File selection observer - now properly responds to changes
    observeEvent(input$file_selection, {
        file_path <- input$file_selection
        req(file.exists(file_path))
        
        if (!is.null(db_conn())) dbDisconnect(db_conn())
        con <- setup_duckdb(file_path)
        db_conn(con)
        
        output$file_info <- renderText({
            req(db_conn())
            get_file_info_duckdb(db_conn())
        })
        
        output$ms_plot <- renderPlotly({
            req(db_conn())
            plot_tic_duckdb(db_conn())
        })
        
        output$eic_plot <- renderPlotly({ NULL })  # Clear EIC on file change
    })
    
    # Show EIC in separate plot when button is clicked
    observeEvent(input$openEIC, {
        req(db_conn())
        
        # Get TIC data
        tic_query <- "
    SELECT rtime, SUM(intensity) as total_intensity
    FROM ms_data
    GROUP BY rtime
    ORDER BY rtime
  "
        tic_data <- dbGetQuery(db_conn(), tic_query)
        
        # Loop through all m/z ranges
        n <- mz_range_count()
        eic_traces <- list()
        
        for (i in 1:n) {
            min_mz <- input[[paste0("min_mass_", i)]]
            max_mz <- input[[paste0("max_mass_", i)]]
            if (!is.null(min_mz) && !is.null(max_mz) && min_mz <= max_mz) {
                query <- sprintf("
        SELECT rtime, SUM(intensity) as filtered_intensity
        FROM ms_data
        WHERE mz BETWEEN %f AND %f
        GROUP BY rtime
        ORDER BY rtime
      ", min_mz, max_mz)
                eic_data <- dbGetQuery(db_conn(), query)
                eic_traces[[i]] <- list(data = eic_data, label = paste0("EIC (", min_mz, "-", max_mz, " m/z)"))
            }
        }
        
        # Render plot
        output$ms_plot <- renderPlotly({
            p <- plot_ly(source = "main_plot") %>%
                add_trace(
                    data = tic_data,
                    x = ~rtime,
                    y = ~total_intensity,
                    type = 'scatter',
                    mode = 'lines',
                    name = "TIC",
                    line = list(color = 'royalblue')
                )
            
            colors <- c("forestgreen", "orange", "purple", "brown", "darkred", "darkblue")
            for (i in seq_along(eic_traces)) {
                eic_data <- eic_traces[[i]]$data
                label <- eic_traces[[i]]$label
                color <- colors[(i %% length(colors)) + 1]
                
                p <- p %>% add_trace(
                    x = eic_data$rtime,
                    y = eic_data$filtered_intensity,
                    type = 'scatter',
                    mode = 'lines',
                    name = label,
                    line = list(color = color)
                )
            }
            
            p %>% layout(
                title = "TIC + EICs",
                xaxis = list(title = "Retention Time (seconds)"),
                yaxis = list(title = "Intensity"),
                font = list(size = 14),
                margin = list(l = 60, r = 40, b = 60, t = 60, pad = 10),
                hovermode = "closest"
            )
        })
    })
    
    # Show current file name - now uses input$file_selection
    output$current_file <- renderText({
        paste("Loaded file:", basename(input$file_selection))
    })
    
    # Mass spectrum on click
    observeEvent(event_data("plotly_click", source = "main_plot"), {
        click_data <- event_data("plotly_click", source = "main_plot")
        req(click_data)
        clicked_rtime <- click_data$x
        
        query <- sprintf("
      SELECT spectrum_id, rtime
      FROM ms_data
      ORDER BY ABS(rtime - %f)
      LIMIT 1
    ", clicked_rtime)
        
        closest <- dbGetQuery(db_conn(), query)
        req(nrow(closest) > 0)
        spectrum_id <- closest$spectrum_id[1]
        
        spectrum_query <- sprintf("
      SELECT mz, intensity
      FROM ms_data
      WHERE spectrum_id = %d
      ORDER BY mz
    ", spectrum_id)
        
        spectrum_data <- dbGetQuery(db_conn(), spectrum_query)
        spectrum_data$norm_intensity <- 100 * spectrum_data$intensity / max(spectrum_data$intensity, na.rm = TRUE)
        
        output$spectrum_plot <- renderPlotly({
            req(nrow(spectrum_data) > 0)
            
            plot_ly(
                data = spectrum_data,
                x = ~mz,
                y = ~norm_intensity,
                type = 'bar',
                marker = list(color = 'red'),
                text = ~paste('m/z: ', round(mz, 2), 
                              '<br>Intensity: ', round(intensity, 0),
                              '<br>Relative Intensity: ', round(norm_intensity, 1), '%'),
                hoverinfo = 'text',
                textposition = 'none'  # This hides text labels on the bars
            ) %>%
                layout(
                    title = paste("Spectrum ID:", spectrum_id, ", Retention Time: ", round(closest$rtime, 2)),
                    xaxis = list(title = "m/z"),
                    yaxis = list(title = "Relative Intensity (%)", range = c(0, 105)),
                    autosize = TRUE,
                    margin = list(l = 60, r = 40, b = 60, t = 60, pad = 10),
                    font = list(size = 14),
                    hovermode = "closest"
                ) %>%
                config(responsive = TRUE)
        })
    })
}

shinyApp(ui, server)
