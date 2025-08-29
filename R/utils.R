#### Functions ####

# Function to setup DuckDB connection with the parquet file
setup_duckdb <- function(parquet_file) {
    tryCatch({
        con <- dbConnect(duckdb::duckdb())
        dbExecute(con, sprintf("CREATE VIEW ms_data AS SELECT * FROM parquet_scan('%s')", parquet_file))
        return(con)
    }, error = function(e) {
        showNotification(paste("Error connecting to database:", e$message), type = "error")
        return(NULL)
    })
}

# Function to plot TIC with plotly using DuckDB
plot_tic_duckdb <- function(con) {
    withProgress(message = 'Generating TIC...', {
        # Query to aggregate intensity by retention time for TIC
        tic_query <- "
            SELECT rtime, SUM(intensity) as total_intensity
            FROM ms_data
            GROUP BY rtime
            ORDER BY rtime
        "
        
        # Execute query
        tic_data <- dbGetQuery(con, tic_query)
        
        # Create plotly figure
        fig <- plot_ly(
            data = tic_data, # or eic_data
            x = ~rtime,
            y = ~total_intensity, # or ~filtered_intensity
            type = 'scatter',
            mode = 'lines',
            source = "main_plot", # <-- Add this line
            line = list(color = 'royalblue'),
            hoverinfo = 'text',
            text = ~paste('RT: ', round(rtime, 2), ' sec<br>Intensity: ', formatC(total_intensity, format = "e", digits = 2))
        ) %>%
            layout(
                title = "Total Ion Chromatogram (TIC)",
                xaxis = list(title = "Retention Time (seconds)"),
                yaxis = list(title = "Intensity"),
                font = list(size = 14),
                margin = list(l = 60, r = 40, b = 60, t = 60, pad = 10),
                hovermode = "closest"
            )
        
        return(fig)
    })
}

# Function to plot EIC with plotly using DuckDB
plot_eic_duckdb <- function(con, min_mass, max_mass) {
    withProgress(message = 'Generating EIC...', {
        # Query to filter by m/z range and aggregate
        eic_query <- sprintf("
            SELECT rtime, SUM(intensity) as filtered_intensity
            FROM ms_data
            WHERE mz BETWEEN %f AND %f
            GROUP BY rtime
            ORDER BY rtime
        ", min_mass, max_mass)
        
        # Execute query
        eic_data <- dbGetQuery(con, eic_query)
        
        # Create plotly figure
        
        fig <- plot_ly(
            data = eic_data,
            x = ~rtime,
            y = ~filtered_intensity, # or ~filtered_intensity
            type = 'scatter',
            mode = 'lines',
            source = "main_plot", # <-- Add this line
            line = list(color = 'forestgreen'),
            hoverinfo = 'text',
            text = ~paste('RT: ', round(rtime, 2), ' sec<br>Intensity: ', formatC(filtered_intensity, format = "e", digits = 2))
        ) %>%
            layout(
                title = paste0("Extracted Ion Chromatogram (EIC): m/z ", min_mass, " - ", max_mass),
                xaxis = list(title = "Retention Time (seconds)"),
                yaxis = list(title = "Intensity"),
                font = list(size = 14),
                margin = list(l = 60, r = 40, b = 60, t = 60, pad = 10),
                hovermode = "closest"
            )
        
        return(fig)
    })
}

# Function to get file information using DuckDB
get_file_info_duckdb <- function(con) {
    tryCatch({
        # Total number of spectra
        num_spectra_query <- "SELECT COUNT(DISTINCT spectrum_id) FROM ms_data"
        num_spectra <- dbGetQuery(con, num_spectra_query)[[1]]
        
        # Retention time range
        rt_range_query <- "SELECT MIN(rtime), MAX(rtime) FROM ms_data"
        rt_range <- dbGetQuery(con, rt_range_query)
        
        # MS levels present
        ms_levels_query <- "SELECT DISTINCT ms_level FROM ms_data ORDER BY ms_level"
        ms_levels <- dbGetQuery(con, ms_levels_query)[[1]]
        
        # Spectra counts by MS level
        ms_level_counts_query <- "SELECT ms_level, COUNT(DISTINCT spectrum_id) FROM ms_data GROUP BY ms_level ORDER BY ms_level"
        ms_level_counts <- dbGetQuery(con, ms_level_counts_query)
        
        # m/z range
        mz_range_query <- "SELECT MIN(mz), MAX(mz) FROM ms_data"
        mz_range <- dbGetQuery(con, mz_range_query)
        
        # Combine all information
        info_text <- paste0(
            "File successfully loaded!\n\n",
            "Number of spectra: ", num_spectra, "\n",
            "Retention time range: ", paste(round(as.numeric(rt_range), 2), collapse = " - "), " seconds\n",
            "MS levels present: ", paste(ms_levels, collapse = ", "), "\n\n",
            "Spectra counts by MS level:\n"
        )
        
        # Add MS level counts
        for (i in 1:nrow(ms_level_counts)) {
            info_text <- paste0(
                info_text,
                "  MS", ms_level_counts$ms_level[i], ": ", ms_level_counts[[2]][i], "\n"
            )
        }
        
        # Add m/z range
        info_text <- paste0(
            info_text,
            "\nApproximate m/z range: ",
            paste(round(as.numeric(mz_range), 2), collapse = " - "), "\n"
        )
        
        return(info_text)
    }, error = function(e) {
        return(paste("Error getting file information:", e$message))
    })
}