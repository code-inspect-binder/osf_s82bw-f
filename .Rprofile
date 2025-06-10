if (interactive() && Sys.getenv("RSTUDIO") == "1") {
    message("Scheduling flowR addin installation...")
    later::later(function() {
        message("Installing flowR addin...")
        try(rstudioaddinflowr:::install_node_addin(), silent = FALSE)
    }, delay = 2)
    }
    