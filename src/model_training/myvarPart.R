myvarPart <-
  function(A, B, C = NA, AB, AC = NA, BC = NA, ABC = NA,
           model.type = NULL, A.name = "Factor A", B.name = "Factor B", 
           C.name = "Factor C", plot = TRUE, plot.digits = 3, cex.names = 1.5, 
           cex.values = 1.2, main = "", cex.main = 2, plot.unexpl = TRUE) {
    
    # version 1.7 (17 Mar 2017)
    
    if (!is.null(model.type)) message ("NOTE: Argument 'model.type' is no longer used.")
    
    partials <- c(A, B, C, AB, BC, AC, ABC)
    if (is.finite(partials[c(1:2, 4)]) && is.na(partials[c(3, 5:7)]))  {
      twofactors <- TRUE
    } else if (all(is.finite(partials))) {
      twofactors <- FALSE
    } else {
      stop ("You must provide numeric values for either A, B and AB (for variation partitioning among two factors) 
            or A, B, C, AB, BC, AC and ABC (for variation partitioning among three factors). See Details.")
    }
    
    if (!all(na.omit(partials) >= 0 & na.omit(partials) <= 1)) stop ("Values must be between 0 and 1.")
    
    totalexpl <- ifelse(twofactors, AB, ABC)
    unexpl <- 1 - totalexpl
    
    if (twofactors) {
      Apure <- totalexpl - B
      Bpure <- totalexpl - A
      ABoverlap <- totalexpl - Apure - Bpure
      output.names <- c(paste("Pure", A.name), paste("Pure", B.name),
                        paste0("Pure ", A.name, "_", B.name, " overlap"),
                        "Unexplained")
      results <- data.frame(c(Apure, Bpure, ABoverlap, unexpl),
                            row.names = output.names)
      
    } else { # end if 2 factors
      
      Apure <- totalexpl - BC
      Bpure <- totalexpl - AC
      Cpure <- totalexpl - AB
      ABoverlap <- totalexpl - Apure - Bpure - C
      BCoverlap <- totalexpl - Bpure - Cpure - A
      ACoverlap <- totalexpl - Apure - Cpure - B
      ABCoverlap <- totalexpl - Apure - Bpure - Cpure - ABoverlap - BCoverlap - ACoverlap
      output.names <- c(paste("Pure", A.name),
                        paste("Pure", B.name),
                        paste("Pure", C.name),
                        paste0("Pure ", A.name, "_", B.name, " overlap"),
                        paste0("Pure ", B.name, "_", C.name, " overlap"),
                        paste0("Pure ", A.name,"_", C.name," overlap"),
                        paste0(A.name,"_",B.name,"_",C.name," overlap"),
                        "Unexplained")
      results <- data.frame(c(Apure, Bpure, Cpure, ABoverlap, BCoverlap,
                              ACoverlap, ABCoverlap, unexpl),
                            row.names = output.names)
    }  # end else
    
    colnames(results) <- "Proportion"
    #n <- nrow(results)
    #if(model.type == "GLM")  results <- results[1:(n-1),]  # deletes "unexplained" line (data unavailable for GLM)
    
    if (plot) {  # adapted from Daniel's http://stackoverflow.com/questions/1428946/venn-diagrams-with-r
      
      circle <- function(x, y, r) {
        ang <- seq(0, 2*pi, length = 100)
        xx <- x + r * cos(ang)
        yy <- y + r * sin(ang)
        polygon(xx, yy)
      }  # end circle funtion (by Daniel)
      
      Apure <- round(Apure, plot.digits)  # shorten values for plotting
      Bpure <- round(Bpure, plot.digits)
      ABoverlap <- round(ABoverlap, plot.digits)
      if(!twofactors) {
        Cpure <- round(Cpure, plot.digits)
        BCoverlap <- round(BCoverlap, plot.digits)
        ACoverlap <- round(ACoverlap, plot.digits)
        ABCoverlap <- round(ABCoverlap, plot.digits)
      }
      
      if (twofactors) {
        plot(0, 0, ylim = c(-1, 10), xlim = c(-1, 7), type = "n", axes = FALSE,
             ylab = "", xlab = "", main = main, cex.main = cex.main)
        circle(3,3,3)
        circle(3,6,3)
        text(x = c(3, 3), y = c(9.5, -0.5), labels = c(A.name, B.name),
             cex = cex.names)
        text(x = c(3, 3, 3), y = c(7, 4.75, 2), c(Apure, ABoverlap, Bpure),
             cex = cex.values)
        
      } else { # end if 2 factors
        
        plot(0, 0, ylim = c(-1, 10), xlim = c(-1, 10), type = "n", axes = FALSE,
             ylab = "", xlab = "", main = main, cex.main = cex.main)
        circle(3, 6, 3)
        circle(6, 6, 3)
        circle(4.5, 3,  3)
        #Cname.loc = ifelse((plot.unexpl), 6, 4.5)
        text(x = c(2.5, 6.5, 4.5), y = c(9.5, 9.5, -0.5),
             labels = c(A.name, B.name, C.name), cex = cex.names, adj = c(0.5, 0.5, 0))
        text(x = c(1.8, 7.2, 4.5, 4.5, 2.8, 6.2, 4.5), y = c(6.6, 6.6, 2, 7, 4, 4, 5), labels = c(Apure, Bpure, Cpure, ABoverlap, ACoverlap, BCoverlap, ABCoverlap), cex = cex.values)
      } # end if 2 factors else
      
      if (plot.unexpl)  {
        rect(-1, -1, 10, 10)
        text(x = -0.9, y = -0.2, label = paste0("Unexplained\n", round(unexpl, plot.digits)), adj = 0, cex = cex.values)
      }
      
    }  # end if plot
    
    if (all.equal(sum(results, na.rm = TRUE), 1)) cat("")
    else warning ("Results don't sum up to 1; are you sure your input data are correct?")  # but this doesn't work because results always sum to 1 anyway
    
    results
    
  }