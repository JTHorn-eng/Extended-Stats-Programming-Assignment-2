




plot_seir = function(seir_results, title = 'SEIR Model Results') {

    par(mfcol=c(1,1),mar=c(4,4,2,1))

    n_days <- length(seir_results$S)
    days <- 1:n_days

    max_y <- max(c(seir_results$S, seir_results$E, seir_results$I), na.rm = TRUE)

    plot(
        days
        , seir_results$S
        , type = "l"
        , col = "black"
        ,ylim = c(0, max_y)
        , xlab = "Day"
        , ylab = "N"
    )
    lines(days, seir_results$E, col = "blue")
    lines(days, seir_results$I, col = "red")
    lines(days, seir_results$R, col = "green")


    legend("topright", legend = c("S (Susceptible)", "E (Exposed)", "I (Infectious)", "R (Recovered)"),
            col = c("black", "blue", "red", "green "), lty = 1, cex = 0.8)

    mtext(title, outer = TRUE, line = -1.5, cex = 1.2)

}

