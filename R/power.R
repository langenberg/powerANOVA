#' @export
#' @importFrom matrixcalc is.positive.definite
get_sums_of_squares <- function(
    n = 10,
    mu = c(1, 1, 2, 2.1),
    Sigma,
    var = 1,
    cov = 0.5,
    effect,
    idata,
    ...
) {

    if (missing(Sigma)) {
        Sigma <- matrix(cov, ncol = length(mu), nrow = length(mu))
        diag(Sigma) <- var
    }

    if (is.vector(Sigma)) {
        Sigma <- matrix(Sigma, ncol = 1)
    }

    if (!isSymmetric(Sigma) || !is.positive.definite(Sigma)) {
        stop("Sigma is not positive definite.")
    }

    if (is.vector(mu)) {
        mu <- matrix(mu, ncol = 1)
    }

    if (nrow(mu) == 1L) {
        return(list(SSH = mu[1,1,drop=F]^2, SSE = Sigma[1,1,drop=F]))
    }

    if (missing(idata)) {
        idata <- expand.grid(A = paste0("A", 1:length(mu)))
        effect <- "A"
    }

    idesign <- as.formula(paste0("~", paste0(names(idata), collapse = "*")))
    contrasts_arg <- as.list(rep_len("contr.poly", ncol(idata)))
    names(contrasts_arg) <- names(idata)
    B_matrix <- model.matrix(idesign, idata, contrasts.arg = contrasts_arg)
    C_matrix <- solve(B_matrix)
    term_index <- which(effect == attr(terms(idesign), "term.labels"))
    term_rows <- which(attr(B_matrix, "assign") == term_index)
    C_effect <- C_matrix[term_rows, , drop = F]

    SSH <- C_effect %*% mu
    SSH <- SSH %*% t(SSH)
    SSE <- C_effect %*% Sigma %*% t(C_effect)

    list(SSH = SSH, SSE = SSE)
}

#' @export
power_cohens_d <- function(n = 10, cohens_d, alpha = 0.05, ...) {
    t <- abs(cohens_d * sqrt(n))
    pt(qt(alpha/2, df = n-1), ncp = t, df = n-1) + (1 - pt(qt(1-alpha/2, df = n-1), ncp = t, df = n-1))
}

#' @export
#' @importFrom matrixcalc is.positive.definite
power_mu_cov_contrast <- function(n = 10, mu, Sigma, contrast, ...) {
    if (!isSymmetric(Sigma) || !is.positive.definite(Sigma)) {
        stop("Sigma is not positive definite.")
    }
    mu <- contrast %*% mu
    Sigma <- contrast %*% Sigma %*% t(contrast)
    power_mu_cov(n = 10, mu = mu, Sigma = Sigma)
}

#' @export
power_mu_cov <- function(n = 10, ...) {

    sums <- get_sums_of_squares(n, ...)

    if (ncol(sums[[1]]) > 1L) {
        power_mu_cov_multi(n, ...)
    } else {
        power_mu_cov_uni(n, ...)
    }
}

#' @keywords internal
power_mu_cov_uni <- function(n = 10, ...) {
    sums <- get_sums_of_squares(n, ...)
    d_mean <- sqrt(sum(diag(sums[[1]])))
    d_sd <- sqrt(sum(diag(sums[[2]])))
    cohens_d <- d_mean / d_sd


    power_cohens_d(n = n, cohens_d = cohens_d, ...)
}

#' @keywords internal
power_mu_cov_multi <- function(n, ...) {
    sums <- get_sums_of_squares(n, ...)
    p_eta_sq <- sum(diag(sums[[1]])) / (sum(diag(sums[[1]])) + sum(diag(sums[[2]])))
    df1 <- ncol(sums[[1]])


    power_petasq(n, p_eta_sq, df1, ...)
}

#' @export
power_petasq <- function(n, p_eta_sq, df1 = 1, alpha = 0.05, ...) {
    df2 <- df1 * n - df1
    # F_value <- p_eta_sq*df2/(df1 - p_eta_sq*df1)
    ncp <- convert_petasq_f2(p_eta_sq)*df1*n
    1 - pf(qf(1-alpha, df1 = df1, df2 = df2), df1 = df1, df2 = df2, ncp = ncp)
}

#' @export
power_gui <- function(browser = FALSE) {
    # Run the application
    if (browser) {
        shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
    } else {
        shinyApp(ui = ui, server = server)
    }
}

#' @export
convert_cohens_d_petasq <- function(cohens_d, n, population = TRUE) {
    if (population) {
        cohens_d^2/(1 + cohens_d^2)
    } else {
        df2 <- n-1
        cohens_d^2*n/(df2 + cohens_d^2*n)
    }
}

#' @export
convert_petasq_f2 <- function(p_eta_sq) {
    p_eta_sq / (1 - p_eta_sq)
}

#' @export
convert_F_petasq <- function(Fval, df1, df2) {
    (Fval * df1)/(Fval * df1 + df2)
}

#' @export
convert_petasq_F <- function(p_eta_sq, df1, df2) {
    p_eta_sq * df2 / (df1 - p_eta_sq * df1)
}

#' @export
convert_petasq_cohens_d <- function(p_eta_sq, n, population = TRUE) {
    if (population) {
        sqrt(p_eta_sq / (1 - p_eta_sq))
    } else {
        df2 <- n - 1
        sqrt(p_eta_sq * df2 / (n - p_eta_sq * n))
    }
}

#' @export
power_plot_cohens_d <- function(n = 10, cohens_d, alpha = 0.05) {
    t <- cohens_d * sqrt(n)

    lower_null <- qt(alpha/2, df = n - 1) - 1
    upper_null <- qt(1-alpha/2, df = n - 1) + 1
    lower_alt <- qt(alpha/2, df = n - 1, ncp = t) - 1
    upper_alt <- qt(1-alpha/2, df = n - 1, ncp = t) + 1

    mylimits <- c(min(c(lower_null,lower_alt)), max(c(upper_null,upper_alt)))

    ggplot() +
        geom_function(
            fun = function(x) dt(x, n - 1),
            mapping = aes(color = factor(
                "null hypothesis",
                levels = c("null hypothesis", "population")
            ))
        ) +
        geom_function(
            fun = function(x) dt(x, n - 1, ncp = t),
            mapping = aes(color = factor(
                "population", levels = c("null hypothesis", "population")
            ))
        ) +
        stat_function(
            fun = function(x) dt(x, n - 1, ncp = t),
            geom = "area",
            alpha = 0.2,
            fill = "steelblue",
            xlim = c(qt(1-alpha/2, df = n - 1), mylimits[2])
        ) +
        stat_function(
            fun = function(x) dt(x, n - 1, ncp = t),
            geom = "area",
            alpha = 0.2,
            fill = "steelblue",
            xlim = c(mylimits[1], qt(alpha/2, df = n - 1))
        ) +
        geom_vline(xintercept = qt(1-alpha/2, df = n - 1)) +
        geom_vline(xintercept = qt(alpha/2, df = n - 1)) +
        # xlim(c(qt(alpha/2, df = n - 1) - 1, qt(1-alpha/2, df = n - 1) + 1)) +
        xlim(mylimits) +
        ylab("probability density") +
        xlab("t") +
        theme_bw() +
        theme(legend.position = "bottom") +
        guides(color = guide_legend(title = "Hypothesis"))
}

#' @export
power_plot_p_eta_sq <- function(n, p_eta_sq, df1 = 1, alpha = 0.05, ...) {
    df2 <- df1 * n - df1
    ncp <- convert_petasq_f2(p_eta_sq)*df1*n

    mylimits <- c(0, qf(0.95, df1, df2, ncp) + 1)

    ggplot() +
        geom_function(
            fun = function(x)
                df(x, df1 = df1, df2 = df2),
            mapping = aes(color = factor(
                "null hypothesis",
                levels = c("null hypothesis", "population")
            ))
        ) +
        geom_function(
            fun = function(x)
                df(
                    x,
                    df1 = df1,
                    df2 = df2,
                    ncp = ncp
                ),
            mapping = aes(color = factor(
                "population", levels = c("null hypothesis", "population")
            ))
        ) +
        stat_function(
            fun = function(x)
                df(
                    x,
                    df1 = df1,
                    df2 = df2,
                    ncp = ncp
                ),
            geom = "area",
            alpha = 0.2,
            fill = "steelblue",
            xlim = c(qf(1-alpha, df1 = df1, df2 = df2), mylimits[2])
        ) +
        geom_vline(xintercept = qf(1-alpha, df1 = df1, df2 = df2)) +
        # xlim(c(0, qf(1-alpha, df1 = df1, df2 = df2) + 1)) +
        xlim(mylimits) +
        ylab("probability density") +
        xlab("F") +
        theme_bw() +
        theme(legend.position = "bottom") +
        guides(color = guide_legend(title = "Hypothesis"))
}

#' @export
power_plot_mu_cov_contrast <- function(n = 10, mu, Sigma, contrast, ...) {
    mu <- contrast %*% mu
    Sigma <- contrast %*% Sigma %*% t(contrast)
    power_plot_mu_cov(n = 10, mu = mu, Sigma = Sigma, ...)
}

#' @export
power_plot_mu_cov <- function(n = 10, ...) {
    sums <- get_sums_of_squares(n, ...)

    if (ncol(sums[[1]]) == 1L) {
        d_mean <- sqrt(sums[[1]][1,1])
        d_sd <- sqrt(sums[[2]][1,1])
        cohens_d <- d_mean / d_sd
        power_plot_cohens_d(n = n, cohens_d = cohens_d)
    } else {
        p_eta_sq <- sum(diag(sums[[1]])) / (sum(diag(sums[[1]])) + sum(diag(sums[[2]])))
        power_plot_p_eta_sq(n = n, p_eta_sq = p_eta_sq, df1 = ncol(sums[[1]]), ...)
    }
}

#' @keywords internal
power_plot_mu_cov_uni <- function(n = 10, ...) {
    sums <- get_sums_of_squares(n, ...)
    d_mean <- sqrt(sums[[1]][1,1])
    d_sd <- sqrt(sums[[2]][1,1])
    cohens_d <- d_mean / d_sd
}
