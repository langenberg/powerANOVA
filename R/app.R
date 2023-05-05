#' @keywords internal
parse_design <- function(input) {
    design <- gsub(" ", "", input)

    if (!grepl("^[\\d]+(x[\\d]+)*$", design, perl = T)) {
        design <- list(valid = FALSE, error = "Invalid design.")
    } else {
        design <- strsplit(design, "x")[[1]]
        design <- as.integer(design)
        if (all(design < 2)) {
            design <- list(valid = FALSE, error = "At least one factor must have more than one level.")
        } else {
            idata <- lapply(1:length(design), function(x) paste0(LETTERS[x], 1:design[x]))
            names(idata) <- LETTERS[1:length(design)]
            idata <- do.call(expand.grid, idata)

            idesign <- as.formula(paste0("~", paste0(names(idata), collapse = "*")))
            term_labels <- attr(terms(idesign), "term.labels")
            term_labels <- paste0("[", paste0(term_labels, collapse = ", "), "]")

            order_idata <- paste0(
                "[",
                paste0(apply(idata, 1, paste0, collapse = ""), collapse = ", "),
                "]"
            )

            design <-
                list(
                    design = design,
                    idesign = idesign,
                    idata = idata,
                    valid = TRUE
                )
        }
    }
    design
}

#' @keywords internal
parse_bool <- function(input) {
    bool <- gsub(" ", "", input)
    bool <- gsub("[", "", bool, fixed = T)
    bool <- gsub("]", "", bool, fixed = T)

    if (!is.logical(input)) {
        if (bool %in% c("TRUE", "T", "FALSE", "F")) {
            list(value = bool %in% c("TRUE", "T"), valid = TRUE)
        } else {
            list(valid = FALSE, error = "Not a boolean.")
        }
    } else {
        list(value = bool, valid = TRUE)
    }
}

#' @keywords internal
parse_int <- function(input, n = 1) {
    int <- gsub(" ", "", input)
    int <- gsub("[", "", int, fixed = T)
    int <- gsub("]", "", int, fixed = T)

    if (!grepl("^[\\d]+(,[\\d]+)*$", int, perl = T)) {
        int <- list(valid = FALSE, error = "Not an integer.")
    } else {
        int <- strsplit(int, ",")[[1]]
        if (is.null(n) || length(int) == n) {
            int <- list(value = as.integer(int), valid = TRUE)
        } else {
            int <- list(valid = FALSE, error = "Wrong number of integers.")
        }
    }
    int
}

#' @keywords internal
parse_dbl <- function(input, n = NULL) {
    dbl <- gsub(" ", "", input)
    dbl <- gsub("[", "", dbl, fixed = T)
    dbl <- gsub("]", "", dbl, fixed = T)

    if (!grepl("^[\\d]+(\\.[\\d]+){0,1}(,[\\d]+(\\.[\\d]+){0,1})*$", dbl, perl = T)) {
        dbl <- list(valid = FALSE, error = "Not a double.")
    } else {
        dbl <- strsplit(dbl, ",")[[1]]
        if (is.null(n) || length(dbl) == n) {
            dbl <- list(value = as.numeric(dbl), valid = TRUE)
        } else {
            dbl <- list(valid = FALSE, error = "Wrong number of doubles.")
        }
    }
    dbl
}

#' @keywords internal
#' @importFrom matrixcalc is.positive.definite
parse_dbl_matrix <- function(input, dim = NULL, ...) {
    Sigma <- gsub(" ", "", input)
    Sigma <- gsub("[", "", Sigma, fixed = T)
    Sigma <- gsub("]", "", Sigma, fixed = T)

    if (!grepl("^[\\d]+(\\.[\\d]+){0,1}(,[\\d]+(\\.[\\d]+){0,1})*$", Sigma, perl = T)) {
        Sigma <- list(valid = FALSE, error = "Not a double.")
    } else {
        Sigma <- strsplit(Sigma, ",")[[1]]
        if (is.null(dim) || length(Sigma) == dim^2) {
            if (is.null(dim)) {
                dim <- 1
            }
            Sigma <- matrix(as.numeric(Sigma), ncol = dim, byrow = TRUE, ...)

            if (!isSymmetric(Sigma) || !is.positive.definite(Sigma)) {
                Sigma <- list(valid = FALSE, error = "Sigma is not positive definite.")
            } else {
                Sigma <- list(value = Sigma, valid = TRUE)
            }

        } else {
            Sigma <- list(valid = FALSE, error = "Wrong number of doubles.")
        }
    }
    Sigma
}

#' @keywords internal
parse_effect <- function(input, design) {
    effect <- gsub(" ", "", input)

    if (!design$valid) {
        effect <- list(valid = FALSE, error = "Invalid design.")
    } else {
        if (effect %in% get_term_labels(design)) {
            effect <- list(value = effect, valid = TRUE)
        } else {
            effect <- list(valid = FALSE, error = "Effect is not part of the design.")
        }
    }
    effect
}

#' @keywords internal
get_term_labels <- function(design, string = FALSE) {
    term_labels <- attr(terms(design$idesign), "term.labels")

    if (string) {
        # term_labels <- paste0("[", paste0(term_labels, collapse = ", "), "]")
        term_labels <- paste0(term_labels, collapse = ", ")
    }

    term_labels

}

#' @keywords internal
get_idata_order <- function(design, string = TRUE) {
    order_idata <- apply(design$idata, 1, paste0, collapse = "")

    if (string) {
        # order_idata <- paste0(
        #     "[",
        #     paste0(order_idata, collapse = ", "),
        #     "]"
        # )
        order_idata <- paste0(order_idata, collapse = ", ")
    }

    order_idata
}

#' @keywords internal
is_valid <- function(parsed) {
    all(sapply(parsed, function(x) x$valid))
}

#' @keywords internal
collect_errors <- function(parsed) {
    errors <- sapply(parsed, function(x) x$error)
    errors <- errors[!sapply(errors, is.null)]
    errors <- paste0(1:length(errors), ". ", names(errors), ": ", errors)
    result <- list("Errors have occured")
    for (error in errors) {
        result <- list(result, br(), error)
    }
    result
}

# input <- list(
#     txt_design_eq = "2x",
#     txt_effect_eq = "A",
#     txt_n_eq = 10,
#     txt_mu_eq = "1,1,1,1",
#     txt_var_eq = "1",
#     txt_cov_eq = "0.5",
#     txt_sigma_uneq = "[[1,0.5,0.5,0.5],[0.5,1,0.5,0.5],[0.5,0.5,1,0.5],[0.5,0.5,0.5,1]]"
# )

#' @keywords internal
parse_input_covert_cohens_d <- function(input) {
    cohens_d <- parse_dbl(input$txt_convert_cohens_d)
    n <- parse_int(input$txt_convert_cohens_d_n)
    population <- parse_bool(input$chbox_convert_cohens_d_pop)

    list(cohens_d = cohens_d, n = n, population = population)
}

#' @keywords internal
parse_input_covert_eta <- function(input) {
    p_eta_sq <- parse_dbl(input$txt_convert_eta)
    n <- parse_int(input$txt_convert_eta_n)
    population <- parse_bool(input$chbox_convert_eta_pop)

    list(p_eta_sq = p_eta_sq, n = n, population = population)
}

#' @keywords internal
parse_input_power_cohens_d <- function(input) {

    n <- parse_int(input$txt_n_cohens_d)
    cohens_d <- parse_dbl(input$txt_cohens_d)
    alpha <- parse_dbl(input$txt_alpha_cohens_d)

    list(n = n, cohens_d = cohens_d, alpha = alpha)
}

#' @keywords internal
parse_input_power_eta <- function(input) {
    n <- parse_int(input$txt_n_eta)
    p_eta_sq <- parse_dbl(input$txt_p_eta_sq)
    df1 <- parse_int(input$txt_df1_eta)
    alpha <- parse_dbl(input$txt_alpha_eta)

    list(n = n, p_eta_sq = p_eta_sq, df1 = df1, alpha = alpha)
}

#' @keywords internal
parse_input_power_eq <- function(input) {

    design <- parse_design(input$txt_design_eq)

    effect <- parse_effect(input$txt_effect_eq, design)

    n <- parse_int(input$txt_n_eq)

    mu <- parse_dbl(input$txt_mu_eq, n = nrow(design$idata))

    var <- parse_dbl(input$txt_var_eq, n = 1)

    cov <- parse_dbl(input$txt_cov_eq, n = 1)

    alpha <- parse_dbl(input$txt_alpha_eq)

    list(design = design, effect = effect, mu = mu, n = n, var = var, cov = cov, alpha = alpha)
}

#' @keywords internal
parse_input_power_uneq <- function(input) {

    design <- parse_design(input$txt_design_uneq)

    effect <- parse_effect(input$txt_effect_uneq, design)

    n <- parse_int(input$txt_n_uneq)

    mu <- parse_dbl(input$txt_mu_uneq, n = nrow(design$idata))

    Sigma <- parse_dbl_matrix(input$txt_sigma_uneq, dim = nrow(design$idata))

    alpha <- parse_dbl(input$txt_alpha_uneq)

    list(design = design, effect = effect, mu = mu, n = n, Sigma = Sigma, alpha = alpha)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Power Calculator for Repeated Measures ANOVA"),

    withMathJax(),

    sidebarLayout(
        sidebarPanel(
            p(
                code("powerANOVA"),
                ' is an R package intended as a companion to the article: '
            ),
            p(
                "Langenberg, B., Janczyk, M., Koob, V., Kliegl, R., & Mayer, A. (2022). A tutorial on using the paired t-test for power calculations in repeated measures ANOVA with interactions. ",
                em(
                    "Behavior Research Methods",
                    .noWS = c("after-begin", "before-end", "outside", "after", "before")
                ),
                ". ",
                tags$a(
                    href="https://doi.org/10.3758/s13428-022-01902-8",
                    "https://doi.org/10.3758/s13428-022-01902-8",
                    .noWS = c("after-begin", "before-end", "outside", "after", "before")
                ),
                ".",
                .noWS = c("after-begin", "before-end", "outside", "after", "before")
            ),
            p("Package Version: ", as.character(packageVersion("powerANOVA")))
        ),
        mainPanel(
            tabsetPanel(
                tabPanel(
                    "Cohen's \\(d\\)",
                    br(),
                    sidebarLayout(
                        sidebarPanel(
                            textInput("txt_n_cohens_d",
                                      label = h3("N"),
                                      value = "20"),
                            textInput(
                                "txt_cohens_d",
                                label = h3("Cohen's \\(d\\)"),
                                value = "0.5"
                            ),
                            textInput(
                                "txt_alpha_cohens_d",
                                label = h3("\\(\\alpha\\)"),
                                value = "0.05"
                            )
                        ),

                        # Show a plot of the generated distribution
                        mainPanel(
                            plotOutput("power_plot_cohens_d"),
                            br(),
                            uiOutput("power_cohens_d")
                        )
                    )
                ),
                tabPanel(
                    "\\(\\eta^2_p\\)",
                    br(),
                    sidebarLayout(
                        sidebarPanel(
                            textInput("txt_n_eta",
                                      label = h3("N"),
                                      value = "20"),
                            textInput(
                                "txt_p_eta_sq",
                                label = h3("\\(\\eta^2_p\\)"),
                                value = as.character(0.5^2/(1+0.5^2)) # d^2*N/(df2+d^2*N)
                            ),
                            textInput(
                                "txt_df1_eta",
                                label = h3("\\(df_1\\)"),
                                value = "1"
                            ),
                            textInput(
                                "txt_alpha_eta",
                                label = h3("\\(\\alpha\\)"),
                                value = "0.05"
                            )
                        ),

                        # Show a plot of the generated distribution
                        mainPanel(
                            plotOutput("power_plot_eta"),
                            br(),
                            uiOutput("power_eta")
                        )
                    )
                ),
                tabPanel(
                    "means and equal (co)variances (\\(\\mu + \\sigma^2_i + \\sigma_{ij}\\))",
                    br(),
                    sidebarLayout(
                        sidebarPanel(
                            textInput("txt_design_eq",
                                      label = h3("Design"),
                                      value = "2x2"),
                            uiOutput("txt_out_effect_eq"),
                            textInput("txt_effect_eq",
                                      label = h3("Effect"),
                                      value = "A"),
                            textInput("txt_n_eq",
                                      label = h3("N"),
                                      value = "20"),
                            uiOutput("txt_out_order_eq"),
                            textInput(
                                "txt_mu_eq",
                                label = h3("\\(\\mu\\)"),
                                value = "0,0.35,0,0.35"
                            ),
                            textInput(
                                "txt_var_eq",
                                label = h3("\\(\\sigma^2_i\\)"),
                                value = "1"
                            ),
                            textInput(
                                "txt_cov_eq",
                                label = h3("\\(\\sigma_{ij}\\)"),
                                value = "0.5"
                            ),
                            uiOutput("txt_out_sigma_eq"),
                            textInput(
                                "txt_alpha_eq",
                                label = h3("\\(\\alpha\\)"),
                                value = "0.05"
                            )
                        ),

                        # Show a plot of the generated distribution
                        mainPanel(
                            plotOutput("power_plot_eq"),
                            br(),
                            uiOutput("power_eq")
                        )
                    )
                ),

                tabPanel(
                    "means and UNequal (co)variances (\\(\\mu + \\Sigma\\))",
                    br(),
                    sidebarLayout(
                        sidebarPanel(
                            textInput("txt_design_uneq",
                                      label = h3("Design"),
                                      value = "2x2"),
                            uiOutput("txt_out_effect_uneq"),
                            textInput("txt_effect_uneq",
                                      label = h3("Effect"),
                                      value = "A"),
                            textInput("txt_n_uneq",
                                      label = h3("N"),
                                      value = "20"),
                            uiOutput("txt_out_order_uneq"),
                            textInput(
                                "txt_mu_uneq",
                                label = h3("\\(\\mu\\)"),
                                value = "0,0.35,0,0.35"
                            ),
                            uiOutput("txt_out_order_sigma_uneq"),
                            textInput(
                                "txt_sigma_uneq",
                                label = h3("\\(\\Sigma\\)"),
                                value = "[[1,0.5,0.5,0.5],[0.5,1,0.5,0.5],[0.5,0.5,1,0.5],[0.5,0.5,0.5,1]]"
                            ),
                            uiOutput("txt_out_sigma_uneq"),
                            textInput(
                                "txt_alpha_uneq",
                                label = h3("\\(\\alpha\\)"),
                                value = "0.05"
                            )
                        ),

                        mainPanel(
                            plotOutput("power_plot_uneq"),
                            br(),
                            uiOutput("power_uneq"),
                        )
                    )
                ),

                tabPanel(
                    "effect size converter",
                    br(),
                    fluidRow(
                        column(
                            4,
                            sidebarPanel(width = 12,
                                textInput(
                                    "txt_convert_cohens_d",
                                    label = h3("Cohen's \\(d\\)"),
                                    value = "0.5"
                                ),
                                textInput(
                                    "txt_convert_cohens_d_n",
                                    label = h3("\\(N\\)"),
                                    value = "20"
                                ),
                                checkboxInput("chbox_convert_cohens_d_pop",
                                    label = "Population",
                                    value = TRUE),
                                hr(),
                                h3("\\(\\eta^2_p\\)"),
                                uiOutput("txt_out_cohens_d_to_eta")
                            ),
                        ),
                        column(
                            4,
                            sidebarPanel(width = 12,
                                textInput(
                                    "txt_convert_eta",
                                    label = h3("\\(\\eta^2_p\\)"),
                                    value = as.character(0.5^2/(1 + 0.5^2)) # d^2*N/(df2+d^2*N)
                                ),
                                textInput(
                                    "txt_convert_eta_n",
                                    label = h3("\\(N\\)"),
                                    value = "20"
                                ),
                                checkboxInput("chbox_convert_eta_pop",
                                              label = "Population",
                                              value = TRUE),
                                hr(),
                                h3("Cohen's \\(d\\)"),
                                uiOutput("txt_out_eta_cohens_d")
                            )
                        )
                    )
                )
            )
        )
    )
)

# Define server logic required to draw a histogram
#' @keywords internal
server <- function(input, output) {

    output$txt_out_sigma_eq <- renderUI({
        parsed <- parse_input_power_eq(input)
        if (is_valid(parsed[c("design", "effect", "cov", "var")])) {
            Sigma <-
                matrix(
                    parsed$cov$value,
                    ncol = nrow(parsed$design$idata),
                    nrow = nrow(parsed$design$idata)
                )
            diag(Sigma) <- parsed$var$value
            Sigma <- apply(Sigma, 1, paste0, collapse = " & ")
            Sigma <- paste0(Sigma, collapse = " \\\\ ")
            withMathJax(
                h3("\\(\\Sigma\\)"),
                sprintf(
                    "$$\\mathbf{\\Sigma} = \\begin{pmatrix} %s \\end{pmatrix}$$",
                    Sigma
                )
            )
        }
    })

    output$txt_out_sigma_uneq <- renderUI({
        parsed <- parse_input_power_uneq(input)
        if (is_valid(parsed["Sigma"])) {
            Sigma <- parsed$Sigma$value
            Sigma <- apply(Sigma, 1, paste0, collapse = " & ")
            Sigma <- paste0(Sigma, collapse = " \\\\ ")
            withMathJax(
                sprintf(
                    "$$\\mathbf{\\Sigma} = \\begin{pmatrix} %s \\end{pmatrix}$$",
                    Sigma
                )
            )
        }
    })

    output$txt_out_cohens_d_to_eta <- renderUI({
        parsed <- parse_input_covert_cohens_d(input)

        if (is_valid(parsed)) {
            p_eta_sq <-
                convert_cohens_d_petasq(
                    cohens_d = parsed$cohens_d$value,
                    n = parsed$n$value,
                    parsed$population$value
                )
            withMathJax(sprintf("\\(\\eta^2_p = %f\\)", p_eta_sq))
        }
    })

    output$txt_out_eta_cohens_d <- renderUI({
        parsed <- parse_input_covert_eta(input)

        if (is_valid(parsed)) {
            cohens_d <-
                convert_petasq_cohens_d(
                    p_eta_sq = parsed$p_eta_sq$value,
                    n = parsed$n$value,
                    parsed$population$value
                )
            withMathJax(sprintf("Cohen's \\(d = %f\\)", cohens_d))
        }
    })

    output$power_plot_cohens_d <- renderPlot({
        parsed <- parse_input_power_cohens_d(input)

        if (is_valid(parsed)) {
            power_plot_cohens_d(
                n = parsed$n$value,
                cohens_d = parsed$cohens_d$value,
                alpha = parsed$alpha$value
            )

        }
    })

    output$power_plot_eta <- renderPlot({
        parsed <- parse_input_power_eta(input)

        if (is_valid(parsed)) {
            if (parsed$df1$value == 1L) {
                power_plot_cohens_d(
                    n = parsed$n$value,
                    cohens_d = convert_petasq_cohens_d(p_eta_sq = parsed$p_eta_sq$value),
                    alpha = parsed$alpha$value
                )
            } else {
                power_plot_p_eta_sq(
                    n = parsed$n$value,
                    p_eta_sq = parsed$p_eta_sq$value,
                    df1 = parsed$df1$value,
                    alpha = parsed$alpha$value
                )
            }

        }
    })

    output$power_plot_eq <- renderPlot({
        parsed <- parse_input_power_eq(input)

        print(parsed)

        if (is_valid(parsed)) {
            power_plot_mu_cov(
                n = parsed$n$value,
                mu = parsed$mu$value,
                var = parsed$var$value,
                cov = parsed$cov$value,
                effect = parsed$effect$value,
                idata = parsed$design$idata,
                alpha = parsed$alpha$value
            )
        }
    })

    output$power_plot_uneq <- renderPlot({
        parsed <- parse_input_power_uneq(input)
        if (is_valid(parsed)) {
            power_plot_mu_cov(
                n = parsed$n$value,
                mu = parsed$mu$value,
                Sigma = parsed$Sigma$value,
                effect = parsed$effect$value,
                idata = parsed$design$idata,
                alpha = parsed$alpha$value
            )
        }
    })

    output$power_cohens_d <- renderUI({

        parsed <- parse_input_power_cohens_d(input)

        if (is_valid(parsed)) {
            pwr <- power_cohens_d(
                n = parsed$n$value,
                cohens_d = parsed$cohens_d$value,
                alpha = parsed$alpha$value
            )
            withMathJax(sprintf("\\(PWR = %f\\)", pwr))
        } else {
            do.call(withMathJax, collect_errors(parsed))
        }
    })

    output$power_eta <- renderUI({

        parsed <- parse_input_power_eta(input)

        if (is_valid(parsed)) {
            pwr <- power_petasq(
                n = parsed$n$value,
                p_eta_sq = parsed$p_eta_sq$value,
                df1 = parsed$df1$value,
                alpha = parsed$alpha$value
            )
            withMathJax(sprintf("\\(PWR = %f\\)", pwr))
        } else {
            do.call(withMathJax, collect_errors(parsed))
        }
    })

    output$power_eq <- renderUI({
        parsed <- parse_input_power_eq(input)
        if (is_valid(parsed)) {
            pwr <- power_mu_cov(
                n = parsed$n$value,
                mu = parsed$mu$value,
                var = parsed$var$value,
                cov = parsed$cov$value,
                effect = parsed$effect$value,
                idata = parsed$design$idata,
                alpha = parsed$alpha$value
            )
            withMathJax(sprintf("\\(PWR = %f\\)", pwr))
        } else {
            do.call(withMathJax, collect_errors(parsed))
        }
    })

    output$power_uneq <- renderUI({
        parsed <- parse_input_power_uneq(input)
        if (is_valid(parsed)) {
            pwr <- power_mu_cov(
                n = parsed$n$value,
                mu = parsed$mu$value,
                Sigma = parsed$Sigma$value,
                effect = parsed$effect$value,
                idata = parsed$design$idata,
                alpha = parsed$alpha$value
            )
            withMathJax(sprintf("\\(PWR = %f\\)", pwr))
        } else {
            do.call(withMathJax, collect_errors(parsed))
        }
    })

    output$txt_out_effect_eq <- renderUI({
        parsed <- parse_input_power_eq(input)
        if (is_valid(parsed)) {
            withMathJax(sprintf(
                "Choose one from the following effects: %s",
                get_term_labels(parsed$design, string = T)
            ))
        } else {
            withMathJax(sprintf("%s", parsed$design$error))
        }
    })

    output$txt_out_order_eq <- renderUI({
        parsed <- parse_input_power_eq(input)
        if (is_valid(parsed)) {
            withMathJax(sprintf(
                "\\(\\mu\\) must conform to this order: %s",
                get_idata_order(parsed$design, string = T)
            ))
        } else {
            withMathJax(sprintf("%s", parsed$design$error))
        }
    })

    output$txt_out_effect_uneq <- renderUI({
        parsed <- parse_input_power_uneq(input)
        if (is_valid(parsed)) {
            withMathJax(sprintf(
                "Choose one from the following effects: %s",
                get_term_labels(parsed$design, string = T)
            ))
        } else {
            withMathJax(sprintf("%s", parsed$design$error))
        }
    })

    output$txt_out_order_uneq <- renderUI({
        parsed <- parse_input_power_uneq(input)
        if (is_valid(parsed)) {
            withMathJax(sprintf(
                "\\(\\mu\\) must conform to this order: %s",
                get_idata_order(parsed$design, string = T)
            ))
        } else {
            withMathJax(sprintf("%s", parsed$design$error))
        }
    })

    output$txt_out_order_sigma_uneq <- renderUI({
        parsed <- parse_input_power_uneq(input)
        if (is_valid(parsed)) {
            withMathJax(sprintf(
                "Rows and columns of \\(\\Sigma\\) must conform to this order: %s",
                get_idata_order(parsed$design, string = T)
            ))
        } else {
            withMathJax(sprintf("%s", parsed$design$error))
        }
    })

}
