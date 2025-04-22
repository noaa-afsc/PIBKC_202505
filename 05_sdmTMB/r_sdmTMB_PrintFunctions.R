getMaxGradient<-function(x){
  grd = x$tmb_obj$gr();
  grd[(abs(grd)==max(abs(grd)))[1]];
}

print_sdmTMB<-function (x, ...){
    sdmTMB:::reinitialize(x)
    lp <- x$tmb_obj$env$last.par.best
    r <- x$tmb_obj$report(lp)
    delta <- isTRUE(x$family$delta)
    print_header(x)
    if (delta)
        cat("\nDelta/hurdle model 1: -----------------------------------\n")
    print_one_model(x, 1)
    if (delta) {
        cat("\nDelta/hurdle model 2: -----------------------------------\n")
        print_one_model(x, 2)
    }
    if (delta)
        cat("\n")
    print_footer(x)
}

print_one_model<-function (x, m = 1) {
    info <- get_model_print_info(x)
    main <- print_main_effects(x, m = m)
    smooth <- print_smooth_effects(x, m = m)
    iid_re <- print_iid_re(x, m = m)
    tv <- print_time_varying(x, m = m)
    range <- print_range(x, m = m)
    other <- print_other_parameters(x, m = m)
    if (m == 1)
        cat(info$family1, "\n")
    if (m == 2)
        cat(info$family2, "\n")
    print(rbind(main, smooth$smooth_effects))
    cat("\n")
    if (!is.null(smooth$smooth_sds)) {
        cat("Smooth terms:\n")
        print(smooth$smooth_sds)
        cat("\n")
    }
    if (!is.null(iid_re)) {
        cat("Random intercepts:\n")
        print(iid_re)
        cat("\n")
    }
    if (!is.null(tv)) {
        cat("Time-varying parameters:\n")
        print(tv)
        cat("\n")
    }
    cat(other$phi)
    cat(other$tweedie_p)
    cat(other$gengamma_par)
    cat(other$rho)
    cat(range)
    cat(other$sigma_O)
    cat(other$sigma_Z)
    cat(other$sigma_E)
}

print_other_parameters<-function (x, m = 1L) {
    b <- tidy(x, "ran_pars", model = m, silent = TRUE)
    get_term_text <- function(term_name = "", pretext = "") {
        b2 <- as.list(x$sd_report, what = "Estimate")
        if (term_name %in% b$term) {
            a <- mround(b$estimate[b$term == term_name], 2L)
            a <- paste0(pretext, ": ", a, "\n")
        }
        else if (term_name %in% names(b2)) {
            a <- mround(b2[[term_name]], 2L)
            a <- paste0(pretext, ": ", a, "\n")
        }
        else {
            a <- ""
        }
        a
    }
    phi <- get_term_text("phi", "Dispersion parameter")
    tweedie_p <- get_term_text("tweedie_p", "Tweedie p")
    gengamma_par <- if ("gengamma" %in% family(x)[[m]]) {
        get_term_text("gengamma_Q", "Generalized gamma Q")
    }
    else ""
    sigma_O <- get_term_text("sigma_O", "Spatial SD")
    xtra <- if (x$spatiotemporal[m] == "ar1")
        "marginal "
    else ""
    sigma_E <- get_term_text("sigma_E", paste0("Spatiotemporal ",
        xtra, toupper(x$spatiotemporal[m]), " SD"))
    rho <- get_term_text("rho", "Spatiotemporal AR1 correlation (rho)")
    if ("sigma_Z" %in% b$term) {
        sigma_Z <- x$tmb_obj$report(x$tmb_obj$env$last.par.best)$sigma_Z
        sigma_Z <- sigma_Z[, m, drop = TRUE]
        a <- mround(sigma_Z, 2L)
        sigma_Z <- paste0("Spatially varying coefficient SD (",
            x$spatial_varying, "): ", a, "\n")
        sigma_Z <- gsub("\\(\\(", "\\(", sigma_Z)
        sigma_Z <- gsub("\\)\\)", "\\)", sigma_Z)
        sigma_Z <- paste(sigma_Z, collapse = "")
    }
    else {
        sigma_Z <- ""
    }
    named_list(phi, tweedie_p, sigma_O, sigma_E, sigma_Z, rho,
        gengamma_par)
}

named_list<-function (...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1L]
    if (is.null(nm <- names(L))) {
        nm <- snm
    }
    if (any(nonames <- nm == "")) {
        nm[nonames] <- snm[nonames]
    }
    stats::setNames(L, nm)
}

print_anisotropy<-function (x, m = 1L, digits = 1L, return_dat = FALSE) {
    aniso_df <- plot_anisotropy(x, return_data = TRUE)
    aniso_df$degree <- aniso_df$angle * 180/pi
    if (isTRUE(x$family$delta)) {
        aniso_df_sp <- aniso_df[aniso_df$random_field == "spatial" &
            aniso_df$model_num == m, ][1, c("a", "b", "degree")]
        aniso_df_st <- aniso_df[aniso_df$random_field == "spatiotemporal" &
            aniso_df$model_num == m, ][1L, c("a", "b", "degree")]
    }
    else {
        if (x$spatial[m] != "off") {
            aniso_df_sp <- aniso_df[aniso_df$random_field ==
                "spatial", ][1L, c("a", "b", "degree")]
        }
        if (x$spatiotemporal[m] != "off") {
            aniso_df_st <- aniso_df[aniso_df$random_field ==
                "spatiotemporal", ][1L, c("a", "b", "degree")]
        }
    }
    if (x$spatial[m] != "off") {
        aniso_df_sp[1:2] <- mround(aniso_df_sp[1:2], digits)
        aniso_df_sp[3] <- mround(aniso_df_sp[3], 0L)
        aniso_df_sp[4] <- paste0("Matérn anisotropic range (spatial): ",
            aniso_df_sp[2], " to ", aniso_df_sp[1], " at ", aniso_df_sp[3],
            " deg.", "\n")
    }
    if (x$spatiotemporal[m] != "off") {
        aniso_df_st[1:2] <- mround(aniso_df_st[1:2], digits)
        aniso_df_st[3] <- mround(aniso_df_st[3], 0L)
        aniso_df_st[4] <- paste0("Matérn anisotropic range (spatiotemporal): ",
            aniso_df_st[2], " to ", aniso_df_st[1], " at ", aniso_df_st[3],
            " deg.", "\n")
    }
    if (x$spatial[m] == "on" && x$spatiotemporal[m] == "off") {
        range_text <- aniso_df_sp[[4]]
    }
    if (x$spatial[m] == "off" && x$spatiotemporal[m] != "off") {
        range_text <- aniso_df_st[[4]]
    }
    if (x$tmb_data$share_range[m] == 1L && x$spatial[m] == "on" &&
        x$spatiotemporal[m] != "off") {
        range_text <- aniso_df_sp[[4]]
    }
    if (x$tmb_data$share_range[m] == 0L && x$spatial[m] == "on" &&
        x$spatiotemporal[m] != "off") {
        range_text <- paste0(aniso_df_sp[[4]], aniso_df_st[[4]])
    }
    if (return_dat)
        return(list(sp = aniso_df_sp, st = aniso_df_st))
    range_text
}

print_range<-function (x, m = 1L, digits = 2L)
{
    b <- tidy(x, effects = "ran_pars", model = m, silent = TRUE)
    range <- b$estimate[b$term == "range"]
    if (is.null(range)) {
        return(NULL)
    }
    range <- mround(range, digits)
    range_text <- if (x$tmb_data$share_range[m]) {
        paste0("Matérn range: ", range[1], "\n")
    }
    else {
        paste0("Matérn range (spatial): ", range[1], "\n", "Matérn range (spatiotemporal): ",
            range[2], "\n")
    }
    if (as.logical(x$tmb_data$anisotropy)) {
        range_text <- print_anisotropy(x = x, m = m)
    }
    if (x$spatial[m] == "off" && x$spatiotemporal[m] == "off") {
        range_text <- NULL
    }
    range_text
}

print_time_varying<-function (x, m = 1) {
    if (!is.null(x$time_varying)) {
        sr <- x$sd_report
        sr_se <- as.list(sr, "Std. Error")
        sr_est <- as.list(sr, "Estimate")
        b_rw_t_est <- sr_est$b_rw_t[, , m]
        b_rw_t_se <- sr_se$b_rw_t[, , m]
        tv_names <- colnames(model.matrix(x$time_varying, x$data))
        mm_tv <- cbind(round(as.numeric(b_rw_t_est), 2L), round(as.numeric(b_rw_t_se),
            2L))
        colnames(mm_tv) <- c("coef.est", "coef.se")
        time_slices <- x$time_lu$time_from_data
        row.names(mm_tv) <- paste(rep(tv_names, each = length(time_slices)),
            time_slices, sep = "-")
        p <- tidy(x, effects = "ran_pars", model = m, silent = TRUE)
        if (any(p$term == "rho_time")) {
            rho_tv <- as.matrix(p[p$term == "rho_time", c("estimate",
                "std.error")])
            colnames(rho_tv) <- colnames(mm_tv)
            row.names(rho_tv) <- paste("rho", tv_names, sep = "-")
            mm_tv <- rbind(mm_tv, rho_tv)
        }
    }
    else {
        mm_tv <- NULL
    }
    mm_tv
}

print_iid_re<-function (x, m = 1) {
    .tidy <- tidy(x, "ran_pars", model = m, silent = TRUE)
    if ("sigma_G" %in% .tidy$term) {
        re_int_names <- x$split_formula[[1]]$barnames
        re_int_mat <- matrix(NA_real_, nrow = length(re_int_names),
            ncol = 1L)
        re_int_mat[, 1L] <- round(.tidy$estimate[.tidy$term ==
            "sigma_G"], 2L)
        rownames(re_int_mat) <- re_int_names
        colnames(re_int_mat) <- "Std. Dev."
    }
    else {
        re_int_mat <- NULL
    }
    re_int_mat
}

print_main_effects<-function (x, m = 1) {
    b <- as.data.frame(tidy(x, model = m, silent = TRUE))
    b$estimate <- round(b$estimate, 2L)
    b$std.error <- round(b$std.error, 2L)
    mm <- cbind(b$estimate, b$std.error)
    colnames(mm) <- c("coef.est", "coef.se")
    row.names(mm) <- b$term
    mm
}

print_smooth_effects<-function (x, m = 1, edf = NULL, silent = FALSE) {
    sr <- x$sd_report
    sr_se <- as.list(sr, "Std. Error")
    sr_est <- as.list(sr, "Estimate")
    if (x$tmb_data$has_smooths) {
        bs_se <- round(sr_se$bs[, m], 2L)
        bs <- round(sr_est$bs[, m], 2L)
        sm <- parse_smoothers(formula = x$formula[[m]], data = x$data)
        sm_names_orig <- unlist(lapply(sm$Zs, function(x) attr(x,
            "s.label")))
        sm_names <- gsub("s\\(", "", sm_names_orig)
        sm_names <- gsub("t2\\(", "", sm_names)
        sm_names <- gsub("\\)$", "", sm_names)
        sm_names_bs <- paste0("s", sm_names)
        sm_classes <- unlist(sm$classes)
        xx <- lapply(sm_names_bs, function(.x) {
            n_sm <- grep(",", .x) + 1
            if (length(n_sm)) {
                .x <- paste(.x, seq_len(n_sm), sep = "_")
                gsub(",", "", .x)
            }
            else {
                .x
            }
        })
        sm_names_bs <- unlist(xx)
        sm_names_sds <- paste0("sds(", sm_names, ")")
        mm_sm <- cbind(bs, bs_se)
        .sm_names_bs <- sm_names_bs[!sm_classes %in% c("cc.smooth.spec",
            "cs.smooth.spec")]
        if (length(.sm_names_bs) != nrow(mm_sm)) {
            if (nrow(mm_sm) > 0L) {
                msg <- c("Smoother fixed effect matrix names could not be retrieved.",
                  "This could be from using t2(), setting m != 2 in s(), or other unanticipated s() arguments.",
                  "This does not affect model fitting.", "We'll use generic covariate names ('scovariate') here instead.",
                  "")
                if (!silent)
                  cli_inform(msg)
                row.names(mm_sm) <- paste0("scovariate-", seq_len(nrow(mm_sm)))
            }
            else {
                mm_sm <- NULL
            }
        }
        else {
            row.names(mm_sm) <- .sm_names_bs
        }
        smooth_sds <- round(exp(sr_est$ln_smooth_sigma[, m]),
            2L)
        re_sm_mat <- matrix(NA_real_, nrow = length(smooth_sds),
            ncol = 1L)
        re_sm_mat[, 1] <- smooth_sds
        rownames(re_sm_mat) <- sm_names_sds
        colnames(re_sm_mat) <- "Std. Dev."
        if (!is.null(edf)) {
            if (is_delta(x)) {
                lp_regex <- paste0("^", m, "LP-s\\(")
                edf <- edf[grepl(lp_regex, names(edf))]
            }
            else {
                edf <- edf[grepl("^s\\(", names(edf))]
            }
            edf <- round(edf, 2)
            re_sm_mat <- cbind(re_sm_mat, matrix(edf, ncol = 1))
            colnames(re_sm_mat)[2] <- "EDF"
        }
    }
    else {
        re_sm_mat <- NULL
        mm_sm <- NULL
    }
    list(smooth_effects = mm_sm, smooth_sds = re_sm_mat)
}

parse_smoothers<-function (formula, data, knots = NULL,
                           newdata = NULL, basis_prev = NULL) {
    terms <- all_terms(formula)
    smooth_i <- get_smooth_terms(terms)
    basis <- list()
    basis_out <- list()
    Zs <- list()
    Xs <- list()
    labels <- list()
    classes <- list()
    if (length(smooth_i) > 0) {
        has_smooths <- TRUE
        smterms <- terms[smooth_i]
        ns <- 0
        ns_Xf <- 0
        for (i in seq_along(smterms)) {
            if (grepl("bs\\=\\\"re", smterms[i]))
                cli_abort("bs = 're' is not currently supported for smooths")
            if (grepl("fx\\=T", smterms[i]))
                cli_abort("fx = TRUE is not currently supported for smooths")
            obj <- eval(str2expression(smterms[i]))
            labels[[i]] <- obj$label
            classes[[i]] <- attr(obj, "class")
            if (is.null(newdata)) {
                basis_out[[i]] <- basis[[i]] <- mgcv::smoothCon(object = obj,
                  data = data, knots = knots, absorb.cons = TRUE,
                  modCon = 3, diagonal.penalty = FALSE)
            }
            else {
                basis[[i]] <- basis_prev[[i]]
            }
            for (j in seq_along(basis[[i]])) {
                ns_Xf <- ns_Xf + 1
                rasm <- mgcv::smooth2random(basis[[i]][[j]],
                  names(data), type = 2)
                if (!is.null(newdata)) {
                  rasm <- s2rPred(basis[[i]][[j]], rasm, newdata)
                }
                for (k in seq_along(rasm$rand)) {
                  ns <- ns + 1
                  Zs[[ns]] <- rasm$rand[[k]]
                }
                Xs[[ns_Xf]] <- rasm$Xf
            }
        }
        sm_dims <- unlist(lapply(Zs, ncol))
        Xs <- do.call(cbind, Xs)
        b_smooth_start <- c(0, cumsum(sm_dims)[-length(sm_dims)])
    }
    else {
        has_smooths <- FALSE
        sm_dims <- 0L
        b_smooth_start <- 0L
        Xs <- matrix(nrow = 0L, ncol = 0L)
    }
    list(Xs = Xs, Zs = Zs, has_smooths = has_smooths, labels = labels,
        classes = classes, basis_out = basis_out, sm_dims = sm_dims,
        b_smooth_start = b_smooth_start)
}

all_terms<-function (x) {
    if (!length(x)) {
        return(character(0))
    }
    if (!inherits(x, "terms")) {
        x <- terms(stats::as.formula(x))
    }
    rm_wsp(attr(x, "term.labels"))
}

rm_wsp<-function (x) {
    out <- gsub("[ \t\r\n]+", "", x, perl = TRUE)
    dim(out) <- dim(x)
    out
}

get_smooth_terms<-function (terms) {
    x1 <- grep("s\\(", terms)
    x2 <- grep("t2\\(", terms)
    c(x1, x2)
}


print_header<-function (x) {
    info <- get_model_print_info(x)
    cat(info$title)
    cat(info$formula)
    cat(info$mesh)
    cat(info$time)
    cat(info$data)
    cat(info$overall_family)
}

get_model_print_info<-function(x){
    delta <- isTRUE(x$family$delta)
    spatial_only <- as.logical(x$tmb_data$spatial_only)
    fit_by <- if (isTRUE(x$reml))
        "REML"
    else "ML"
    if (all(spatial_only)) {
        title <- paste0("Spatial model fit by ", fit_by, " ['sdmTMB']\n")
    }
    else {
        title <- paste0("Spatiotemporal model fit by ", fit_by,
            " ['sdmTMB']\n")
    }
    if (all(x$spatial == "off") && all(x$spatiotemporal == "off")) {
        title <- paste0("Model fit by ", fit_by, " ['sdmTMB']\n")
    }
    aniso <- as.logical(x$tmb_data$anisotropy)
    covariance <- paste0(if (aniso)
        "anisotropic"
    else "isotropic")
    formula <- paste0("Formula: ", deparse(x$call$formula), "\n")
    if (deparse(x$call$time) != "NULL") {
        time <- paste0("Time column: ", deparse(x$call$time),
            "\n")
        time <- gsub("\\\"", "", time)
        time <- gsub("\\'", "", time)
    }
    else {
        time <- NULL
    }
    mesh <- paste0("Mesh: ", deparse(x$call$mesh), " (", covariance,
        " covariance)\n")
    data <- paste0("Data: ", deparse(x$call$data), "\n")
    if (length(mesh) > 1L)
        mesh <- NULL
    if (length(data) > 1L)
        data <- NULL
    if ("clean_name" %in% names(x$family)) {
        overall_family <- x$family$clean_name
    }
    else {
        overall_family <- paste0(x$family$family[1], "(link = '",
            x$family$link[1], "')")
    }
    overall_family <- paste0("Family: ", overall_family, "\n")
    if (delta) {
        family1 <- paste0("Family: ", x$family$family[1], "(link = '",
            x$family$link[1], "')")
        family2 <- paste0("Family: ", x$family$family[2], "(link = '",
            x$family$link[2], "')")
    }
    else {
        family1 <- NULL
        family2 <- NULL
    }
    criterion <- paste0(fit_by, " criterion at convergence: ",
        mround(x$model$objective, 3), "\n")
    sdmTMB:::named_list(delta, spatial_only, title, time, mesh, formula,
                        data, family1, family2, overall_family, criterion, covariance)
}

mround<-function (x, digits) {
    sprintf(paste0("%.", digits, "f"), round(x, digits))
}

print_footer<-function (x,verbose=FALSE) {
    info <- sdmTMB:::print_model_info(x)
    cat(info$criterion);
    if (verbose){
      cat("\nSee ?tidy.sdmTMB to extract these values as a data frame.\n")
      if (as.logical(x$tmb_data$anisotropy)) {
          cat("See ?plot_anisotropy to plot the anisotropic range.\n")
      }
      sink(tempfile())
      suppressWarnings(suppressMessages(s <- sanity(x)))
      sink()
      if (!all(unlist(s))) {
          cat("\n**Possible issues detected! Check output of sanity().**\n")
      }
    }
}
