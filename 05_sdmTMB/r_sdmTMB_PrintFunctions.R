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

print_header<-function (x) {
    info <- print_model_info(x)
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
