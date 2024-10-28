
PEACE_harmonize <- function(data, model, nocovar = FALSE, return_avg = FALSE) {

    # use fitted PEACE model to harmonize data;
    # can be a test set

    softmax_stable <- function(x) {
        max_x <- max(x)
        exp_x <- exp(x - max_x)
        return(exp_x / sum(exp_x))
    }

    log_phi <- log(model$phi)
    log_1mphi <- log(1-model$phi)

    n_iter <- dim(model$phi)[1]
    n_data_points <- dim(data$y)[1]

    iter_list <- list()
    y_combat_avg <- NULL
    for (i_iter in 1:n_iter) {
        row_list <- list()
        for (i in 1:n_data_points) {

            b <- data$batch[i]
            if (b == 1) {
                row_list[[i]] <- data$y[i,]
            } else {
                y_i <- data$y[i,]

                beta_x1 <- 0
                beta_x2 <- 0
                if (!nocovar) {
                    for (m in 1:ncol(data$x)) {
                        beta_x1 <- beta_x1 + (data$x[i,m] * model$beta[i_iter,1,m,])
                        beta_x2 <- beta_x2 + (data$x[i,m] * model$beta[i_iter,2,m,])
                    }
                    
                }

                alpha1_beta_x <- beta_x1 + model$alpha[i_iter,1,]
                alpha2_beta_x <- beta_x2 + model$alpha[i_iter,2,]

                tmp1 <- log_phi[i_iter,b] + sum(dnorm(y_i, alpha1_beta_x + model$gamma[i_iter,b-1,], model$delta_sigma1[i_iter,b-1,], log=TRUE))
                tmp2 <- log_1mphi[i_iter,b] + sum(dnorm(y_i, alpha2_beta_x + model$gamma[i_iter,b-1,], model$delta_sigma2[i_iter,b-1,], log=TRUE))
                r <- softmax_stable(c(tmp1,tmp2))

                y_combat1 <- (y_i - alpha1_beta_x - model$gamma[i_iter,b-1,]) / model$delta[i_iter,b-1,] + alpha1_beta_x
                y_combat2 <- (y_i - alpha2_beta_x - model$gamma[i_iter,b-1,]) / model$delta[i_iter,b-1,] + alpha2_beta_x
                y_combat <- r[1] * y_combat1 + r[2] * y_combat2

                row_list[[i]] <- y_combat
            }
        }
        y_combat_df <- do.call(rbind, row_list) %>% as_tibble()

        if (return_avg) {
            if (is.null(y_combat_avg)) {
                y_combat_avg <- y_combat_df
            } else {
                y_combat_avg <- y_combat_avg + (y_combat_df - y_combat_avg) / i_iter  # update average with new value
            }
        } else {
            iter_list[[i_iter]] <- y_combat_df
        }
    }

    if (return_avg) {
        return(y_combat_avg)
    } else {
        return(iter_list)
    }

}

get_PEACE_data <- function(df, covars) {

    # filter out NA
    if (covars == "age") {
        df <- df %>% drop_na(age)
    } else if (covars == "age_sex_apoe") {
        df <- df %>% drop_na(age, sex, apoe)
    }

    features <- df %>% select(all_of(feature_names)) %>% as.matrix
    batches <- df %>% pull(tracer) %>% case_match("PIB"~1,"AV45"~2)
    age <- c(scale(df %>% pull(age)))  # z-score, convert to vector
    sex <- df %>% pull(sex) %>% case_match("M"~1,"F"~0)
    apoe <- df %>% pull(apoe)

    data <- list(
        I = nrow(features),
        K = ncol(features),
        B = length(unique(batches)),
        y = features,
        batch = batches
    )

    if (covars == "age_sex_apoe") {
        covar <- cbind(age, sex, apoe)
        data[["M"]] <- ncol(covar)
        data[["x"]] <- covar
    } else if (covars == "age") {
        covar <- cbind(age)
        data[["M"]] <- ncol(covar)
        data[["x"]] <- covar
    }

    return(data)

}

apply_PEACE_wrapper <- function(model, data, df, harmonization_method_str, nocovar = FALSE) {
    # memory-efficient computation of average of y_combat over PEACE iterations

    # extract model parameters: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
    model_extract <- extract(model, permuted = TRUE)

    # apply to crossover data
    y_peace_avg <- PEACE_harmonize(data, model_extract, nocovar, return_avg = TRUE)
    peace_df <- df %>%
        select(-all_of(feature_names)) %>%
        bind_cols(y_peace_avg) %>%
        mutate(harmonization_method = harmonization_method_str)

    return(peace_df)
}
