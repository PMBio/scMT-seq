library(weights)
library(dplyr)
library(tidyr)

log_counts <- function(x) {
  return (log10(x + 1))
}

# weighted correlation
wtd_cor <- function(x, y, weights=NULL, method='pearson',
  alternative='two.sided', n_min=3) {
  r <- data.frame(r=NA, r_lo=NA, r_up=NA, p=NA, n=NA, n_wtd=NA)
  o <- !is.na(x) & !is.na(y)
  x <- x[o]
  y <- y[o]
  weights <- weights[o]
  r$n <- sum(o)
  if (is.null(weights)) {
    r$n_wtd <- r$n
  } else {
    r$n_wtd <- sum(weights)
  }
  if (r$n >= n_min & var(x) > 0 & var(y) > 0) {
    ct <- NULL
    ct_wtd <- NULL
    if (is.null(weights)) {
      ct <- cor.test(x, y, alternative=alternative, method=method)
    } else {
      if (method == 'spearman') {
        rx <- wtd.rank(x, weights=weights)
        ry <- wtd.rank(y, weights=weights)
        ct <- cor.test(rx, ry, alternative=alternative, method=method)
      } else {
        ct_wtd <- wtd.cor(x, y, weight=weights)
      }
    }
    if (!is.null(ct)) {
      r$r <- ct$estimate[1]
      r$p <- ct$p.value
      if (!is.null(ct$conf.int)) {
        r$r_lo <- ct$conf.int[1]
        r$r_up <- ct$conf.int[2]
      }
    } else {
      # TODO: replace 1.8
      r$r <- ct_wtd[1]
      r$p <- ct_wtd[4]
      r$r_lo <- max(-1, ct_wtd[1] - 1.8 * ct_wtd[2])
      r$r_up <- min(1, ct_wtd[1] + 1.8 * ct_wtd[2])
    }
  }
  return (r)
}

# mean imputation by columns
impute <- function(d) {
  means <- colMeans(d, na.rm=T)
  if (any(is.na(means))) {
    stop('Insufficient data for mean imputation!')
  }
  for (i in 1:length(means)) {
    d[is.na(d[,i]), i] <- means[i]
  }
  return (d)
}

# adjust batch effect
adjust_lm <- function(y, x) {
  m <- lm(y ~ x)
  return (residuals(m))
}

adjust_df <- function(d) {
  d$yr <- adjust_lm(d$y, d$x)
  return (d)
}

pca <- function(d) {
  # columns are samples
  d <- scale(d, center=T, scale=F)
  d <- t(d)
  s <- svd(d)
  vec <- s$u
  rownames(vec) <- rownames(d)
  val <- s$d**2
  val <- val / sum(val)
  return (list(vec=vec, val=val))
}


# rbinds data frames by common columns
rbind_frames <- function(d) {
  cols <- colnames(d[[1]])
  for (dd in d) {
    cols <- intersect(cols, colnames(dd))
  }
  dcols <- lapply(d, function(x) subset(x, select=cols))
  d <- do.call(rbind.data.frame, dcols)
  return (d)
}

# reads several files by fun and concats them
read_all <- function(filenames, fun) {
  d <- list()
  for (filename in filenames) {
    h <- fun(filename)
    h$file <- basename(filename)
    d[[length(d) + 1]] <- h
  }
  h <- rbind_frames(d) %>% mutate(file=factor(file))
  return (h)
}
