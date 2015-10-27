library(weights)
library(dplyr)
library(tidyr)
library(stringr)

theme_pub <- function() {
  p <- theme(
    axis.text=element_text(size=rel(1.2), color='black'),
    axis.title=element_text(size=rel(1.5)),
    # axis.title.y=element_text(vjust=1.5),
    # axis.title.x=element_text(vjust=-0.2),
    legend.position='top',
    legend.text=element_text(size=rel(1.2)),
    legend.title=element_text(size=rel(1.2)),
    legend.key=element_rect(fill='transparent'),
    panel.border=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour="black", size=1),
    axis.ticks.length = unit(.3, 'cm'),
    axis.ticks.margin = unit(.3, 'cm')
    )
  return (p)
}

format_sample <- function(s) {
  l <- str_split(s, '_')
  l <- sapply(l, function(x) x[length(x)])
  return (l)
}

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

pca <- function(d, center=T, scale=F) {
  # columns are samples
  d <- scale(d, center=center, scale=scale)
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

read_report_meta <- function(filename, n=NULL) {
  if (!is.null(n)) {
  h <- 'cut -f 2-5,7,8,12'
    h <- sprintf('head -n %d %s | %s', n, filename, h)
  } else {
    h <- sprintf('%s %s', h, filename)
  }
  h <- 'cut -f 2-5,7,8,12'
  if (!is.null(n)) {
    h <- sprintf('head -n %d | %s', n, h)
  }
  h <- read.table(pipe(h), head=T, sep='\t')
  names(h) <- tolower(names(h))
  h <- h %>% rename(chromo=chromosome)
  h <- h %>% tbl_df
  return (h)
}

read_report_values <- function(filename, samples=NULL, n=NULL) {
  h <- 'cut -f 13-'
  if (!is.null(n)) {
    h <- sprintf('head -n %d %s | %s', n, filename, h)
  } else {
    h <- sprintf('%s %s', h, filename)
  }

  h <- read.table(pipe(h), head=T, sep='\t')
  if (!is.null(samples)) {
    h <- subset(h, select=intersect(colnames(h), samples))
  }
  h <- h %>% tbl_df
  return (h)
}

read_samples_list <- function(sample_file) {
  d <- read.table(sample_file, head=F) %>% unlist %>% as.vector
  return (d)
}

read_samples_stats <- function(stats_file, samples=NULL) {
  d <- read.table(opts$samples_stats, sep='\t', head=T) %>%
    rename(sample=id, cpg_rate=CpG.rate, chh_rate=CHH.rate, cph_rate=CpH.rate) %>%
    droplevels %>% tbl_df
  if (!is.null(samples)) {
    d <- d[d$sample %in% samples,] %>% droplevels
  }
  return (d)
}

plot_pca_vec <- function(pc_vec, x=1, y=2) {
  t <- data.frame(sample=factor(rownames(pc_vec)),
    pcx=pc_vec[,x], pcy=pc_vec[,y])
  p <- ggplot(t, aes(x=pcx, y=pcy)) + geom_point() +
    geom_text(aes(label=sample), vjust=-.4, hjust= .3, size=3) +
    xlab(sprintf('pc%d', x)) + ylab(sprintf('pc%d', y)) +
    guides(color=F) + theme_pub()
  return (p)
}

plot_pca_val <- function(pc_val) {
  t <- data.frame(pc=1:length(pc_val), val=pc_val)
  p <- ggplot(t, aes(x=pc, y=val)) +
    geom_bar(stat='identity', fill='salmon', color='black') +
    xlab('principle component') +
    ylab('% variance explained') + theme_pub()
  return (p)
}
