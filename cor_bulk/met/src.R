read_meta <- function(filename, n=NULL) {
  h <- 'cut -f 2-5,7,8,12'
  if (!is.null(n)) {
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

read_all_meta <- function(filenames) {
  m <- read_meta(filenames[1])
  if (length(filenames) > 1) {
    for (i in 2:length(filenames)) {
      s <- read_meta_quick(filenames[i])
      stopifnot(all(m$start == s$start))
    }
  }
  return (m)
}

read_values <- function(filename, samples=NULL, n=NULL) {
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

read_all_values <- function(filenames, samples=NULL, n=NULL) {
  d <- lapply(filenames, function(x) read_values(x, samples=samples, n=n))
  e <- list()
  for (dd in d) {
    if (ncol(dd) > 0) {
      e[[length(e) + 1]] <- dd
    }
  }
  d <- e
  stopifnot(length(d) > 0)
  h <- d[[1]]
  if (length(d) > 1) {
    for (i in 2:length(d)) {
      if (ncol(d[[i]]) > 0) {
        h <- cbind.data.frame(h, d[[i]])
      }
    }
  }
  h <- h %>% tbl_df
  return (h)
}
