data_select <- function(n, by='met_var') {
  h <- dat$s %>% select_('name', 'id_.x', 'id_.y', sel=by)
  h <- h %>% group_by(name, id_.x) %>% top_n(1, sel) %>% ungroup
  h <- h %>% group_by(name) %>% top_n(n, sel) %>% ungroup
  h <- h %>% select(name, id_.x, id_.y, sel)
  return (h)
}

data_clust <- function(s) {
  em <- dat$em %>% semi_join(s, by=c('name', 'id_.x', 'id_.y')) %>%
    select(id_.x, id_.y, gene_id, sample, expr, met)
  e <- em %>% select(id_.y, gene_id, sample, expr) %>% spread(sample, expr)
  m <- em %>% select(id_.y, gene_id, sample, met) %>% spread(sample, met)
  rs <- dat$s %>% semi_join(s, by=c('name', 'id_.x', 'id_.y'))
  return (list(r=rs, e=e, m=m, em=em))
}

samples_colors <- function(x) {
  h <- rep(opts$clust_colors['default'], length(x))
  for (n in names(opts$clust)) {
    h[x %in% opts$clust[[n]]] <- opts$clust_colors[n]
  }
  names(h) <- x
  return (h)
}

plot_heat <- function(d, xlab='value', col=NULL, col_colors=NULL,
  labRow=NA, ...) {
  d <- as.matrix(d)
  if (is.null(col)) {
    col <- rev(brewer.pal(9, 'Spectral'))
    col <- colorRampPalette(col)(50)
  }

  if (nrow(d) > 500) {
    dendro='column'
  } else {
    dendro = 'both'
  }

  col_colors <- samples_colors(colnames(d))

  p <- heatmap.2(d, density.info='none', trace='none', col=col,
    keysize=0.5, dendro=dendro, labRow=labRow,
    ColSideColors=col_colors,
    lhei=c(2,9),
    lwid=c(2, 5), key.title='', srtCol=45, key.xlab=xlab, ...)
  return (p)
}

plot_heat_expr <- function(dclust, col=NULL, ...) {
  d <- dclust$e %>% select(-id_.y) %>% to_matrix(rowcol='gene_id')
  if (is.null(col)) {
    col <- brewer_cols('OrRd')
  }
  pe <- plot_heat(d, col=col, xlab='Expression', ...)
  return (pe)
}

plot_heat_met <- function(dclust, col=NULL, ...) {
  d <- dclust$m %>% select(-id_.y) %>% to_matrix(rowcol='gene_id')
  if (is.null(col)) {
    col <- brewer_cols('PuBuGn', rev=F)
  }
  pm <- plot_heat(d, col=col, xlab='Methylation', ...)
  return (pm)
}

plot_tracks <- function(d, r=F) {
  te <- plot_track(filter(d, param=='expr_var')) +
    scale_fill_gradient(low='white', 'high'='red', name='var(expr)')
  tm <- plot_track(filter(d, param=='met_var')) +
    scale_fill_gradient(low='white', 'high'='royalblue2', name='var(met)')
  if (r) {
    tr <- plot_track(filter(d, param=='r')) +
      scale_fill_gradient2(low='royalblue', mid='white', 'high'='red2', name='r')
    grid.arrange(tm, tr, te, ncol=3)
  } else {
    grid.arrange(tm, te, ncol=2)
  }
}

brewer_cols <- function(pal='Spectral', rev=F) {
  ncol <- brewer.pal.info[pal, 'maxcolors']
  col <- colorRampPalette(brewer.pal(ncol, pal))(50)
  if (rev) {
    col <- rev(col)
  }
  return (col)
}

to_matrix <- function(d, rowcol=NULL) {
  d <- as.data.frame(d)
  if (!is.null(rowcol)) {
    rownames(d) <- d[[rowcol]]
    d <- d[,setdiff(names(d), rowcol)]
  }
  d <- as.matrix(d)
  return (d)
}

data_track <- function(dclust, p=NA) {
  h <- dclust$e
  if (!is.na(p)) {
    h <- h[p$rowInd,]
  }
  h <- h %>% select(id_.y)
  rc <- h %>% inner_join(dclust$r, by='id_.y')
  d <- rc %>% select(r, expr_var, met_var) %>% mutate(i=1:n()) %>%
    gather(param, value, -i)
  return (d)
}

plot_track <- function(d, low='white', high='red') {
  p <- ggplot(d, aes(x=param, y=i, fill=value)) + geom_tile() +
    theme(line=element_blank(), axis.text=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      panel.background=element_blank(),
      plot.margin=unit(rep(0, 4), 'mm'),
      legend.position='top', legend.direction='vertical')
  return (p)
}
