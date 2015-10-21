plot_heat <- function(d, xlab='value', col_colors=NULL, ...) {
  d <- as.matrix(d)
  colors <- rev(brewer.pal(9, 'Spectral'))
  colors <- colorRampPalette(colors)(50)

  if (nrow(d) > 500) {
    dendro='column'
  } else {
    dendro = 'both'
  }
  p <- heatmap.2(d, density.info='none', trace='none', col=colors,
    keysize=1.0, dendro=dendro,
    lwid=c(2, 5), key.title='', srtCol=45, key.xlab=xlab, ...)
  return (p)
}

get_from_sample <- function(x, what='cond') {
  d <- dat$samples
  h <- d[match(x, as.vector(d$sample)),][[what]]
  return (h)
}

plot_heat_col <- function(d) {
  cond_colors <- c('2i'='green', 'serum'='yellow')
  cond <- as.vector(get_from_sample(colnames(d), 'cond'))
  col_colors <- cond_colors[cond]
  p <- plot_heat(d, col_colors=col_colors)
  return (p)
}


## Clustering

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
  p <- heatmap.2(d, density.info='none', trace='none', col=col,
    keysize=1.0, dendro=dendro, labRow=labRow,
    lwid=c(2, 5), key.title='', srtCol=45, key.xlab=xlab, ...)
  return (p)
}

# select n genes from context with strongest correlation
data_clust <- function(name_='prom_non_cgi', n=500, genes=F) {
  rs <- dat$rs %>% filter(name == name_) %>% group_by(name, id_.x) %>%
    top_n(1, abs(r)) %>% ungroup
  if (genes) {
    rs <- rs %>% semi_join(select(dat$genes, ens_id), by=c('ens_id'))
  }
  rs <- rs %>% group_by(name) %>% top_n(n, abs(r)) %>% ungroup
  em <- dat$em %>% semi_join(rs, by=c('name', 'id_.x', 'id_.y')) %>%
    select(id_.x, id_.y, gene_id, sample, expr, met)
  e <- em %>% select(id_.y, gene_id, sample, expr) %>% spread(sample, expr)
  m <- em %>% select(id_.y, gene_id, sample, met) %>% spread(sample, met)
  return (list(r=rs, e=e, m=m, em=em))
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

plot_tracks <- function(d) {
  tr <- plot_track(filter(d, param=='r')) +
    scale_fill_gradient2(low='royalblue', mid='white', 'high'='red2', name='r')
  te <- plot_track(filter(d, param=='expr_var')) +
    scale_fill_gradient(low='white', 'high'='red', name='var(expr)')
  tm <- plot_track(filter(d, param=='met_var')) +
    scale_fill_gradient(low='white', 'high'='royalblue2', name='var(met)')
  grid.arrange(te, tr, tm, ncol=3)
}

plot_heat_expr <- function(dclust, col=NULL) {
  d <- dclust$e %>% select(-id_.y) %>% to_matrix(rowcol='gene_id')
  if (is.null(col)) {
    col <- brewer_cols('OrRd')
  }
  pe <- plot_heat(d, col=col, xlab='Expression')
  return (pe)
}

plot_heat_met <- function(dclust, pe, col=NULL) {
  d <- dclust$m %>% select(-id_.y) %>% to_matrix(rowcol='gene_id')
  if (is.null(col)) {
    col <- brewer_cols('PuBuGn', rev=F)
  }
  pm <- plot_heat(d, col=col, Rowv=pe$RowDendrogram, xlab='Methylation')
  return (pm)
}

