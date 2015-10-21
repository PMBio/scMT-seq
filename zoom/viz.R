theme_pub <- function() {
  p <- theme(
    axis.text=element_text(size=rel(1.2), color='black'),
    axis.title=element_text(size=rel(1.5)),
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

plot_scatter <- function(name_, id_.x_, id_.y_, legend=F) {
  d <- dat$em$all %>% filter(name == name_, id_.x == id_.x_, id_.y == id_.y_)
  dm <- d[1,]
  rc <- cmp$r %>% filter(id_.x==id_.x_, id_.y==id_.y_)
  r <- wtd_cor(d$expr, d$met, d$weight)
  pos <- d %>% select(start.y, end.y) %>% distinct
  title <- sprintf('%s (%d-%d) r=%.2f p=%5g',
    name_, rc$start.y, rc$end.y,
    rc$r, rc$p_adj)
  p <- ggplot(d, aes(x=met, y=expr)) +
    stat_smooth(method=lm, color='grey34', aes(weight=weight), size=1) +
    geom_point(aes(color=sample, size=weight), alpha=0.7) +
    scale_size(range=c(4, 8)) +
    xlab('\nMethylation rate') + ylab('Expression rate\n') +
    ggtitle(title)
 if (legend) {
    p <- p +
      guides(
        color=guide_legend(title='Sample ', ncol=7),
        size=F) +
      theme(legend.direction='horizontal')
   } else {
     p <- p + guides(color=F, size=F)
   }
  p <- p + theme_pub()
  return (p)
}


plot_var <- function(d, pa=NULL, xlab=F) {
  p <- ggplot(d, aes(x=0.5*(start.y+end.y), y=wtd_var))
  if (!is.null(pa)) {
    p <- p + pa
  }
  p <- p +
    geom_line(color='blue', lwd=1) +
  # geom_point(color='blue') +
    xlab('') + ylab('Variance\n') +
    scale_x_continuous(labels=comma) +
    theme_pub() +
    theme(axis.title.x=element_blank())
  if (xlab == F) {
    p <- p + theme(axis.text.x=element_blank())
  }
  return (p)
}

plot_cor <- function(d, pa=NULL, legend=F, xlab=F) {
  p <- ggplot(d, aes(x=0.5*(start.y+end.y), y=r)) +
    ylim(-1, 1) + geom_hline(yintercept=0, color='darkgrey')
  if (!is.null(pa)) {
    p <- p + pa
  }
  p <- p +
    geom_ribbon(aes(ymin=r_lo, ymax=r_up), alpha=0.2) +
    geom_line(aes(color=-log10(p)), lwd=1) +
    # geom_point(aes(color=-log10(p)), size=3) +
    scale_color_gradient(low='black', high='red') +
    scale_x_continuous(labels=comma) +
    xlab('') + ylab('Correlation\n') +
    theme_pub() +
    theme(legend.position='top', axis.title.x=element_blank())
  if (xlab == F) {
    p <- p + theme(axis.text.x=element_blank())
  }
  if (!legend) {
    p <- p + guides(color=F)
  }
  return (p)
}

plot_met <- function(d, pa=NULL, legend=F, mean=T) {
  p <- ggplot(d, aes(x=0.5*(start.y+end.y), y=met, color=sample))
  if (!is.null(pa)) {
    p <- p + pa
  }
  if (mean) {
    m <- d %>% group_by(id_.x, id_.y) %>% summarise(
      start.y=unique(start.y), end.y=unique(end.y),
      mean=weighted.mean(met, weight, na.rm=T)
      ) %>% ungroup
    p <- p +
      geom_line(data=m, aes(x=0.5*(start.y+end.y), y=mean),
        color='grey27', size=0.8)
  }

  p <- p +
    geom_point(aes(color=sample, size=log2(weight), alpha=weight)) +
    theme_pub() +
    guides(color=F, alpha=F, size=F) +
    xlab('') + ylab('Methylation rate\n') +
    scale_x_continuous(labels=comma) +
    scale_size(range=c(1, 4)) +
    scale_alpha(range=c(0.5, 1))
  if (legend) {
    p <- p + guides(size=guide_legend(title='# CpGs in window'))
  }
  return (p)
}
# id__ <- 2487
# d <- make_data(id__)
# p <- plot_met(d$em)
# print(p)

make_data <- function(id__, max_dist=opts$max_dist) {
  e <- dat$expr$rate %>% filter(id_ == id__)
  e_meta <- e[1,] %>% select(-c(sample, expr))

  m <- dat$win$met %>% filter(
    chromo == as.character(e_meta$chromo),
    start >= e_meta$start - max_dist,
    end <= e_meta$end + max_dist)
  em <- e %>% inner_join(m, by='sample')

  r <- em %>% group_by(id_.x, id_.y) %>%
    do(wtd_cor(.$met, .$expr, weights=.$weight)) %>% ungroup
  r <- r %>%
    inner_join(dat$expr$meta, by=c('id_.x'='id_')) %>%
    inner_join(dat$win$meta, by=c('id_.y'='id_'))

  v <- em %>% group_by(id_.x, id_.y) %>%
    summarise(
      wtd_mean=weighted.mean(met, weight, na.rm=T),
      wtd_var=wtd.var(met, weight, na.rm=T)
      ) %>% ungroup
  v <- v %>%
    inner_join(dat$expr$meta, by=c('id_.x'='id_')) %>%
    inner_join(dat$win$meta, by=c('id_.y'='id_'))

  stopifnot(nrow(v) == nrow(r))
  stopifnot(all(r$start.x == v$start.x))
  stopifnot(all(r$end.x == v$end.x))

  return (list(em=em, r=r, v=v))
}

data_anno <- function(a) {
  d <- list()
  for (n in names(a)) {
    an <- a[[n]]
    d[[n]]$name <- n
    d[[n]]$start <- an[1]
    d[[n]]$end <- an[2]
  }
  d <- do.call(rbind.data.frame, d) %>% mutate(name=factor(name)) %>%
    gather(pos, x, -name) %>% tbl_df
  return (d)
}

plot_anno <- function(d) {
  p <- geom_vline(data=d, aes(xintercept=x, linetype=name))
  return (p)
}

make_plots_scatter <- function(name_, id__) {
  em <- dat$em$all %>% filter(name == name_, id_.x == id__)
  em_meta <- em %>% select(name, id_.x, id_.y) %>% distinct
  p <- list()
  for (i in 1:nrow(em_meta)) {
    emi <- em_meta[i,]
    p[[i]] <- plot_scatter(em %>% filter(id_.y == emi$id_.y))
  }
  return (p)
}
# name_ <- 'active_enhancer'
# id__ <- 4373
# p <- make_plots_scatter(name_, id__)

make_plots_track <- function(name_, id__) {
  em_meta <- dat$em$meta %>% filter(name == name_, id_.x == id__) %>%
    arrange(start.y)

  p <- list()
  d <- data_anno(em_meta[1,]$start.x, em_meta[1,]$end.x,
    em_meta$start.y, em_meta$end.y)
  p$anno <- plot_anno(d)

  pd <- make_data(id__)
  p$var <- plot_var(pd$v, p$anno)
  p$cor <- plot_cor(pd$r, p$anno)
  p$met <- plot_met(pd$em, p$anno)

  return (p)
}
# name_ <- 'p300'
# id__ <- 2647
# p <- make_plots_track(name_, id__)

make_plots <- function(name_, id__) {
  p <- make_plots_track(name_, id__)
  p$scatter <- make_plots_scatter(name_, id__)
  return (p)
}

plot_to_file <- function(name, id__, out_dir='plots') {
  m <- dat$expr$meta %>% filter(id_ == id__)
  fb <- sprintf('%s/%s_id%d_%s_%s', out_dir, name, id__, m$ens_id, m$gene_id)
  p <- make_plots(name, id__)

  fn <- sprintf('%s_01.pdf', fb)
  nrow <- ceil(length(p$scatter) / 2)
  pdf(file=fn, width=10, height=5 * nrow)
  h <- p$scatter
  h[['ncol']] <- 2
  do.call(grid.arrange, h)
  dev.off()

  fn <- sprintf('%s_02.pdf', fb)
  pdf(file=fn, width=10, height=8)
  h <- list(p$var, p$cor, p$met)
  h[['heights']] <- c(0.2, 0.2, 0.6)
  do.call(grid.arrange, h)
  dev.off()
}
# name_ <- 'active_enhancer'
# id__ <- 4373
# plot_to_file(name_, id__)

plot_tracks <- function(name, id__) {
  p <- make_plots(name, id__)
  grid.arrange(p$var, p$cor, p$met, heights=c(0.2, 0.2, 0.6))
}
# name_ <- 'active_enhancer'
# id__ <- 55
# plot_tracks(name_, id__)

