plot_raw = function(dataset)
{
  options(warn=-1)

  raw = cowplot::plot_grid(
    plot_tumour_data(dataset),
    plot_normal_data(dataset),
    nrow = 2,
    ncol = 1
  )

  p = cowplot::plot_grid(plot_joint_data(dataset),
                     raw,
                     nrow = 1,
                     ncol = 2)

  options(warn=0)
  p
}

plot_tumour_data = function(dataset)
{
  tumour = as_tumour(dataset$data)
  VAF_range_tumour = dataset$params$VAF_range_tumour

  ggplot(tumour %>%
           filter(VAF > 0),
         aes(VAF, fill = OK_tumour)) +
    geom_vline(
      xintercept = VAF_range_tumour,
      color = 'forestgreen',
      size = .3,
      linetype = 'dashed'
    ) +
    ggrepel::geom_label_repel(
      data = data.frame(x = VAF_range_tumour, label = round(VAF_range_tumour, 2)),
      aes(x = x, y = Inf, label = label),
      fill = 'forestgreen',
      color = 'white',
      size = 3
    ) +
    geom_histogram(binwidth = 0.01) +
    guides(fill = FALSE) +
    xlim(-0.01, 1) +
    labs(y = "n", title = 'Tumour sample') +
    # labs(title = paste0('Tumour (n = ', nrow(tumour %>% filter(VAF > 0)),')')) +
    scale_fill_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray')) +
    my_ggplot_theme() +
    # geom_text(
    #   x = 0.85,
    #   y = h * .8,
    #   label = 'Tumour',
    #   size = 5
    # ) +
    coord_cartesian(clip = 'off')
}

plot_normal_data = function(dataset)
{
  ggplot(
    as_normal(dataset$data),
    aes(VAF)) +
    geom_histogram(binwidth = 0.01, aes(fill = OK_tumour)) +
    xlim(-0.01, 1) +
    my_ggplot_theme() +
    labs(y = "n", title = 'Normal sample') +
    coord_cartesian(clip = 'off') +
    scale_fill_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray')) +
    guides(fill = FALSE)
}

plot_normal_data_inline = function(dataset)
{
  germline = as_normal(dataset$data)

  sorted_data = germline %>% arrange(VAF) %>% mutate(x = row_number())
  x_nonz = which.max(sorted_data$VAF > 0)

  n = sorted_data %>% nrow
  nz = sum(sorted_data$VAF == 0)

  ggplot(sorted_data, aes(x = x, y = VAF)) +
    geom_point(size = .5) +
    # xlim(0, 1) +
    labs(x = 'Mutation') +
    scale_fill_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray'))  +
    my_ggplot_theme() +
    geom_vline(
      xintercept = x_nonz,
      color = 'red',
      size = .3,
      linetype = 'dashed'
    ) +
    ggrepel::geom_label_repel(
      data = data.frame(
        label = paste0('VAF = 0\n(', ((nz / n) * 100) %>% round(2), '%)'),
        x = x_nonz,
        y = Inf
      ),
      size = 3,
      aes(x = x, y = y, label = label)
    )
}

plot_joint_data = function(dataset)
{
  ggplot(
    as_joint(dataset$data) %>%
           mutate(used = ifelse(OK_tumour, "Included", "Excluded")),
         aes(x = n_VAF, y = t_VAF)) +
    geom_point(size = .5, alpha = .5, aes(color = used)) +
    my_ggplot_theme() +
    xlim(0, 1) + ylim(0, 1) + labs(
      title = bquote(bold(TIN) ~ '   Tumour vs normal (n = ' ~ .(sum(dataset$data$OK_tumour)) * ')'),
      y = 'Tumour VAF',
      x = 'Normal VAF'
    ) +
    scale_color_manual(values = c(`Included` = 'black', `Excluded` = 'gray'))  +
    coord_cartesian(clip = 'off') +
    theme(legend.position = c(0.8, 0.2)) +
    guides(color = guide_legend(""))
}
