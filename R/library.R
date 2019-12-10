# R_LIB = "~/re_gecip/shared_allGeCIPs/R_dm/R_library/"
#
# .libPaths(R_LIB)
# library(tidyverse, lib.loc = R_LIB)
# library(mobster, lib.loc = R_LIB)
# library(BMix, lib.loc = R_LIB)
# library(cowplot, lib.loc = R_LIB)
# library(pio, lib.loc = R_LIB)
# library(ggpubr, lib.loc = R_LIB)
# library(ggforce, lib.loc = R_LIB)


# Plot a figure with layout 1x3 reporting raw data
# plot_raw = function(dataset,
#                     label,
#                     percentage = .3,
#                     VAF_range_tumour = c(0, 0.7))
# {
#   S =  percentage * nrow(dataset$tumour)
#
#   tumour = dataset$tumour %>% mutate(PASS = (filters == 'PASS'))
#   germline = dataset$germline %>% mutate(PASS = (filters == 'PASS'))
#   joint = dataset$joint %>% sample_n(S) %>% mutate(PASS = (filters == 'PASS'))
#
#   f1 = ggplot(tumour %>%
#                 filter(VAF > 0, PASS) %>%
#                 mutate(
#                   fill = ifelse(
#                     VAF > VAF_range_tumour[1] &  VAF < VAF_range_tumour[2],
#                   "Included",
#                   "Excluded")
#                   ),
#               aes(VAF, fill = fill)) +
#     geom_vline(xintercept = VAF_range_tumour, color = 'forestgreen', size = .3, linetype = 'dashed') +
#     geom_histogram(binwidth = 0.01) +
#     guides(fill = guide_legend('')) +
#     xlim(0, 1) +
#     labs(title = 'Tumour / Normal') +
#     scale_fill_manual(values = c(`Included` = 'black', `Excluded` = 'gray')) +
#     # facet_wrap(~ PASS, nrow = 2, scales = 'free_y') +
#     mobster:::my_ggplot_theme()
#
#   f2 = ggplot(germline  %>% filter(VAF > 0, PASS), aes(VAF)) +
#     geom_histogram(binwidth = 0.01, fill = 'black') +
#     xlim(0, 1) +
#     scale_fill_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray'))  +
#     mobster:::my_ggplot_theme()
#
#   f12 = ggpubr::ggarrange(f1,f2, legend = 'bottom', common.legend = T, nrow = 2, ncol = 1)
#
#   f3 = ggplot(joint %>% filter(PASS),
#               aes(x = VAF.germline, y = VAF.tumour)) +
#     #   stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
#     #   scale_fill_viridis_c()
#     geom_point(size = .5, alpha = .5) +
#     mobster:::my_ggplot_theme() +
#     xlim(0, 1) + ylim(0, 1) + labs(
#       title = paste0('Tumour versus germline (n = ', nrow(joint %>% filter(PASS)), ')'),
#       y = 'Tumour VAF',
#       x = 'Normal VAF',
#       subtitle = label,
#       caption = paste0('S = ', S, ' (', percentage * 100, '% of n)')
#     ) +
#     scale_color_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray'))
#
#
#   pl = cowplot::plot_grid(f1,f2, nrow = 2, ncol = 1, align = 'v')
#   pl = cowplot::plot_grid(f3,pl, nrow = 1, ncol = 2, align = 'h')
#
#
#   # cowplot::plot_grid(f1,
#   #                    f2,
#   #                    f3,
#   #                    nrow = 1,
#   #                    ncol = 3,
#   #                    align = 'h')
#   # f1 = ggplot(tumour %>% filter(VAF > 0), aes(VAF, fill = PASS)) +
#   #   geom_rect(
#   #     xmin = VAF_range_tumour[1],
#   #     xmax = VAF_range_tumour[2],
#   #     ymin = -Inf,
#   #     ymax = Inf,
#   #     fill = alpha("forestgreen", 0.2)
#   #   ) +
#   #   geom_histogram(binwidth = 0.01) +
#   #   xlim(0, 1) +
#   #   labs(title = 'Tumour sample', subtitle = label) +
#   #   scale_fill_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray')) +
#   #   facet_wrap(~ PASS, nrow = 2, scales = 'free_y') +
#   #   mobster:::my_ggplot_theme()
#   #
#   # f2 = ggplot(germline  %>% filter(VAF > 0), aes(VAF, fill = PASS)) +
#   #   geom_histogram(binwidth = 0.01) + xlim(0, 1) +
#   #   labs(title = 'Germline sample', subtitle = label) +
#   #   scale_fill_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray'))  +
#   #   facet_wrap(~ PASS, nrow = 2, scales = 'free_y') +
#   #   mobster:::my_ggplot_theme()
#   #
#   #
#   # f3 = ggplot(joint,
#   #             aes(x = VAF.germline, y = VAF.tumour, color = PASS)) +
#   #   #   stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
#   #   #   scale_fill_viridis_c()
#   #   geom_point(size = .5, alpha = .1) +
#   #   mobster:::my_ggplot_theme() +
#   #   xlim(0, 1) + ylim(0, 1) + labs(
#   #     title = 'Tumour versus germline',
#   #     y = 'VAF tumour',
#   #     x = 'VAF germline',
#   #     subtitle = label,
#   #     caption = paste0('S = ', S, '(', percentage * 100, '%)')
#   #   ) +
#   #   scale_color_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray'))  +
#   #   facet_wrap(~ PASS, nrow = 2)
#   #
#   #
#   # cowplot::plot_grid(f1,
#   #                    f2,
#   #                    f3,
#   #                    nrow = 1,
#   #                    ncol = 3,
#   #                    align = 'h')
# }

#
# plot_raw = function(dataset,
#                     label,
#                     VAF_range_tumour = c(0, 0.7))
# {
#   tumour = dataset$tumour %>% filter(filters == 'PASS')
#   germline = dataset$germline %>% filter(filters == 'PASS')
#   joint = dataset$joint %>% filter(filters == 'PASS')
#
#   h = hist(tumour%>%
#              filter(VAF > 0) %>%
#              pull(VAF), plot = FALSE, breaks = seq(0, 1, 0.01))$counts %>% max
#
#   f1 = ggplot(tumour %>%
#                 filter(VAF > 0) %>%
#                 mutate(
#                   fill = ifelse(
#                     VAF > VAF_range_tumour[1] &  VAF < VAF_range_tumour[2],
#                     "Included",
#                     "Excluded")
#                 ),
#               aes(VAF, fill = fill)) +
#     geom_vline(xintercept = VAF_range_tumour, color = 'forestgreen', size = .3, linetype = 'dashed') +
#     ggrepel::geom_label_repel(
#       data = data.frame(x=VAF_range_tumour, label = round(VAF_range_tumour, 2)),
#       aes(x = x, y = Inf, label = label),
#       fill = 'forestgreen', color = 'white', size = 3
#     ) +
#     geom_histogram(binwidth = 0.01) +
#     guides(fill = FALSE) +
#     xlim(-0.01, 1) +
#     labs(y = "n", title =' ') +
#     # labs(title = paste0('Tumour (n = ', nrow(tumour %>% filter(VAF > 0)),')')) +
#     scale_fill_manual(values = c(`Included` = 'black', `Excluded` = 'gray')) +
#     # facet_wrap(~ PASS, nrow = 2, scales = 'free_y') +
#     mobster:::my_ggplot_theme() +
#     geom_text(x = 0.85, y= h * .8, label = 'Tumour', size = 5) +
#     coord_cartesian(clip = 'off')
#
#
#   h = hist(germline$VAF, plot = FALSE, breaks = seq(0, 1, 0.01))$counts %>% max
#
#   f2 = ggplot(germline, aes(VAF)) +
#     geom_histogram(binwidth = 0.01, fill = 'black') +
#     xlim(-0.01, 1) +
#     scale_fill_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray'))  +
#     mobster:::my_ggplot_theme()+
#     geom_text(x = 0.85, y= h * .8, label = 'Normal', size = 5) +
#     labs(y = "n") +
#     coord_cartesian(clip = 'off')
#     # labs(title = paste0('Normal (n = ', nrow(germline),')'))
#
#
#   sorted_data = germline %>% arrange(VAF) %>% mutate(x=row_number())
#   x_nonz = which.max(sorted_data$VAF > 0)
#   q_x = quantile(sorted_data$VAF, probs = c(0.5))
#
#   f2b = ggplot(sorted_data, aes(x=x,y=VAF)) +
#     geom_point(size = .5) +
#     # xlim(0, 1) +
#     scale_fill_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray'))  +
#     mobster:::my_ggplot_theme() +
#     geom_vline(xintercept = x_nonz, color = 'red', size = .3, linetype = 'dashed') +
#     ggrepel::geom_label_repel(
#       data = data.frame(
#         label = paste('VAF = 0'),
#         x = x_nonz,
#         y = Inf),
#       size = 3,
#       aes(x=x,y=y,label=label)
#     )
#
#
#   # f12 = ggpubr::ggarrange(f1,f2, legend = 'bottom', common.legend = T, nrow = 2, ncol = 1)
#   # cowplot::plot_grid(f12,f2b, align = 'v',  nrow = 2, ncol = 1, rel_heights = c(.75, 0.25))
#   f12 = cowplot::plot_grid(f1,f2, nrow = 2, ncol = 1)
#
#
#
#   f3 = ggplot(joint,
#               aes(x = VAF.germline, y = VAF.tumour)) +
#     #   stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
#     #   scale_fill_viridis_c()
#     geom_point(size = .5, alpha = .5) +
#     mobster:::my_ggplot_theme() +
#     xlim(0, 1) + ylim(0, 1) + labs(
#       title = paste0('Tumour vs germline (n = ', nrow(joint), ')'),
#       y = 'Tumour VAF',
#       x = 'Normal VAF'
#       # caption = paste0('S = ', S, ' (', percentage * 100, '% of n)')
#     ) +
#     scale_color_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray'))  +
#     # geom_text(x = 0.65, y= .95,
#     #           label = paste0('Tumour versus germline (n = ', nrow(joint), ')'),
#     #           size = 5) +
#     coord_cartesian(clip = 'off')
#
#
#
#   # pl = cowplot::plot_grid(f12,f2, nrow = 2, ncol = 1, align = 'v')
#   pl = cowplot::plot_grid(f3,f12, nrow = 1, ncol = 2)
#
#   return(pl)
#
#   # cowplot::plot_grid(f1,
#   #                    f2,
#   #                    f3,
#   #                    nrow = 1,
#   #                    ncol = 3,
#   #                    align = 'h')
#   # f1 = ggplot(tumour %>% filter(VAF > 0), aes(VAF, fill = PASS)) +
#   #   geom_rect(
#   #     xmin = VAF_range_tumour[1],
#   #     xmax = VAF_range_tumour[2],
#   #     ymin = -Inf,
#   #     ymax = Inf,
#   #     fill = alpha("forestgreen", 0.2)
#   #   ) +
#   #   geom_histogram(binwidth = 0.01) +
#   #   xlim(0, 1) +
#   #   labs(title = 'Tumour sample', subtitle = label) +
#   #   scale_fill_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray')) +
#   #   facet_wrap(~ PASS, nrow = 2, scales = 'free_y') +
#   #   mobster:::my_ggplot_theme()
#   #
#   # f2 = ggplot(germline  %>% filter(VAF > 0), aes(VAF, fill = PASS)) +
#   #   geom_histogram(binwidth = 0.01) + xlim(0, 1) +
#   #   labs(title = 'Germline sample', subtitle = label) +
#   #   scale_fill_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray'))  +
#   #   facet_wrap(~ PASS, nrow = 2, scales = 'free_y') +
#   #   mobster:::my_ggplot_theme()
#   #
#   #
#   # f3 = ggplot(joint,
#   #             aes(x = VAF.germline, y = VAF.tumour, color = PASS)) +
#   #   #   stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
#   #   #   scale_fill_viridis_c()
#   #   geom_point(size = .5, alpha = .1) +
#   #   mobster:::my_ggplot_theme() +
#   #   xlim(0, 1) + ylim(0, 1) + labs(
#   #     title = 'Tumour versus germline',
#   #     y = 'VAF tumour',
#   #     x = 'VAF germline',
#   #     subtitle = label,
#   #     caption = paste0('S = ', S, '(', percentage * 100, '%)')
#   #   ) +
#   #   scale_color_manual(values = c(`TRUE` = 'black', `FALSE` = 'gray'))  +
#   #   facet_wrap(~ PASS, nrow = 2)
#   #
#   #
#   # cowplot::plot_grid(f1,
#   #                    f2,
#   #                    f3,
#   #                    nrow = 1,
#   #                    ncol = 3,
#   #                    align = 'h')
# }

# Analyze tumour data (x is input)
# - carries out a fit with up to 3 mixture components fit by the BIC score
# - determines the clonal cluster
# -
# analyse_MOBSTER = function(x,
#                            cutoff_miscalled_clonal,
#                            cutoff_lv_assignment,
#                            ...)
# {
#   #  MOBSTER fit ad plot
#   mobster_fit_tumour = mobster_fit(x,
#                                    ...)
#
#   # Clonal cluster is guessed to be the highest peak below cutoff_miscalled_clonal = 510% VAF - seems reasonable without CNA
#   all_clusters = mobster_fit_tumour$best$Clusters %>%
#     filter(type == 'Mean', fit.value < cutoff_miscalled_clonal) %>%
#     arrange(desc(fit.value)) %>%
#     pull(cluster)
#
#   # We also check that the clonal cluster has higher dimension of the others
#   cl_map = sapply(all_clusters, mobster:::is_reasonable_clonal_cluster,
#          x = mobster_fit_tumour$best)
#
#   clonal_cluster_fit = NULL
#   if(all(!cl_map)) {
#     warning("Check clonal cluster")
#     clonal_cluster_fit = all_clusters[1]
#   }
#   else{
#     if(!cl_map[1]) message("Not using C1 as clonal cluster")
#     clonal_cluster_fit = names(which.max(cl_map))
#   }
#
#   # # top-2 if Betas
#   # if (length(all_clusters) > 1)
#   # {
#   #   top2 = all_clusters[1:2]
#   #
#   #   if (top2[2] != 'Tail')
#   #   {
#   #     m_top2 = mobster_fit_tumour$best$Clusters %>%
#   #       filter(type == 'Mean', fit.value < cutoff_miscalled_clonal) %>%
#   #       arrange(desc(fit.value))  %>%
#   #       pull(fit.value)
#   #
#   #     d_m1 = mobster::ddbpmm(mobster_fit_tumour$best, data = m_top2[1])
#   #     d_m2 = mobster::ddbpmm(mobster_fit_tumour$best, data = m_top2[2])
#   #
#   #     if (d_m2 < d_m1)
#   #       clonal_cluster_fit = top2
#   #   }
#   # }
#   #
#   # clonal_cluster_fit = estimated_tumour_purity = clonal_tumour = NULL
#   # m = 1
#   #
#   # repeat
#   # {
#   #   # Candidate clonal cluster
#   #   clonal_cluster_fit = all_clusters[1:m]
#
#   # Tumour purity is therefore 2 * the Beta peak of the clonal cluster
#   estimated_tumour_purity = mobster_fit_tumour$best$Clusters %>%
#     filter(cluster %in% clonal_cluster_fit, type == 'Mean') %>%
#     pull(fit.value) %>% mean * 2
#
#   # List of clonal mutations in the tumour, with LV > cutoff_lv_assignment
#   clonal_tumour = mobster::Clusters(mobster_fit_tumour$best, cutoff_assignment = cutoff_lv_assignment) %>%
#     filter(cluster %in% clonal_cluster_fit) %>%
#     pull(id)
#
#   # # Stop with > 5 clonal mutations
#   # if(length(clonal_tumour) > 5) break;
#   #
#   # Stop if we tried all
#   #   if(m == length(all_clusters)) stop("No clonal mutations with this LV value, try lowering.")
#   #
#   #   message("Augmenting clonal cluster to achieve more mutations")
#   #   m = m + 1
#   # }
#
#   # plot_input = ggplot(x, aes(VAF)) + geom_histogram(binwidth = 0.01) + xlim(0, 1)
#   plot_fit_tumour = plot.dbpmm(mobster_fit_tumour$best, cutoff_assignment = cutoff_lv_assignment) +
#     labs(title = paste0('MOBSTER fit (clonal: ', clonal_cluster_fit, ')'), caption = NULL)
#
#   plot_lv = mobster::plot_latent_variables(mobster_fit_tumour$best, cutoff_lv_assignment)
#   # plot_fs = mobster::plot_entropy(mobster_fit_tumour$best)
#
#   figure = cowplot::plot_grid(
#     # plot_input,
#     plot_fit_tumour,
#     plot_lv,
#     # plot_fs,
#     ncol = 2,
#     nrow = 1,
#     align = 'h',
#     axis = 'bt'
#   )
#
#
#
#   return(
#     list(
#       input = mobster::Clusters(mobster_fit_tumour$best) %>% rename(tumour.cluster = cluster),
#       fit = mobster_fit_tumour,
#       plot = figure,
#       clonal_cluster = clonal_cluster_fit,
#       clonal_mutations = clonal_tumour,
#       estimated_purity = estimated_tumour_purity
#     )
#   )
# }

# Analyse germline samples with BMix
# analyse_germline = function(x)
# {
#   df_input = x %>%
#     select(NV, DP) %>%
#     as.data.frame
#
#
#   fit_normal = bmixfit(df_input,
#                        K.BetaBinomials = 0)
#
#   Binomial_peaks = BMix::Parameters(fit_normal) %>% pull(mean)
#   Binomial_pi = fit_normal$pi
#
#   clonal = Binomial_peaks %*% Binomial_pi
#   nzero = sum(df_input$NV == 0)
#
#   #  Plots
#   # plot_input = ggplot(x, aes(VAF)) + geom_histogram(binwidth = 0.01) + xlim(0, 1)
#   plot_fit_normal = BMix::plot_clusters(fit_normal, data = df_input) +
#     xlim(-0.01, 1) +
#     labs(
#       title = paste0("BMix fit - n = ", nrow(df_input), " (", nzero, ' zeros)'),
#       subtitle = paste0(names(fit_normal$pi), ' (', round(fit_normal$pi, 2), '%)', collapse = ', ')
#       ) +
#     geom_vline(xintercept = Parameters(fit_normal) %>% pull(mean), linetype = 'dashed', size = .3) +
#     ggrepel::geom_label_repel(
#       data = BMix::Parameters(fit_normal),
#       aes(x = mean, y = Inf, label = round(mean, 2)),
#       fill = 'black', color = 'white', size = 3
#     ) +
#     scale_fill_brewer(palette = "Set2")
#
#     # xlim(-0.01, min(1, max(df_input$NV/df_input$DP) * 2.5))
#
#   # plot_fit_normal_d = BMix::plot_density(fit_normal, data = df_input)
#   plot_fit_normal_d = bin_heatmap(C = median(df_input$DP), p = round(clonal, 2))
#
#   figure = cowplot::plot_grid(
#     # plot_input,
#     plot_fit_normal,
#     plot_fit_normal_d,
#     ncol = 2,
#     nrow = 1,
#     align = 'h'
#   )
#
#   # TIN estimation
#
#   clonal_cluster = names(fit_normal$pi)
#   ccl =   paste0(round(fit_normal$pi, 2), ' x ', round(fit_normal$B.params, 2), collapse = ' + ')
#   # clonal_cluster = names(which.max(fit_normal$B.params))
#   # highest_Binomial_peak = fit_normal$B.params[clonal_cluster]
#   highest_Binomial_peak = as.vector(clonal)
#
#
#   x$germline.cluster = fit_normal$labels
#
#   # estimated_normal_purity = 1 - (0.5 - highest_Binomial_peak) * 2
#   estimated_normal_purity =  highest_Binomial_peak * 2
#
#   return(
#     list(
#       input = x,
#       fit = fit_normal,
#       plot = figure,
#       clonal_cluster = clonal_cluster,
#       clonal_mutations = x %>% filter(germline.cluster == clonal_cluster) %>% pull(id),
#       estimated_purity = estimated_normal_purity
#     )
#   )
#
# # }
#
# # Assemble results from tumour and germlinegermline
# assemble_results = function(joint,
#                             MOBSTER_analysis_results,
#                             germline_analysis_results)
# {
#   data_table = germline_analysis_results$input %>%
#     full_join(
#       MOBSTER_analysis_results$input,
#       by = c("chr", "from", "to", "ref", "alt", "id", "filters"),
#       suffix = c('.germline', '.tumour')
#     ) %>%
#     mutate(
#       is_clonal = (
#         tumour.cluster %in% MOBSTER_analysis_results$clonal_cluster
#       ) & (
#         germline.cluster %in% germline_analysis_results$clonal_cluster
#       ),
#       VAF.germline = ifelse(is.na(germline.cluster), 0, VAF.germline)
#     ) %>%
#     select(
#       chr,
#       from,
#       to,
#       ref,
#       alt,
#       id,
#       filters,
#       ends_with('count'),
#       ends_with('tumour'),
#       ends_with('germline'),
#       ends_with('cluster'),
#       is_clonal
#     )
#
#   TIT = MOBSTER_analysis_results$estimated_purity / 2
#   TIN = germline_analysis_results$estimated_purity / 2
#
#   # S =  percentage * nrow(joint)
#
#   ggplot(data_table %>% mutate(is_clonal = ifelse(is_clonal, 'Clonal', 'Subclonal')),
#                   aes(x = VAF.germline, y = VAF.tumour, color = is_clonal)) +
#     geom_point(
#       size = .3,
#       alpha = 0.1,
#       color = 'black'
#     ) +
#     # scale_x_continuous(limits=c(0, 1)) +
#     facet_zoom(xlim = c(0, 0.05)) +
#     ylim(0, 1) +
#     mobster:::my_ggplot_theme()
#
#
#   figure = ggplot(data_table %>% mutate(is_clonal = ifelse(is_clonal, 'Clonal', 'Subclonal or CNA')),
#                   aes(x = VAF.germline, y = VAF.tumour, color = is_clonal)) +
#     geom_point(
#       data = joint,
#       size = .3,
#       alpha = 0.1,
#       color = 'black'
#     ) +
#     geom_point(size = .7, alpha = .3) +
#     geom_abline(linetype = 'dashed',
#                 color = 'red',
#                 size = .3) +
#     geom_hline(
#       yintercept = TIT,
#       linetype = 'dashed',
#       color = 'black',
#       size = .3
#     ) +
#     geom_vline(
#       xintercept = TIN,
#       linetype = 'dashed',
#       color = 'black',
#       size = .3
#     ) +
#     geom_point(
#       data = data.frame(x = TIN, y = TIT),
#       aes(x, y),
#       color = 'black',
#       fill = 'white',
#       size = 3,
#       shape = 21
#     ) +
#     xlim(0, 1) +
#     ylim(0, 1) +
#     guides(color = guide_legend('')) +
#     mobster:::my_ggplot_theme() +
#     labs(
#       title = paste0("Contamination"),
#       x = 'VAF germline',
#       y = "VAF tumour"
#       # subtitle = paste0(
#       #   "Cancer cells in normal sample: ",
#       #   round(TIN * 200, 3),
#       #   '%\n',
#       #   "Cancer cells in tumour sample: ",
#       #   round(TIT * 200, 3),
#       #   '%'
#       # )
#       # caption = paste0('S = ', S, '(', percentage * 100, '%)')
#     ) +
#     scale_color_manual(values = c(`Clonal` = 'cyan4', `Subclonal or CNA` = 'black'))
#
#     # geom_mark_ellipse(aes(fill = is_clonal, label = is_clonal))
#
#
#   data_table = data_table %>%
#     mutate(
#       group = paste(tumour.cluster, '~', germline.cluster)
#       )
#
#   group_labels = data_table$group %>% unique
#   group_labels_colors = pio:::nmfy(group_labels, rep("black", length(group_labels)))
#
#   for(cl in MOBSTER_analysis_results$clonal_cluster)
#     group_labels_colors[grepl(group_labels, pattern = cl)] = 'cyan4'
#
#
#   f2 = ggplot(data_table,
#               aes(
#                 x = VAF.germline,
#                 y = VAF.tumour,
#                 shape = group,
#                 color = group
#               )) +
#     # geom_point(
#     #   data = joint,
#     #   aes(x = VAF.germline,
#     #       y = VAF.tumour),
#     #   size = .3,
#     #   alpha = 0.1,
#     #   color = 'black',
#     #   inherit.aes = FALSE
#     #
#     # ) +
#     # geom_density_2d(data = data_table %>% filter(VAF.germline > 0), size = .1) +
#     geom_point(size = 2, alpha = .7) +
#     geom_abline(linetype = 'dashed',
#                 color = 'red',
#                 size = .3) +
#     geom_hline(
#       yintercept = TIT,
#       linetype = 'dashed',
#       color = 'black',
#       size = .3
#     ) +
#     geom_vline(
#       xintercept = TIN,
#       linetype = 'dashed',
#       color = 'black',
#       size = .3
#     ) +
#     geom_point(
#       data = data.frame(x = TIN, y = TIT),
#       aes(x, y),
#       color = 'black',
#       fill = 'white',
#       size = 3,
#       shape = 21
#     ) +
#     guides(shape = guide_legend("Clusters"), color = FALSE) +
#     mobster:::my_ggplot_theme() +
#     labs(
#       title = paste0(""),
#       x = 'VAF germline',
#       y = "VAF tumour"
#       # subtitle = paste0(
#       #   "Clonal cluster (tumour) : ",
#       #   paste(MOBSTER_analysis_results$clonal_cluster, collapse = ', '),
#       #   " (n = ", length(MOBSTER_analysis_results$clonal_mutations), ')',
#       #   '\n',
#       #   "Germline weigth : ",
#       #   paste(germline_analysis_results$clonal_cluster, collapse = ', ')
#       # )
#     ) +
#     scale_color_manual(values = group_labels_colors)
#
#   vtin = data.frame(
#     variable = c('Tumour cells', 'Normal cells', 'Tumour cells', 'Normal cells'),
#     value = c(TIN * 2, 1 - (TIN * 2), TIT * 2, 1 - (TIT * 2)),
#     sample = c(
#       'Normal sample',
#       'Normal sample',
#       'Tumour sample',
#       'Tumour sample'
#     )
#   )
#
#   f3 = ggplot(data = vtin,
#               aes(x='1', y = value, fill = variable)) +
#     geom_bar(stat = 'identity') +
#     mobster:::my_ggplot_theme() +
#     facet_wrap(~sample, ncol = 2, nrow = 1) +
#     # coord_flip() +
#     coord_polar("y", start = 0, clip = 'off') +
#     labs(
#       title = paste0("Sample composition"),
#       x = element_blank(),
#       y = "Percentage",
#       subtitle = paste0(
#         "TIN: ",
#         round(vtin$value[1] * 100, 2),
#         '%',
#         " (",
#         round(vtin$value[2] * 100, 2),
#         '% germline)\n',
#         "TIT : ",
#         round(vtin$value[3] * 100, 2),
#         '% (',
#         round(vtin$value[4] * 100, 2),
#         '% germline)'
#       )
#     ) +
#     scale_fill_manual(values = c(`Tumour cells` = 'plum4', `Normal cells` = 'plum2')) +
#     guides(fill = guide_legend(''))
#
#
#   return(list(data_table = data_table,
#               summary_plot = figure,
#               zoom_plot = f2,
#               sample_composition = f3,
#               TIN = TIN, TIT = TIT))
#
#   #
#   # figure = cowplot::plot_grid(
#   #   # plot_input,
#   #   figure,
#   #   f2,
#   #   f3,
#   #   ncol = 3,
#   #   nrow = 1,
#   #   align = 'h', axis = 'x'
#   # )
#   #
#   #
#   # return(list(data_table = data_table, plot = figure, TIN = TIN, TIT = TIT))
#   #
# }



