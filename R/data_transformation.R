as_normal = function(x) {
  x %>%
    dplyr::mutate(DP = n_ref_count + n_alt_count,
                  NV = n_alt_count,
                  VAF = NV / DP) %>%
    dplyr::select(-dplyr::ends_with('count'))
}

as_tumour = function(x) {
  x %>%
    dplyr::mutate(DP = t_ref_count + t_alt_count,
                  NV = t_alt_count,
                  VAF = NV / DP) %>%
    dplyr::select(-dplyr::ends_with('count'))
}

as_joint = function(x) {
  x %>%
    mutate(
      n_VAF = n_alt_count / (n_alt_count + n_ref_count),
      t_VAF = t_alt_count / (t_alt_count + t_ref_count)
    )
}
