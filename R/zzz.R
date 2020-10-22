.onLoad <- function(libname, pkgname)
{
  # =-=-=-=-=-=-
  # Required packages will be listed here
  # =-=-=-=-=-=-
  #requirements = c('dplyr', 'crayon', 'pio',  'ggpubr',
  #                 'cowplot', 'RColorBrewer', 'mobster', 'BMix', 'VIBER')
  #
  #ip = installed.packages()
  #sapply(
  #  requirements,
  #  function(r)
  #  {
  #    if(!(r %in% ip[, 'Package'])) stop("Missing package ", r, " - you should install it to run TINC.")
  #  }
  #)
  #
  #suppressMessages(sapply(requirements, require, character.only = TRUE))

  # =-=-=-=-=-=-
  # Package options
  # =-=-=-=-=-=-
  options(pio.string_fg_colour = crayon::bgYellow$black)

  # =-=-=-=-=-=-
  # Header
  # =-=-=-=-=-=-

  TINC_welcome_message =  getOption('TINC_welcome_message', default = TRUE)

  if(TINC_welcome_message)
  {
    # pio::pioHdr('TINC - Tumour in Normal contamination')
    # pio::pioStr("Author : ", "Giulio Caravagna <gcaravagn@gmail.com>", suffix = '\n')
    # pio::pioStr("GitHub : ", "caravagn/TINC", suffix = '\n')
    # pio::pioStr("   WWW : ", "https://caravagn.github.io/TINC/", suffix = '\n')
    #
    #
    # cat(
    #   "\n > TINC is part of the", crayon::green("\"evoverse\""),
    #   crayon::blue("[https://bit.ly/2orn94e]"),
    #   "- a collection of packages to implement Cancer Evolution analyses from cancer sequencing data.\n"
    # )

    pk = 'TINC'
    pk_l = 'Tumour In Normal Contamination'
    www = "https://caravagn.github.io/TINC/"
    em = "gcaravagn@gmail.com"

    cli::cli_alert_success(
      'Loading {.field {pk}}, {.emph \'{pk_l}\'}. Support : {.url { www}}' )

    # pio::pioStr("GitHub : ", "caravagn/TINC", suffix = '\n')
    # pio::pioStr("   WWW : ", "https://caravagn.github.io/TINC/", suffix = '\n')

    options(TINC_welcome_message = FALSE)
  }

  invisible()
}
