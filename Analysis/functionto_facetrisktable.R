#' A workaround function in addition to survminer package
#' by Bingxin S.
#' 
#' @examples
#' library(survminer)
#' require("survival")
#' source("https://raw.githubusercontent.com/BingxinS/survminer-fix/master/ggsurvplot_facet_risktable.R")
#'#::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#' data <- lung
#' fit <- survfit(Surv(time, status) ~ sex + ph.ecog, data)
#' ggsurvplot_facet_risktable(fit, data, risktable.ylab="Number at Risk")
#'#::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#' data <- lung
#' fit <- survfit(Surv(time, status) ~  sex + ph.ecog + ph.karno, data)
#' ggsurvplot_facet_risktable(fit, data)
#'#::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#' https://github.com/BingxinS/survminer-fix/blob/master/README_facet_risktable.md


# use the fixed faceted plot along this workaround function, or use the most recent experimental version of survminer
source("https://raw.githubusercontent.com/BingxinS/survminer-fix/master/ggsurvplot_facet_fix_12072018.R")
# report available at https://github.com/BingxinS/survminer-fix/blob/master/README.md

ggsurvplot_facet_risktable <- function (fit, data = NULL, 
                                        risktable.title = NULL,
                                        risktable.ylab = NULL,
                                        text.size = 10,
                                        legend_code = guide_legend(nrow = 2),
                                        ...) 
{
  # get strata.by and facet.by from fit$call
  fit.str <- sub('.*~', '', as.character(fit$call)[2])
  # facet.by = 
  fit.factor <-  strsplit(fit.str, " ")
  fit.factor.length <- length(fit.factor[[1]])/2
  if (fit.factor.length > 3) 
    stop("facet.by has more than 2 factors. Please use 1 or 2 faceting factors.")
  if (fit.factor.length <= 1) 
    stop("facet.by missing. Please use 1 or 2 faceting factors.")
  strata.by <- fit.factor[[1]][2]
  # Generate risk table for each facet plot item
  if(fit.factor.length == 2) {
    facet.by <- fit.factor[[1]][4]
    #================get the plot ================
    plt_fct <- ggsurvplot_facet(fit, data, 
                                    facet.by = facet.by, ...) +
      guides(color = legend_code) + 
      labs(x=element_blank()) 
    #================get the tables================
    ggsurv <- ggsurvplot(fit, data,
                         risk.table = TRUE, 
                         ggtheme = theme_bw())
    facet.formula <- paste0("~", facet.by) %>% stats::as.formula()
    tbl_fct <- 
      ggplot(ggsurv$table$data, ggplot2::aes_string("time", strata.by)) + 
      geom_text(aes(label = n.risk), size = max (2, (text.size-6)) ) +
      facet_wrap(facet.formula, labeller = label_both) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(color="black", size = text.size),
            axis.title = element_text(color="black", size = text.size+4),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.text.y = ggtext::element_markdown(colour = ggpubr::get_palette(palette = "jco", data[,strata.by] %>% levels %>% length(.))
))+
      ggtitle(risktable.title) +
      xlab("Follow-up (Years)") + ylab(risktable.ylab) +
      scale_x_continuous(limits = c(-1, 15))
  }
  
  if(fit.factor.length == 3) {   
    facet.by = fit.factor[[1]][c(4,6)]
    #================get the plot ================
    plt_fct <- ggsurvplot_facet(fit, data, 
                                    facet.by = facet.by, ...) #come back
    #================get the tables================
    ggsurv <- ggsurvplot(fit, data,
                         risk.table = TRUE, 
                         ggtheme = theme_bw(), ...)
    facet.formula <- paste(facet.by, collapse = " ~ ") %>% stats::as.formula()
    tbl_fct <- 
      ggplot(ggsurv$table$data, ggplot2::aes_string("time", strata.by)) + 
      geom_text(aes(label = n.risk), size = max (2, (text.size-6)) ) +
      facet_grid(facet.formula, labeller = label_both) +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill="white", size = 1),
            axis.text.y = element_text(color="black", size = text.size),
            axis.title = element_text(color="black", size = text.size+4),
            strip.text = element_text(color="black", size = text.size))+
      ggtitle(risktable.title) +
      xlab("Time") + ylab(risktable.ylab) 
  }
  #========== combine plot and risk tables ===========
  g_plt_fct <- ggplotGrob(plt_fct)
  g_tbl_fct <- ggplotGrob(tbl_fct)
  min_ncol <- min(ncol(g_plt_fct), ncol(g_tbl_fct))
  g_plt_tbl <- gridExtra::gtable_rbind(g_plt_fct[, 1:min_ncol], g_tbl_fct[, 1:min_ncol], size="last")
  g_plt_tbl$widths <- grid::unit.pmax(g_plt_fct$widths, g_tbl_fct$widths)
  grid::grid.newpage() 
  grid::grid.draw(g_plt_tbl)
}
