#' Correlation Hex Plot Function
#'
#' This function generates correlation hexagon plots for a given data frame.
#'
#' @importFrom stats cor
#' @importFrom dplyr any_of filter mutate select distinct group_split left_join
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_hex scale_fill_viridis_c theme_bw theme element_rect margin element_text element_line element_blank labs facet_wrap geom_text
#' @importFrom patchwork plot_layout wrap_plots
#' @importFrom ggtext element_textbox
#' @importFrom reshape2 melt
#' @importFrom purrr map map_chr walk2
#' @importFrom utils globalVariables
#'
#' @param dt A data frame containing the data to be analyzed.
#' @param id.col The name of the ID column, default is "protein".
#' @param cor.method The correlation method to use, default is "pearson".
#' @param savefile The directory to save plots, default is "outputfile".
#' @param singleplotsize Size of each individual plot as a numeric vector of width and height, default is c(3, 2.5).
#' @param facetplotsize Size of facet plot as a numeric vector of width and height, default is c(9, 7.5).
#' @param bin Number of bins for hex plot, default is 50.
#'
#' @return A list containing facet plot and individual plots.
#' @export
#'
cor_hex=function(dt=dtplot,
                 id.col="protein",
                 cor.method=c("pearson", "kendall", "spearman")[1],
                 savefile="outputfile",
                 singleplotsize=c(3,2.5),#width height
                 facetplotsize=c(3*3,2.5*3),#width height
                 bin=50
){



  dt=dt |> dplyr::select(any_of(id.col),dplyr::everything())
  if (!is.data.frame(dt)) {
    dt <- as.data.frame(dt)
  }

  dt[-1] <- lapply(dt[-1], function(x) {
    if (is.numeric(x)) {
      x[is.infinite(x) | x == 0 | x == "" | is.nan(x)] <- NA
    }
    return(x)
  })


  dt[-1] <- lapply(dt[-1], function(x) as.numeric(as.character(x)))

  corvalue<- round(cor(as.matrix(dt[-1]), method = cor.method,use = "pairwise.complete.obs"), 3)



  corInofrs=reshape2::melt(corvalue, na.rm = TRUE) |> #dplyr::rename("xlab"="Var1","ylab"="Var2") |>
    dplyr::mutate(compars=paste0(Var1,"_vs_",Var2))



  dt_filtered=dt |>
    tidyr::pivot_longer(!id.col, names_to = "xlab", values_to = "x") |>
    dplyr::left_join(dt |> tidyr::pivot_longer(!id.col, names_to = "ylab", values_to = "y"),by = id.col)  |>
    dplyr::filter(xlab > ylab) |>
    dplyr::mutate(compars=paste0(xlab,"_vs_",ylab),
                  x=log10(x+1),
                  y=log10(y+1))  |>
    dplyr::left_join(corInofrs, by = "compars")


  #  dt_filtered <- dt.plot1  |> dplyr::filter(xlab > ylab)

  dt_label <- dt_filtered  |>
    dplyr::distinct(compars, .keep_all = TRUE)

  plt=ggplot(dt_filtered, aes(x = x, y = y)) +
    geom_hex(bins = bin) +
    scale_fill_viridis_c(option = "viridis", guide = "none") +
    theme_bw() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "white"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      axis.title.x = element_text(size = 8, face = 'bold', margin = margin(t = 0.5)),
      axis.title.y = element_text(size = 8, face = 'bold', margin = margin(r = 0.5)),
      axis.text = element_text(size = 8, color = "black"),
      axis.line.x = element_line(color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 8, face = "italic"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right",
      strip.text = ggtext::element_textbox(
        size = 8, face = 'italic', color = "grey20",
        hjust = 0.5, halign = 0.5, r = unit(5, "pt"),
        width = unit(5.5, "npc"),
        padding = margin(3, 0, 3, 0),
        margin = margin(1, 1, 1, 1)
      ),
      plot.title = element_text(hjust = 0.5, size = 8, face = 'bold'),
      panel.spacing = unit(1, 'lines')
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    facet_wrap(~ compars, scales = "free", drop = TRUE) +
    geom_text(data = dt_label, aes(label = paste0("R^2 = ", round(value, 3))),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, size = 2, color = "#9A2631")




  plots <- dt_filtered |>
    dplyr::group_split(compars) |>  # 按照 compars 列分割数据
    purrr::map(~ ggplot2::ggplot(.x, ggplot2::aes(x = x, y = y)) +
                 ggplot2::geom_hex(bins = bin) +
                 ggplot2::scale_fill_viridis_c(option = "viridis", guide = "none") +
                 ggplot2::theme_bw() +
                 ggplot2::theme(
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                   panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                   plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
                   axis.title.x = ggplot2::element_text(size = 8, face = 'bold', margin = margin(t = 0.5)),
                   axis.title.y = ggplot2::element_text(size = 8, face = 'bold', margin = margin(r = 0.5)),
                   axis.text = ggplot2::element_text(size = 8, color = "black"),
                   axis.line.x = ggplot2::element_line(color = "black"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8, face = "italic"),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   legend.position = "right",
                   strip.text = ggtext::element_textbox(
                     size = 8, face = 'bold', color = "grey20",
                     hjust = 0.5, halign = 0.5, r = unit(5, "pt"),
                     width = unit(5.5, "npc"),
                     padding = margin(3, 0, 3, 0),
                     margin = margin(1, 1, 1, 1)
                   ),
                   plot.title = element_text(hjust = 0.5, size = 8, face = 'bold'),
                   panel.spacing = unit(1, 'lines')
                 ) +
                 ggplot2::labs(
                   title = NULL,
                   x = .x$xlab[1],  # 使用每组数据的 xlab
                   y = .x$ylab[1]   # 使用每组数据的 ylab
                 ) +
                 ggplot2::geom_text(
                   aes(label = paste0("R^2 = ", round(.x$value[1], 3))),
                   x = -Inf, y = Inf,  # 左上角的位置
                   hjust = -0.1, vjust = 1.1,  # 调整对齐方式以放置在框内的左上角
                   size = 2, color = "#9A2631"
                 )
    )



  if (savefile == "" || is.null(savefile)) {
    savepath <- getwd()
  } else {
    if (dir.exists(savefile)) {
      savepath <- savefile
    } else {
      dir.create(savefile, recursive = TRUE)
      savepath <- savefile
    }
  }



  plot_names <- dt_filtered |>
    dplyr::group_split(compars) |>
    purrr::map_chr(~ unique(.x$compars))

  #  savefile
  # singleplotsize
  #  savepath


  purrr::walk2(plots, plot_names, ~ {
    ggplot2::ggsave(filename = paste0(savepath,"/",.y, ".png"), plot = .x, width = singleplotsize[1], height = singleplotsize[2], dpi = 300)
    ggplot2::ggsave(filename = paste0(savepath,"/",.y, ".pdf"), plot = .x, width = singleplotsize[1], height = singleplotsize[2],device = "pdf")
    export::graph2ppt(x =.x, file =paste0(savepath,"/",.y, ".ppt"),
                      vector.graphic = TRUE,width = singleplotsize[1], height = singleplotsize[2], aspectr = sqrt(2), append = FALSE)
  })

  .save_zcp <- function(Fig,FigName,outputfile,widths,heights){
    Filepaths=paste0(outputfile,"/",FigName,c(".pdf",".png",".ppt"))
    ggplot2::ggsave(Filepaths[1], width =widths, plot=Fig,height = heights,device = "pdf")
    ggplot2::ggsave(Filepaths[2], width =widths,  plot=Fig,height = heights,device = "png")
    export::graph2ppt(x = Fig, file =Filepaths[3],
                      vector.graphic = TRUE, width =widths, height =heights, aspectr = sqrt(2), append = FALSE)
  }

  .save_zcp(Fig =plt,FigName = "facetwrapPlot_cor",outputfile =savepath,widths = facetplotsize[1],heights = facetplotsize[2])

  result=list(facet_plt=plt,singleplot=plots)
  return(result)

}
