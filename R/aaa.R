#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter select case_when pull
#' @importClassesFrom tibble tbl_df
#' @import Matrix
#' @import ggplot2
NULL

#' CytoSignal Color Palette
#' @description
#' A vector of colors used in visualizing discrete values. Originally taken from
#' \code{ggsci::pal_igv(palette = 'default')(51)}
#' @export
csColors <- c(
    '#5050FF', '#CE3D32', '#749B58', '#F0E685', '#466983', '#BA6338',
    '#5DB1DD', '#802268', '#6BD76B', '#D595A7', '#924822', '#837B8D',
    '#C75127', '#D58F5C', '#7A65A5', '#E4AF69', '#3B1B53', '#CDDEB7',
    '#612A79', '#AE1F63', '#E7C76F', '#5A655E', '#CC9900', '#99CC00',
    '#A9A9A9', '#CC9900', '#99CC00', '#33CC00', '#00CC33', '#00CC99',
    '#0099CC', '#0A47FF', '#4775FF', '#FFC20A', '#FFD147', '#990033',
    '#991A00', '#996600', '#809900', '#339900', '#00991A', '#009966',
    '#008099', '#003399', '#1A0099', '#660099', '#990080', '#D60047',
    '#FF1463', '#00D68F', '#14FFB1'
)
