#' barplot by ggplot2
#' @description barplot with number in top.
#' @param height data for plot
#' @param fill,xlab,ylab parameters pass to ggplot.
#' @param draw plot or not
#' @param postfix Postfix of text labled in top of bar.
#' @return ggplot object.
#' @importFrom ggplot2 ggplot geom_bar theme_classic geom_text aes xlab ylab
#' @importFrom ggfittext geom_bar_text
#' @examples
#' ribosomeProfilingQC:::ggBar(sample.int(100, 3))
ggBar <- function(height, fill="gray80", draw=TRUE, xlab, ylab, postfix){
  if(length(names(height))!=length(height)){
    n <- seq_along(height)
  }else{
    n <- names(height)
  }
  height <- as.numeric(height)
  label <- formatC(height)
  x <- factor(n, levels = unique(n))
  y <- height
  df <- data.frame(x=x, y=y,
                   label=label, fill=fill)
  if(!missing(postfix)){
    if(!is.logical(postfix)){
      df$label <- paste(as.character(df$label), postfix)
      postfix <- TRUE
    }
  }else{
    postfix <- FALSE
  }
  plot <- ggplot(data=df, aes(x=x, y=y)) +
    geom_bar(stat = "identity", fill=fill) +
    theme_classic()
  if(!missing(postfix)) suppressWarnings(
    plot <- plot +
      geom_bar_text(aes(label=label), grow = FALSE, contrast = TRUE))
  if(!missing(xlab)) plot <- plot+xlab(xlab)
  if(!missing(ylab)) plot <- plot+ylab(ylab)
  if(draw) {
    print(plot)
  }else{
    plot
  }
}
