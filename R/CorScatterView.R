#' Plot scatter plot and label correlations
#'
#' @docType methods
#' @name CorScatterView
#' @rdname CorScatterView
#'
#' @param dat a data.frame object
#' @param x X-axis for plotting
#' @param y Y-axis for plotting
#' @param facet Facet groups
#' @param nrow The number of rows for facet_wrap
#' @param cor.method pearson or spearman
#' 
#' @import ggplot
#' @export
CorScatterView <- function(dat, x, y, facet = NULL, nrow = 2, cor.method = "pearson"){
    gg <- dat[, c(x, y)]
    colnames(gg) <- c("X", "Y")
    if(!is.null(facet)){ gg$Group <- dat[, facet] }else{ gg$Group <- "" }
    gg <- cbind(gg, dat[, !colnames(dat)%in%colnames(gg)])
    lab1 <- lapply(unique(gg$Group), function(x){
        idx <- gg$Group==x
        if(sum(idx)<3) return(NULL)
        cor1 <- cor.test(gg$X[idx], gg$Y[idx], method = cor.method)
        lab1 <- paste0("r = ", round(cor1$estimate, 3), 
                       "\np = ", format(cor1$p.value, digits = 2, scientific = TRUE))
        if(cor.method=="spearman"){
            lab1 <- paste0("rho = ", round(cor1$estimate, 3), 
                           "\np = ", format(cor1$p.value, digits = 2, scientific = TRUE))
        }
        data.frame(X = min(gg$X[idx]), Y = max(gg$Y[idx]), label = lab1, Group = x)
    })
    lab1 <- do.call(rbind, lab1)
    p1 <- ggplot(gg, aes(X, Y)) + geom_point() + geom_smooth(method = "lm") + theme_classic(base_size = 14) + 
            labs(x = x, y = y) + geom_text(aes(X, Y, label = label), hjust = 0, vjust = 1, data = lab1)
    if(!is.null(facet)){ 
        p1 <- p1 + facet_wrap(~Group, scales = "free", nrow = nrow) + theme(strip.background = element_blank())
    }
    p1
}
