hpd <- function (object, level = 0.95, allowSplit = FALSE, ...) 
{
    if (is.na(level) || length(level) != 1 || level <= 0 || level >= 1) 
        stop("level must be between 0 and 1")
    sorted = sort(object$y, decreasing = TRUE)
    heightIdx = min(which(cumsum(sorted) >= sum(object$y) * level))
    height = sorted[heightIdx]
    indices = which(object$y >= height)
    gaps <- which(diff(indices) > 1)
    if (length(gaps) > 0 && !allowSplit) {
        warning("The HPD is discontinuous but allowSplit = FALSE;\n    the result is a valid CrI but not HDI.")
        cumul <- cumsum(object$y)/sum(object$y)
        upp.poss <- low.poss <- which(cumul < 1 - level)
        for (i in low.poss) upp.poss[i] <- min(which(cumul > 
            cumul[i] + level))
        width <- upp.poss - low.poss
        best <- which(width == min(width))
        result <- c(lower = mean(object$x[low.poss[best]]), upper = mean(object$x[upp.poss[best]]))
    } else {
        begs <- indices[c(1, gaps + 1)]
        ends <- indices[c(gaps, length(indices))]
        result <- cbind(begin = object$x[begs], end = object$x[ends])
        if(!allowSplit) {
            result <- as.vector(result)
            names(result) <- c("lower", "upper")
        }
    }
    attr(result, "level") <- level
    attr(result, "height") <- height
    return(result)
}
