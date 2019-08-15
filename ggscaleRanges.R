setGeneric("ggscaleRanges",
           function(gm, transFun, flanking = 1000,  ...) standardGeneric("ggscaleRanges"))

## If it is a GRanges/GRangesList, convert to IRanges and call again  
## NOTE: GenomicRangesORGRangesList AKA GenomicRanges_OR_GenomicRangesList was also called GenomicRanges_OR_GRangesList for a while
## Use whichever class name applies for current env
GRGRL_AKAs = c("GenomicRanges_OR_GRangesList", "GenomicRanges_OR_GRangesList","GenomicRangesORGRangesList") # new,old,deprecated 
GRGRL_ClassName = GRGRL_AKAs[isClass(GRGRL_AKAs)][1]
setMethod("ggscaleRanges", GRGRL_ClassName, 
          function(gm, transFun, flanking = 1000, ...) 
          {
            tryCatch(gm <- unlist(gm), error = function(e) {
              print(paste0('cant unlist gm; attempting without unlisting it. Original error: ', e))
              }
              )
            if (length(unique(seqnames(gm))) > 1){
              stop("Don't know how to transform position coords from different contigs.")
            }
            gm <- IRanges(start = start(gm), end = end(gm))
            ggscaleRanges(gm, transFun, flanking=flanking, ...)
          })

setMethod("ggscaleRanges", "IRanges",
          function(gm, transFun, flanking = 1000, transFunGaps)
          {
            if (missing(transFunGaps)){ 
              transFunGaps <- transFun
            }
            gm.exons <- disjoin(gm)
            gm.introns <- gaps(gm.exons) # If this is done on a GRange, the first gap is extraneous and should be deleted
            df.exons <- as.data.frame(gm.exons)
            df.exons$type <- factor('range', levels = c('range','gap','external'))
            df.introns <- as.data.frame(gm.introns)
            df.introns$type <- factor('gap', levels = c('range','gap','external'))
            df.flat <- rbind(df.exons, df.introns)
            df.flat <- df.flat[order(df.flat$start), ]
            if(flanking > 0){
              df.flat <- rbind(df.flat[1, ], df.flat) # duplicate top row
              df.flat <- rbind(df.flat, df.flat[nrow(df.flat), ]) # duplicate bottom row
              df.flat[1, 1:3] = list(start=df.flat$start[2]-flanking, end=df.flat$start[2] - 1, width = flanking)
              df.flat[nrow(df.flat), 1:3] = list(start = df.flat$end[nrow(df.flat)] + 1, end = df.flat$end[nrow(df.flat)] + flanking, width = flanking)
              df.flat$type[c(1, nrow(df.flat))] = 'external'
            }
            df.flat$twidth[df.flat$type == 'range'] <- transFun(df.flat$width[df.flat$type == 'range'])
            df.flat$twidth[df.flat$type != 'range'] <- transFunGaps(df.flat$width[df.flat$type != 'range'])
            df.flat$twidth <- df.flat$twidth * sum(df.flat$width) / sum(df.flat$twidth) # normalize to total width of original gene model
            df.flat$cum.twidth <- cumsum(df.flat$twidth)
            pre <- c(df.flat$start[1] - 1, df.flat$end)
            post <- c(0, df.flat$cum.twidth) + df.flat$start[1] - 1
            return(data.frame(pre = pre, post = post))
         })

setGeneric("makeSLinkTrack",
           function(mappings, ...) standardGeneric("makeSLinkTrack"))

setMethod("makeSLinkTrack", "data.frame",
          function(mappings, invert = F, bottom_margin = 0.15, inflexions = 'none', ...)
          {
            invert <- as.numeric(invert) * (1 - bottom_margin)
            if (inflexions != 'none'){
              if ('inflexions' %in% colnames(mappings)){
                if(inflexions == 'all'){
                  mappings <- mappings[mappings$inflexions != 0, , drop = F]
                } else if (inflexions == 'positive') {
                  mappings <- mappings[mappings$inflexions > 0, , drop = F]
                } else if (inflexions == 'negative') {
                  mappings <- mappings[mappings$inflexions < 0, , drop = F]
                }
              } else {
                warning('Inflexions requested in s.link plot but no inflexions column in input dataframe.')
              }
            }
            ggplot(mappings) + geom_point(data = mappings[1, , drop = F], aes(x = pre), y = -1) +
              ggplot2::geom_segment(aes(x = pre, xend = post), y = bottom_margin + invert, yend = 1 - invert, ...) +
              theme_null()  
          })
            

setGeneric("approxfun_trans",
           function(x, y, ...) standardGeneric("approxfun_trans"))

setMethod("approxfun_trans", "numeric",
          eval(
            function(x, y, rule = 2:2, ...) {
              filt <- !(is.na(x) | is.na(y))
              if (!all(rank(x[filt]) == rank(y[filt]))){
                warning('Transformation generated for scaling is not monotonic.')
              }
              trans_new("zoned_scale", 
                        transform = approxfun(x = x, y = y, rule = rule, ...), 
                        inverse = approxfun(x = y, y = x, rule = rule, ...)
              )} 
            )
          )

setMethod("approxfun_trans", "data.frame",   #if dataframe, split into 2 vectors and call self
            function(x, rule = 2:2, ...) {
              if(missing(x)) {
                stop('x- and y- data must be provided to the scale transform object to generate the interpolation function, e.g. trans=approxfun_trans(x)')
              }
              if(nrow(x) < 2 | !(all(c('pre', 'post') %in% colnames(x)))  ){
                stop('Dataframe provided to the scale transform object must have > 1 row and contain columns "pre" and "post".')
              }
              approxfun_trans(x = x$pre, y = x$post, rule = rule, ...)
            }
          )


setGeneric("ggscaleDistToRanges",
           function(gm, transFun, flanking = 1000,  ...) standardGeneric("ggscaleDistToRanges"))


setMethod("ggscaleDistToRanges", GRGRL_ClassName, #"GenomicRangesORGRangesList", ## convert GRanges to IRanges, call again
          function(gm, transFun, flanking = 1000, ...) 
          {
            tryCatch(gm <- unlist(gm), error = function(e) {
              print(paste0('cant unlist gm; attempting without unlisting it. Original error: ', e))
              }
              )
            if (length(unique(seqnames(gm))) > 1){
              stop("Don't know how to transform position coords from different contigs.")
            }
            gm <- IRanges(start = start(gm), end = end(gm))
            "ggscaleDistToRanges"(gm, transFun, flanking = flanking, ...)
          })

setMethod("ggscaleDistToRanges", "IRanges",
          function(gm, flanking = 1000, maxprox = 200, minprox = 1)
          {
            distdf <- data.frame(pos = (start(range(gm)) - flanking):(end(range(gm)) + flanking), dist = minprox)
            for(i in 1:length(gm)){
              rn <- gm[i]
              tmp <- data.frame(pos = (start(rn) - maxprox) : (end(rn) + maxprox),
                              dist = c(1:maxprox , rep(maxprox + 1, width(rn)), maxprox:1)
              )
              distdf <- rbind(distdf, tmp)
            }
            distdf <- aggregate(distdf, by = list(distdf$pos), FUN = max)
            prerange <- distdf$pos[nrow(distdf)] - distdf$pos[1]
            postrange <- sum(distdf$dist)
            distdf <- distdf[order(distdf$pos), ]
            distdf$intgrl <- cumsum(distdf$dist)
            distdf$intgrl_trans <- distdf$pos[1] + distdf$intgrl * prerange / postrange
            distdf$deriv1 <- c(0, distdf$dist[2:nrow(distdf)] - distdf$dist[1:(nrow(distdf) - 1)]) # slope of scale value
            distdf$deriv2 <- c(0, distdf$deriv1[2:nrow(distdf)] - distdf$deriv1[1:(nrow(distdf) - 1)]) # 2nd deriv
            distdf <- distdf[ (distdf$deriv1 != 0 ) | (distdf$deriv2 != 0 ) | c(T, rep(F, nrow(distdf) - 2), T), ] # only include rows where scale changes; plus first and last rows
            return(data.frame(pre = distdf$pos, post = distdf$intgrl_trans, inflexions = distdf$deriv2))
            # TODO: probably better stored as an Rle
            })
