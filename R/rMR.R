biglm::biglm
#'@import biglm
#'
background.resp <-
    function(data, DO.var.name, time.var.name = "std.time",
             start.time, end.time, ...){
        
        orig = "1970-01-01 00:00:00 UTC"
        data$y <- eval(parse(text = paste("data$", DO.var.name, sep = "")))
        data$x <- eval(parse(text = paste("data$", time.var.name, sep = "")))
        if(class(start.time) != class(end.time)) {
            stop ("start time and end time must be of same atomic class")
        }
        if(is.character(start.time)==T){
            data <- data[data$x >= as.POSIXct(start.time,
                                              origin = orig) &
                             data$x<= as.POSIXct(end.time,
                                                 origin = orig),]  
        }else if (class(start.time)[1]==("POSIXct") &
                      class(start.time)[2]==("POSIXt")){
            data <- data[data$x >= start.time &
                             data$x <= end.time] 
        }
        
        
        m1 <- biglm(y ~ x, data)
        
        MR1 <- coef(m1)[2]
        if (MR1>=0){warning("slope control 1 negative")}
        plot(y ~ x, data,...)
        abline(coefficients(m1), col="red", lwd=2.5)
        return(summary(m1))
    }

#'@export
biglm::biglm


Barom.Press <-
    function(elevation.m, units = "atm"){
        if(units == "atm") {fact <-1}
        else if(units == "kpa") {
            fact <- 101.325
        }else if(units == "mmHg") {
            fact <- 760
        }else{
            stop("invalid pressure units, must be
             'atm', 'kpa', or 'mmHg'")
        }
        P<-exp(-9.80665 * 0.0289644 * elevation.m / 
                   (8.31447 * 288.15))*fact
        return(P)
    }

#'@export
#'@import biglm
#'
DO.saturation <-
    function(DO.mgl, temp.C, elevation.m = NULL,
             bar.press = NULL, bar.units = "atm"){
        DO.sat<- DO.mgl / Eq.Ox.conc(temp.C, elevation.m,
                                     bar.press, bar.units)
        return(DO.sat)
    }
#'@export
#'@import biglm 
#'


DO.unit.convert <-
    function(x, DO.units.in, DO.units.out, 
             bar.units.in, bar.press, temp.C,
             bar.units.out){
        if(bar.units.in == bar.units.out){
            bar.press <- bar.press
        }else if(bar.units.in == "atm"){
            if(bar.units.out == "mmHg"){
                bar.press <- bar.press * 760
            }else if(bar.units.out == "kpa"){
                bar.press <- bar.press * 101.32501
            }else{
                stop("invalid 'bar.units.out' -- must be 'atm', 'mmHg', 'kpa'")
            }
        }else if(bar.units.in == "mmHg"){
            if(bar.units.out == "atm"){
                bar.press <- bar.press / 760
            }else if(bar.units.out == "kpa"){
                bar.press <- bar.press * 101.32501 / 760
            }else{
                stop("invalid 'bar.units.out' -- must be 'atm', 'mmHg', 'kpa'")
            }
        }else if(bar.units.in == "kpa"){
            if(bar.units.out == "atm"){
                bar.press <- bar.press / 101.32501
            }
            if(bar.units.out == "mmHg"){
                bar.press <- bar.press * 760 / 101.32501
            }
        }
        
        
        if (DO.units.in == "pct"){
            DO.pct <- x / 100
        }else{
            eq.o2.in <- Eq.Ox.conc(temp.C, bar.units = bar.units.out, 
                                   bar.press =bar.press,
                                   out.DO.meas = DO.units.in)
            DO.pct <- x / eq.o2.in
        }
        
        if (DO.units.out == "pct"){
            DO.conc <- DO.pct * 100
        }else{
            eq.o2.out <- Eq.Ox.conc(temp.C, bar.units = bar.units.out,
                                    bar.press = bar.press,
                                    out.DO.meas = DO.units.out)
            DO.conc <- DO.pct * eq.o2.out
        }
        
        return(DO.conc)
    }
#'@export
#'@import biglm
#'
Eq.Ox.conc <-
    function(temp.C, elevation.m = NULL,
             bar.press = NULL, bar.units = "atm",
             out.DO.meas = "mg/L"){
        
        if( out.DO.meas == "PP"){
            
            if(is.null(bar.press) == FALSE && is.null(elevation.m) == TRUE){
                bar.press <- bar.press
            } else if(is.null(bar.press) == TRUE &&
                          is.null(elevation.m) == FALSE){
                bar.press <- Barom.Press (elevation.m, units = bar.units)
            }else{
                stop("EITHER 'elevation.m' or 'barom.press' must be assigned
                     a value. The other argument must be NULL.")
            }
            
            Cp <-bar.press*0.20946 
            
            }else if (out.DO.meas == "mg/L"){
                std.Eq.Ox.conc <- exp(7.7117 - 1.31403*(log(temp.C + 45.93)))
                if(is.null(bar.press) == FALSE && is.null(elevation.m) == TRUE){
                    if(bar.units == "atm"){
                        bar.press <- bar.press
                    }else if(bar.units == "kpa"){
                        bar.press <- bar.press / 101.325
                    }else if(bar.units == "mmHg"){
                        bar.press <- bar.press / 760
                    }else{
                        stop("invalid pressure units, must be
                         'atm', 'kpa', or 'mmHg'")
                    }
                } else if(is.null(bar.press) == TRUE &&
                              is.null(elevation.m) == FALSE){
                    bar.press <- Barom.Press (elevation.m, units = bar.units)
                }else{
                    stop("EITHER 'elevation.m' or 'barom.press' must be assigned
                         a value. The other argument must be NULL.")
                }
                
                temp.K <- 273.15 + temp.C
                P.H2Ovap <- exp(11.8571 - (3840.7 / temp.K) - (216961/(temp.K^2)))
                theta <- 0.000975 - (1.426e-5 * temp.C) +(6.436e-8 * temp.C^2) 
                
                Cp <- std.Eq.Ox.conc*bar.press * 
                    ((1-(P.H2Ovap/bar.press))*1-(theta*bar.press))/
                    ((1-P.H2Ovap)*(1-theta))
            }else{
                stop("must specify 'out.DO.meas' as 'mg/L' or 'PP'")
            }
        
        return(Cp)
    }

#'@export
#'@import biglm
#'

get.pcrit <-
    function(data, DO = NULL, MR = NULL, Pcrit.below,
             idx.interval = NULL, index.var = "std.time",
             start.idx = NULL, stop.idx = NULL, ...){
        
        data$DO <- eval(parse(text = paste("data$", DO, sep = "")))
        
        if(any(is.na(data$DO)==T)){
            warning("DO variable contains missing values")
        } 
        
        if(is.null(MR) == TRUE){
            
            if(is.null(idx.interval)){
                stop("if 'MR' = NULL, 'interval' must be specified")
            }
            
            
            data$idx <- eval(parse(text =
                                       paste("data$", index.var, sep = "")))
            
            
            data <- data[data$idx >= start.idx
                         & data$idx <= stop.idx,]
            
            
            data$idx <- as.numeric(data$idx) - min(as.numeric(data$idx))
            
            calc.MRs <- data.frame()
            
            i<-0
            
            while(i + idx.interval < length(data$idx)){
                the.interval <- data[data$idx >= i &
                                         data$idx < i + idx.interval,]
                m <- lm(DO ~ idx, data = the.interval)
                MR <- m$coefficients[2]*(-1)
                DO <- mean(the.interval$DO, na.rm = TRUE)
                time <- round(the.interval$idx[1])
                df.row <- t(c(DO, MR))
                calc.MRs <- rbind(calc.MRs, df.row)
                i <- i + idx.interval
                
            }
            
            
            data <- calc.MRs
            names(data) <- c("DO", "MR")
            
        }else if(is.null(MR) == FALSE){
            data$MR <- eval(parse(text = paste("data$", "MR", sep = "")))
        }
        
        
        
        
        minimum.DO <- 
            as.numeric(min(data[data$DO <= Pcrit.below, 1], na.rm=TRUE))
        
        zone.of.interest <- data[data$DO <= Pcrit.below,]
        zone.of.interest <- zone.of.interest[order(zone.of.interest$DO, 
                                                   decreasing=TRUE),]
        
        RSS.table<-data.frame()
        for(w in 2:length(zone.of.interest[,1])-1){
            if(zone.of.interest$DO[w] != minimum.DO){
                break.point <- zone.of.interest$DO[w]
                RSS <- tot.rss(data = data,
                               break.pt = break.point,
                               xvar = "DO", yvar = "MR")
                RSS.row <- cbind(RSS, zone.of.interest$DO[w])
                RSS.table <- rbind(RSS.table,RSS.row)    
            }
        }
        
        DO.idx <- RSS.table[RSS.table$RSS == min(RSS.table$RSS), 2]
        
        
        
        ##creating linear models##
        dat1 <- data[data$DO > DO.idx,]
        mod.1 <- lm(dat1$MR~dat1$DO)
        
        dat2 <- data[data$DO <= DO.idx,]
        mod.2 <- lm(dat2$MR ~ dat2$DO)
        
        sm1 <- summary(mod.1)
        sm2 <- summary(mod.2)
        
        adjr2.pre <- sm1$adj.r.squared
        adjr2.post <- sm2$adj.r.squared
        
        ## Plot ##
        
        plot(MR~DO, data, type= "n", ...)
        points(x = c(dat1$DO, dat2$DO), y = c(dat1$MR, dat2$MR), 
               cex = .7, ...)
        intersect<-(mod.2$coefficients[1] - mod.1$coefficients[1]) /
            (mod.1$coefficients[2] - mod.2$coefficients[2])
        names(intersect)<-NULL
        if(is.na(mod.2$coefficients[2])==T){
            abline(mod.1$coefficients[1], 0)
        }else{abline(coef = mod.1$coefficients, col="red", lwd = 1.5, ...)}
        abline(coef = mod.2$coefficients, col = "red", lwd = 1.5, ...)
        points(x = intersect, y = mod.1$coefficients[1] +
                   mod.1$coefficients[2]*intersect, pch=21,
               col = "blue", cex=1.5)
        dat.pre<-data[data$DO>=(2*intersect),]
        
        P.crit<-as.data.frame(cbind(intersect, adjr2.pre, adjr2.post))
        names(P.crit)<-c("Pcrit", "Adj.r2.above", "Adj.r2.below")
        above.Pc <- mod.1
        names(above.Pc) <- paste("above",names(above.Pc),
                                 sep=".")
        below.Pc <- mod.2
        names(below.Pc) <- paste("below",names(below.Pc),
                                 sep=".")
        Pc<-as.list(c(P.crit, above.Pc,below.Pc))
        return(Pc)
        names(Pc)
    }
#'@export
#'@import biglm
#'
#'
get.witrox.data <-
    function(data.name, lines.skip, delimit="tab", choose.names = F,
             chosen.names = NULL,
             format){
        if(delimit == "tab"){
            separate = "\t"
        }else if(delimit == "space"){
            separate = ""
        }else if(delimit == "comma"){
            separate = ","
        }
        
        if(choose.names==F){
            d<-read.table(file.choose(), sep=separate, skip=lines.skip,
                          header =T,check.names=F)
            invalid.names<-colnames(d)
            valid.names<-make.names(colnames(d))
            var.names<-NULL
            for(i in 1:length(invalid.names)){
                if(invalid.names[i] == valid.names[i]){
                    var.names[i] <- valid.names[i]
                }else{
                    split.name.period <- as.vector(
                        strsplit(valid.names[i], fixed = T,
                                 split = ".")[[1]])
                    split.name.period <- paste(split.name.period, sep="",
                                               collapse="")
                    split.name <- as.vector(
                        strsplit(split.name.period, fixed = F,
                                 split = c("Â"))[[1]])
                    split.name <- paste(split.name, sep=".", collapse=".")
                    var.names[i] <-  as.vector(
                        strsplit(split.name, fixed = F,
                                 split = c("â"))[[1]])
                }
            }
        }else if(choose.names==T){
            d<-read.table(data.name, sep = separate, skip = lines.skip,
                          header = F)
            var.names<-chosen.names
        }
        
        colnames(d) <- var.names
        d[,1]<-as.character(d[,1]) 
        d$std.time <- as.POSIXct(strptime(d[,1], format = format),
                                 origin = "1970-01-01 00:00:00 UTC")
        return(d)
        
    }
#'@export
#'@import biglm

MR.loops <-
    function(data, DO.var.name, time.var.name = "std.time",
             in.DO.meas = "mg/L", out.DO.meas = "mg/L",
             start.idx, stop.idx, system.vol = 1,
             background.consumption = NULL,
             background.start = NULL, background.end = NULL,
             temp.C, elevation.m = NULL,
             bar.press = NULL, bar.units = "atm", ...){
        ## format the time vectors into POSIX ##
        orig = "1970-01-01 00:00:00 UTC"
        start.idx <- as.POSIXct((start.idx), origin = orig)
        stop.idx <- as.POSIXct((stop.idx), origin = orig)
        
        if (length(start.idx) != length(stop.idx)){
            stop ("number of start times not equal
                  to number of stop times")
        }
        
        ## set response variable ##
        data$y <- eval(parse(text = paste("data$", DO.var.name, sep = "")))
        ## set time variable ##
        data$x <- eval(parse(text = paste("data$", time.var.name, sep = "")))
        
        ## set background DO consumption rate ##
        if(is.null(background.consumption) == TRUE){
            bgd.info <- background.resp(DO.var.name = "y", data, 
                                        time.var.name = "x",
                                        background.start, background.end)
            bgd.slope <- bgd.info$mat[2,1]
        }else if(is.null(background.consumption) == FALSE) {
            bgd.slope = background.consumption
        }
        
        if(is.null(bar.press) == FALSE && is.null(elevation.m) == FALSE){
            stop("Either 'bar.press' or 'elevation.m' should be NULL")
        }
        
        ## barometric pressure ##
        if(is.null(bar.press) == FALSE){
            if(is.character(bar.press) == TRUE){
                bar.press <- eval(parse(
                    text = paste("data$",
                           bar.press, sep = "")))
            }else if(is.numeric(bar.press) == TRUE){
                bar.press <- bar.press
            }else{
               stop("'bar.press' must be 'NULL', numeric, or
                      the col.name for barometric pressure")
            }
            }
        
        
        
        ## Temperature ##
        if (is.character(temp.C) == TRUE){
            temp.C <- eval(parse(
                text = paste("data$", temp.C, sep = "")))
        }else if(is.numeric(temp.C) == TRUE){
            temp.C <- temp.C
        }else{
            stop("invalid temp.C argument")
        }
        
        # DO sat conversions #
        if (in.DO.meas == "pct"){
            data$y <- (data$y /100) * 
                Eq.Ox.conc(temp.C, elevation.m,
                           bar.press, bar.units,
                           out.DO.meas)
        }else if (in.DO.meas == "PP"){
            fraction <- data$y / 
                Eq.Ox.conc(temp.C, elevation.m,
                           bar.press, bar.units,
                           out.DO.meas = "PP")
            data$y <- fraction *
                Eq.Ox.conc(temp.C, elevation.m,
                           bar.press, bar.units,
                           out.DO.meas)
        }else if(in.DO.meas == "mg/L"){
            data$y <- data$y
        }else{
            stop("invalid 'in.DO.meas' argument:
                 must be 'pct', 'PP', or 'mg/L' ")
        }
        
        data$y <- system.vol * data$y
        # Now data$y in units of mg #
        
        data$adj.y <- data$y - (as.numeric(data$x) -
                      as.numeric(start.idx[1]))*bgd.slope
        
        plot(adj.y ~ x,
             data = data[(data$x >= start.idx[1] - 600) &
                    data$x <= (tail(stop.idx,1) + 600),],
                    type="n",...)
        
        name.num<-as.character(c(1:length(start.idx)))
        ms<-list()
        MR.summary<-data.frame()
        for(i in 1:length(start.idx)){
            dat <- data[data$x >= start.idx[i]
                        & data$x <= stop.idx[i],]
            
            dat$y <- dat$y - (as.numeric(dat$x-
                                             as.numeric(start.idx[i]))*bgd.slope )
            
            mk <- biglm(y ~ x, dat)
            ms[[i]] <- mk
            
            points(dat$x, dat$y)
            names(ms[[i]])<-paste(names(ms[[i]]), name.num[i], sep=".")
            abline(coef(ms[[i]]),
                   col="red",  lwd = 2)
            
            MR <- coef(mk)[2]*-1
            sds <- summary(mk)$mat[2,4]*sqrt(length(dat[,1]))
            rsquare <- summary(mk)$rsq
            mrrow <- t(c(MR, sds, rsquare))
            MR.summary <- rbind(MR.summary,mrrow)
        }
        names(MR.summary) <- c("MR", "sd.slope", "r.square")
        ofthejedi <- list(MR.summary, ms)
        names(ofthejedi) <- c("MR.summary", "lm.details")
        return(ofthejedi)
        }

#'@export
#'@import biglm

plot.raw <-
    function(data, DO.var.name, time.var.name = "std.time",
             start.time = data$x[1],
             end.time = data$x[length(data$x)], ...){
        orig = "1970-01-01 00:00:00 UTC"
        
        data$x <- eval(parse(text = paste("data$", time.var.name, sep = "")))
        data <- data[data$x >= start.time & data$x <= end.time,]
        data$y <- eval(parse(text = paste("data$", DO.var.name, sep = "")))
        
        plot(y~std.time, data, ...)
    }

#'@export
#'@import biglm

sumsq <-
    function(x){
        mx<-mean(x)
        sq.dev<-((x-mx)^2)
        sumsquares<-sum(sq.dev)
        return(sumsquares)
    }
#'@export
#'@import biglm
#'

tot.rss <-
    function(data, break.pt, xvar, yvar){
        data$x <- eval(parse(text = paste("data$", xvar, sep="")))
        data$y <- eval(parse(text = paste("data$", yvar, sep="")))
        d1 <- data[data$x >= break.pt,]
        m1 <- lm(d1$y~d1$x)
        
        d2 <- data[data$x < break.pt,]
        m2 <- lm(d2$y~d2$x)
        
        trss <- (sumsq(m1$residuals)+sumsq(m2$residuals))
        return(trss)
    }
#'@export
#'@import biglm
#'

