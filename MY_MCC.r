#Purpose: Estimating the Burden of Recurrent Events in the Presence of Competing 
#         Risks:The Method of Mean Cumulative Count
#Reference: Dong et al. 2015 (DOI: 10.1093/aje/kwu289)
#Author: Grace Zhou
#Date: 2024-11-10
#Code reference: MCCfunctions_09132016_GH.r
#Note: This code has been verified
################################################################################

#Load library
library(cmprsk) #cumulative incidence for right censoring only + competing risk
library(dplyr)  #data management
library(mstate) #cumulative incidence for left truncation and/or right censoring + competing risk

#revised crprep function
{
  # Function to create weighted data set for competing risks analyses
  # Author: Ronald Geskus
  # Reference: Geskus RB, 2011. Cause-Specific Cumulative Incidence Estimation
  # and the Fine and Gray Model Under Both Left Truncation and Right Censoring.
  # Biometrics, 67, 39--49.
  # Link: https://www.rdocumentation.org/packages/mstate/versions/0.3.3/topics/crprep.default
  # Revised by Grace Zhou (grace.zhou@stjude.org) in November 2024 to ensure accurate weights

  crprep2 <- function(Tstop, status, data, trans=1, cens=0, Tstart=0, id, strata, 
             keep, shorten=TRUE, rm.na=TRUE, origin=0,
             prec.factor=1000, ...) {
    
    # Coerce data to data.frame if given (from tibble/data.table)
    if (!missing(data)) data <- as.data.frame(data)
    
    ## Extract Tstop data if given by column name
    if (!(is.numeric(Tstop))) {
      if (!is.character(Tstop) | length(Tstop)!=1)
        stop("argument \"Tstop\" should be a numeric vector or a character string")
      if (missing(data))
        stop("missing \"data\" argument not allowed when \"Tstop\" argument is a character string")
      tcol <- match(Tstop, names(data))
      if (is.na(tcol))
        stop("\"Tstop\" not found in data")
      Tstop <- data[, tcol]
    } else {
      if (!is.vector(Tstop))
        stop("argument should be a numeric vector or a character string")
    }
    nn <- length(Tstop)
    
    ## Extract Tstart data if given by column name
    if (is.numeric(Tstart)&length(Tstart) == 1) {
      Tstart <- rep(Tstart, nn)
    } else {
      if (!(is.numeric(Tstart))) {
        if (!is.character(Tstart) | length(Tstart)!=1)
          stop("argument \"Tstart\" should be a numeric vector or a character string")
        if (missing(data))
          stop("missing \"data\" argument not allowed when \"Tstart\" argument is a character string")
        tcol <- match(Tstart, names(data))
        if (is.na(tcol))
          stop("\"Tstart\" not found in data")
        Tstart <- data[, tcol]
      } else {
        if (!is.vector(Tstart))
          stop("argument should be a numeric vector or a character string")
      }
    }
    if (length(Tstart) != nn)
      stop("Tstop and Tstart have different lengths")
    ## Check whether Tstart is needed
    calc.trunc <- any(Tstart[!is.na(Tstart)] != 0)
    
    ## Select rows without missing time value
    sel <- !is.na(Tstart) & !is.na(Tstop)
    if (any(Tstart[sel] >= Tstop[sel]))
      stop("Tstop must be greater than Tstart")
    
    ## Extract status data if given by column name
    if (length(status) == 1) {
      if (!is.character(status))
        stop("argument \"status\" should be a vector or a character string")
      if (missing(data))
        stop("missing \"data\" argument not allowed when \"status\" argument is a character string")
      tcol <- match(status, names(data))
      if (is.na(tcol))
        stop("\"status\" not found in data")
      status <- data[ ,tcol]
    }
    if (length(status) != nn)
      stop("Tstop and status have different lengths")
    
    ## Extract strata data; value 1 if not specified
    if (missing(strata)) {
      strata.val <- rep(1,nn)
    } else {
      if (is.matrix(strata) | is.data.frame(strata))
        stop("only one variable is allowed in \"strata\"")
      if (!(is.vector(as.numeric(factor(strata))) & length(strata) > 1)) {
        if (!is.character(strata))
          stop("argument \"strata\" should be a character string")
        if (missing(data))
          stop("missing \"data\" argument not allowed when \"strata\" argument is a character string")
        tcol <- match(strata, names(data))
        if (is.na(tcol))
          stop("\"strata\" not found in data")
        strata.name <- strata
        strata.val <- data[ ,tcol]
      } else {
        if (length(strata) != nn)
          stop("Tstop and strata have different lengths")
        strata.name <- names(strata)
        strata.val <- strata
      }
    }
    strata.num <- as.numeric(factor(strata.val))
    
    ## Extract id data; values 1:nn if not specified
    if (missing(id)) {
      id.name <- "id"
      id  <- num.id <- 1:nn
    } else {
      if (is.matrix(id) | is.data.frame(id))
        stop("only one variable is allowed in \"id\"")
      if (!(is.vector(id) & length(id) > 1)) { # by name
        if (!is.character(id))
          stop("argument \"id\" should be a character string")
        if (missing(data))
          stop("missing \"data\" argument not allowed when \"id\" argument is a character string")
        tcol <- match(id, names(data))
        if (is.na(tcol))
          stop("\"id\" not found in data")
        id.name <- id
        num.id <- 1:nn
        id <- data[, tcol]
      } else {                                 # by value
        if (length(id) != nn)
          stop("Tstop and id have different lengths")
        id.name <- names(id)
        num.id <- 1:nn
      }
    }
    
    ## Eliminate records with missings in status if rm.na=TRUE
    if(rm.na) sel <- sel & !is.na(status)
    Tstart <- Tstart[sel]
    Tstop <- Tstop[sel]
    status <- status[sel]
    strata.val <- strata.val[sel]
    strata.num <- strata.num[sel]
    id <- id[sel]
    num.id <- num.id[sel]
    n <- length(Tstop)
    
    ## Extract covariate data
    if(!missing(keep)) {
      if (!(is.matrix(keep) | is.data.frame(keep))) {
        if (is.character(keep)) {  # if given by column name
          if (missing(data))
            stop("missing \"data\" argument not allowed when \"keep\" argument is a character vector")
          nkeep <- length(keep)
          kcols <- match(keep, names(data))
          if (any(is.na(kcols)))
            stop("at least one element of \"keep\" not found in data")
          keep.name <- keep
          keep <- data[, kcols]
        } else {                     # if one column, given by value
          nkeep <- 1
          ##        keep.name <- names(as.data.frame(keep))
          keep.name <- names(keep)
          ##        if(is.null(keep.name)) keep.name <- "V1"
          if (length(keep) != nn)
            stop("Tstop and keep have different lengths")
        }
      } else {                         # if a matrix/data.frame
        nkeep <- ncol(keep)
        if(is.data.frame(keep)){
          keep.name <- names(keep)
        } else {
          keep.name <- colnames(keep)
          if(is.null(keep.name)) keep.name <- paste("V",1:nkeep,sep="")
        }
        if (nrow(keep) != nn)
          stop("length Tstop and number of rows in keep are differents")
        if (nkeep == 1)
          keep <- keep[, 1]
      }
    }
    
    Tstart <- Tstart - origin
    Tstop <- Tstop - origin
    
    
    ## Start calculations
    prec <- .Machine$double.eps*prec.factor
    
    ## Calculate product-limit time-to-censoring distribution, "event" not included in case of ties
    surv.cens <- survival::survfit(Surv(Tstart,Tstop+ifelse(status==cens,prec,0),status==cens)~strata.num)
    
    ## Calculate time to entry (left truncation) distribution at t-, use 2*prec in order to exclude censorings at same time
    if(calc.trunc) surv.trunc <- survival::survfit(Surv(-Tstop,-(Tstart+2*prec),rep(1,n))~strata.num)
    ## trunc.dist <- summary(surv.trunc)
    ## trunc.dist$time <- rev(-trunc.dist$time)-prec
    ## trunc.dist$surv <- c(rev(trunc.dist$surv)[-1],1)
    
    ## Create weighted data set for each event type as specified in trans
    data.out <- vector("list",length(trans))
    i.list <- 1
    strat <- sort(unique(strata.num),na.last=TRUE)
    len.strat <- length(strat)
    ## Start weight calculation per event type
    for(failcode in trans) {
      if(len.strat==1){ # no strata
        data.weight <- mstate:::create.wData.omega(Tstart, Tstop, status, num.id, 1, failcode, cens)
        tmp.time <- data.weight$Tstop
        #data.weight$weight.cens[order(tmp.time)] <- summary(surv.cens, times=tmp.time-prec)$surv #Grace modify
        data.weight$weight.cens <- summary(surv.cens, times=tmp.time-prec)$surv #Grace modify
        #if(calc.trunc) data.weight$weight.trunc[order(-tmp.time)] <- summary(surv.trunc, times=-tmp.time)$surv #Grace modify
        if(calc.trunc) data.weight$weight.trunc <- summary(surv.trunc, times=-tmp.time)$surv #Grace modify
      } else {
        data.weight <- vector("list",len.strat)
        if(is.na(strat[len.strat])) {
          tmp.sel <- is.na(strata.num)
          data.weight[[len.strat]] <- data.frame(id=num.id[tmp.sel], Tstart=Tstart[tmp.sel], Tstop=Tstop[tmp.sel], status=status[tmp.sel], strata=NA,  weight.cens=NA)
          if(calc.trunc) data.weight[[len.strat]]$weight.trunc <- NA
        }
        for(tmp.strat in 1:(len.strat-is.na(strat[len.strat]))){
          tmp.sel <- !is.na(strata.num) & strata.num==tmp.strat
          data.weight[[tmp.strat]] <- create.wData.omega(Tstart[tmp.sel], Tstop[tmp.sel], status[tmp.sel], num.id[tmp.sel], tmp.strat, failcode, cens)
          tmp.time <- data.weight[[tmp.strat]]$Tstop
          #data.weight[[tmp.strat]]$weight.cens[order(tmp.time)] <- summary(surv.cens[tmp.strat], times=tmp.time-prec)$surv #Grace modify
          data.weight[[tmp.strat]]$weight.cens <- summary(surv.cens[tmp.strat], times=tmp.time-prec)$surv #Grace modify
          #if(calc.trunc) data.weight[[tmp.strat]]$weight.trunc[order(-tmp.time)] <- summary(surv.trunc[tmp.strat], times=-tmp.time)$surv #Grace modify
          if(calc.trunc) data.weight[[tmp.strat]]$weight.trunc <- summary(surv.trunc[tmp.strat], times=-tmp.time)$surv #Grace modify
        }
        data.weight <- do.call("rbind", data.weight)
      }
      ## Calculate omega-censoring weights
      data.weight <- data.weight[order(data.weight$id,data.weight$Tstop), ]
      data.weight$weight.cens <- unlist(tapply(data.weight$weight.cens, data.weight$id, FUN=function(x) if(length(x)==1&!is.na(x[1])) 1 else x/x[1]))
      
      ## Calculate omega-truncation weights
      if(calc.trunc) {
        data.weight$weight.trunc <- unlist(tapply(data.weight$weight.trunc, data.weight$id, FUN=function(x) if(length(x)==1&!is.na(x[1])) 1 else x/x[1]))
      }
      
      tbl <- table(data.weight$id)
      
      ## Add covariates
      if(!missing(keep)) {
        ## Extract covariate name from function
        if (is.null(keep.name)) {
          m <- match.call(expand.dots = FALSE)
          m <- m[match("keep", names(m))]
          if(!is.null(m)) {
            keep.name <- as.character(m[1])
            keep.name.split <- strsplit(keep.name, '')[[1]]
            tag <- which(keep.name.split == '$')
            if(length(tag) != 0) {
              keep.name <- substring(keep.name, tag[length(tag)]+1)
            } else {
              tag <- which(keep.name.split == '"')
              if(length(tag) != 0) {
                keep.name <- substring(keep.name, tag[1]+1, tag[2]-1)
              }
            }
          }
        }
        
        ## Add covariates to the resultset
        if (nkeep > 0) {
          if (nkeep == 1) {
            keep <- keep[sel]
            ddcovs <- rep(keep, tbl)
            ddcovs <- as.data.frame(ddcovs)
            names(ddcovs) <- as.character(keep.name)
          } else {
            keep <- keep[sel, ]
            ddcovs <- lapply(1:nkeep, function(i) rep(keep[, i], tbl))
            ddcovs <- as.data.frame(ddcovs)
            names(ddcovs) <- keep.name
          }
          data.weight <- cbind(data.weight, ddcovs)
        }
      }
      
      ## Shorten data set by combining rows with event types without censoring or truncation time in between
      if (shorten) {
        if(calc.trunc) {
          keep.rows <- with(data.weight, c(diff(id)!=0 | diff(weight.cens)!=0 | diff(weight.trunc)!=0, TRUE))
        } else {
          keep.rows <- with(data.weight, c(diff(id)!=0 | diff(weight.cens)!=0, TRUE))
        }
        ## First record always included as separate row in order to allow for CSH analysis
        keep.rows[unlist(mapply(seq,1,tbl))==1] <- TRUE
        keep.start <- data.weight$Tstart[unlist(tapply(keep.rows, data.weight$id, FUN=function(x) if(length(x)==1) x else c(TRUE,x[-length(x)])))]
        data.weight <- data.weight[keep.rows,]
        data.weight$Tstart <- keep.start
      }
      
      ## Recalculate tbl after shorten
      tbl <- table(data.weight$id)
      ## Add count
      data.weight$count <- unlist(mapply(seq,1,tbl))
      data.weight$failcode <- failcode
      ## Return to original id
      data.weight$id <- rep(id,tbl)
      
      data.out[[i.list]] <- data.weight
      i.list <- i.list+1
    }
    
    out <- do.call("rbind", data.out)
    
    if(!missing(strata)) {
      ## Extract strata name if given by value
      if (is.null(strata.name)) {
        m <- match.call(expand.dots = FALSE)
        m <- m[match("strata", names(m))]
        if(!is.null(m)) {
          strata.name <- as.character(m[1])
          strata.name.split <- strsplit(strata.name, '')[[1]]
          tag <- which(strata.name.split == '$')
          if(length(tag) != 0) {
            strata.name <- substring(strata.name, tag[length(tag)]+1)
          } else {
            tag <- which(strata.name.split == '"')
            if(length(tag) != 0) {
              strata.name <- substring(strata.name, tag[1]+1, tag[2]-1)
            }
          }
        }
      }
      ## Use original stratum values
      if(is.factor(strata.val)) {
        out$strata <- factor(out$strata, labels=levels(strata.val))
      } else {
        out$strata <- levels(factor(strata.val))[out$strata]
        if(is.numeric(strata.val)) out$strata <- as.numeric(out$strata)
      }
      ## Use original name of column
      tmp.sel <- match("strata", names(out))
      names(out)[tmp.sel] <- strata.name
    } else {
      out$strata <- NULL
    }
    
    if (is.null(id.name)) {
      m <- match.call(expand.dots = FALSE)
      m <- m[match("id", names(m))]
      if(!is.null(m)) {
        id.name <- as.character(m[1])
        id.name.split <- strsplit(id.name, '')[[1]]
        tag <- which(id.name.split == '$')
        if(length(tag) != 0) {
          id.name <- substring(id.name, tag[length(tag)]+1)
        } else {
          tag <- which(id.name.split == '"')
          if(length(tag) != 0) {
            id.name <- substring(id.name, tag[1]+1, tag[2]-1)
          }
        }
      }
    }
    
    row.names(out) <- as.character(1:nrow(out))
    names(out)[1] <- id.name
    class(out) <- c("crprep","data.frame")
    return(out)
  }
}


#Sum of cumulative incidence
GZ.SCI <- function(id, time, status, Tstart=0){
  
  ##Note: this function won't remove duplicated records 
  ##id: subject id
  ##time: individual time to event, time to competing risk, time to censor
  ##status: status at individual ending time (0=censor,1=event,2=competing risk)
  ##Tstart: age at study entry (to solve left-truncated problem)
  
  indata <- data.frame(id=id,
                       time=time,
                       status=status,
                       tstart=Tstart) %>% 
     mutate(time=ifelse(time==tstart,time+exp(-13),time)) #if event time = cohort entry, add a small amount 
  
  calc.trunc <- any(indata$tstart != 0)
  
  indata2 <- indata %>%
    arrange(id,time) %>% 
    group_by(id) %>%
    mutate(row_num = row_number(),
           n=n()) %>%
    ungroup()
  
  Event <- indata2 %>% filter(status==1) %>% tally()
  ftime <- sort(unique(indata2$time))
  
  if(Event$n==0){
    MCC.out <- data.frame(time=as.character(ftime),MCC=rep(0,length(ftime)))
  } else if(Event$n>0){
    
    M <- indata2 %>% filter(status==1) %>% group_by(id) %>% tally() %>% pull(n) %>% max() #Max num.of events per id
    CumI.list <- list()
    
    for (i in 1:M){
      
      in.i.th <- indata2 %>% filter(row_num==i)
      rest <- indata2 %>% filter(n<i & row_num==n) %>% 
        mutate(status=ifelse(status==1,0,status)) #because it contributed in data (i-1) already
      
      out.i.th <- rbind(in.i.th,rest) %>% select(-row_num,-n) 
      
      if (calc.trunc){
        #crprep2: convert data from a short format into a counting process format with
        #time varying weights. These weights correct for right-censored and left-truncated data,
        #allowing for analyses based on the subdistribution hazard
        
        count.data <- crprep2(Tstop="time", status="status", data=out.i.th, trans = 1, cens = 0, Tstart="tstart", id="id")
        fit.i.th <- survfit(Surv(Tstart,Tstop,status==1)~1,data=count.data,weight=weight.cens*weight.trunc)
        output <- summary(fit.i.th)
        output2<- data.frame(time=output$time,CumI=1-output$surv)
        CumI.i.th <- data.frame(time=ftime) %>% 
          left_join(output2,by='time') %>% 
          mutate(CumI=ifelse(row_number()==1 & is.na(CumI),0,CumI),
                 CumI=zoo::na.locf(CumI),
                 M=i)
        
      } else if (!calc.trunc){
        
        fit.i.th <- with(out.i.th,cmprsk::cuminc(time,status))
        CumI.i.th <- cmprsk::timepoints(fit.i.th,ftime)$est[1,] #row 1: event of interest
        
      }
      
      CumI.list[[i]] <- CumI.i.th
      
    }
    
    MCC.raw <- do.call(rbind,CumI.list)
    
    if(!calc.trunc){
      
      MCC.fill <- apply(MCC.raw, 1, zoo::na.locf, na.rm = FALSE)#To fill with previous MCC if it is NA  
      MCC <- rowSums(MCC.fill)
      MCC.out <- data.frame(MCC) %>% mutate(time=colnames(MCC.raw))
      
      
    } else if(calc.trunc){
      
      MCC.fill <- MCC.raw %>% tidyr::pivot_wider(names_from=M,values_from=CumI,names_prefix='M')
      MCC <- rowSums(MCC.fill[,-1])
      MCC.out <- data.frame(MCC) %>% mutate(time=MCC.fill$time) %>% select(time,MCC)
      
    }
  }
  
  MCC.out <- MCC.out %>% mutate_at('time',as.numeric)
  rownames(MCC.out) <- seq(nrow(MCC.out))
  
  return(MCC.out)
  
}

#Sum of cumulative incidence function with bootstrap 95% confidence interval
GZ.SCI.95CI <- function(id, time, status, Tstart, niter=NULL){
  
  ##id: subject id
  ##time: individual time to event, time to competing risk, time to censor
  ##status:status at individual ending time (0=censor,1=event,2=competing risk)
  ##niter=number of bootstrap iterations
  ##Tstart: age at study entry (to solve left-truncated problem)
  
  MCC.out <- GZ.SCI(id, time, status, Tstart)
  
  MCC.event <- MCC.out %>% arrange(time) %>% group_by(MCC) %>% slice(1) %>% ungroup()
  
  indata <- data.frame(id=id,
                       time=time,
                       status=status,
                       tstart=Tstart)
  
  uid <- unique(id)
  seed <- 2016
  
  for (boot in 1:niter){
    
    if (boot %% 100==0) cat(paste0("iteration: ", boot, "\n"))
    
    set.seed(boot+seed)
    sid <- sample(uid,replace = TRUE) 
    sid2 <- data.frame(new.id=seq(1,length(sid),1),id=sid)  
    bootdata <- merge(indata,sid2,by="id") %>% mutate(id=new.id)
    boot.out <- GZ.SCI(id=bootdata$id,time=bootdata$time,status=bootdata$status,Tstart=bootdata$tstart) 
    boot.out2 <- rename_with(boot.out,
                             ~paste0(.x,boot,recycle0=T),
                             starts_with('MCC'))
    
    MCC.event <- MCC.event %>% left_join(boot.out2,by="time") #only capture CI for time point that MCC changed
    
  }
  
  ##If the first row is missing, then MCC=0.
  MCC.event[1,][is.na(MCC.event[1,])] <- 0
  ##For points not in a bootstrapped data, MCC is NA. Need to fill in by its previous MCC value
  MCC.fill <- MCC.event %>% mutate_all(zoo::na.locf)
  
  quantiles <- function(x) {
    return(quantile(x, probs = c(0.025, 0.975)))
  }
  
  MCC.result <- t(apply(MCC.fill[,-c(1:2)], 1, quantiles))
  MCC.result2 <- cbind(MCC.fill[,c(1:2)], MCC.result)
  MCC.output <- MCC.out %>% 
    left_join(MCC.result2,by=c('time','MCC')) %>% 
    group_by(MCC) %>% 
    reframe(time=time,
            lci=zoo::na.locf(`2.5%`),
            uci=zoo::na.locf(`97.5%`)) %>% 
    ungroup() %>% 
    select(time,MCC,lci,uci)
  
  MCC.95CI.list <- list('MCC.output'=MCC.output,
                        'MCC.fill'=MCC.fill)
  
  return(MCC.95CI.list)
  
}

#Wrap up
MY.MCC <- function(id, time, status, Tstart, ci=TRUE, niter=1000){
  
  if (ci){
    
    out.95CI.list <- GZ.SCI.95CI(id, time, status, Tstart, niter) 
    out.95CI.list[['MCC.output']] %>% return()
    
  } else {
    
    GZ.SCI(id,time,status,Tstart) %>% return()
  }
}
