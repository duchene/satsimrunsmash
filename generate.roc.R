# This function takes two vectors of values. One is a binary variable determining whether each case is a positive or not, while the other is a variable of statistic that aims to identify these positives. Assumes that statistics are mostly positive.

generate.roc <- function(positives, stat, Nthresh = 100){
	nonnas <- which(!is.na(stat))
	positives <- positives[nonnas]
	stat <- stat[nonnas]
	statran <- range(stat, na.rm = T)
	threshs <- seq(statran[1], statran[2], length.out = Nthresh)
	rocdata <- matrix(NA, Nthresh, 2)
	pos <- sum(positives)
	neg <- sum(!positives)
	
	for(i in 1:Nthresh){
		rocdata[i, 1] <- sum(!positives[which(stat < threshs[i])]) / neg
		rocdata[i, 2] <- sum(positives[which(stat < threshs[i])]) / pos
	}
	rocdata <- cbind(rocdata, rocdata[,2] - rocdata[,1])
	rownames(rocdata) <- threshs
	colnames(rocdata) <- c("FPR", "TPR", "subt")
	return(rocdata)
}

# The following function generalizes for any number of positives vectors and any number of statistics

general.roc.plotting <- function(positives.vectors, stats.vectors, dat){
	rocslist <- list()
	#dev.new()
	#par(mfrow = c(1, length(positives.vectors)))

	for(i in 1:length(positives.vectors)){
	      thresholds <- bestsubtraction <- tprs <- fprs <- vector()
	      for(j in 1:length(stats.vectors)){
		    roctab <- try(generate.roc(dat[!is.na(dat[,positives.vectors[i]]),positives.vectors[i]], dat[!is.na(dat[,positives.vectors[i]]),stats.vectors[j]]))
		    if(class(roctab) == "try-error") next
		    if(j == 1){
		    	 plot(roctab[,1:2], type = "l", lwd = 2, main = paste(positives.vectors[i], "ROC"))
		    } else {
		      	 lines(roctab[,1:2], col = j)
		    }
		    indbest <- which(roctab[,3] == max(roctab[,3]))[1]
		    thresholds[length(thresholds) + 1] <- round(as.numeric(rownames(roctab)[indbest]), 3)
		    bestsubtraction[length(bestsubtraction) + 1] <- round(roctab[indbest,3], 3)
		    tprs[length(tprs) + 1] <- round(roctab[indbest,2], 3)
		    fprs[length(fprs) + 1] <- round(roctab[indbest,1], 3)
		    rocslist[[length(rocslist) + 1]] <- roctab
		    names(rocslist)[length(rocslist)] <- names(thresholds)[length(thresholds)] <- names(bestsubtraction)[length(bestsubtraction)] <- names(tprs)[length(tprs)] <- names(fprs)[length(fprs)] <- paste0(positives.vectors[i], "_", stats.vectors[j])
	      	    if(j < length(stats.vectors)) next
		    legend(x = "bottomright", legend = c(paste("Best TPR-FPR: ", paste(bestsubtraction, collapse = ", ")), paste("Thresh: ", paste(thresholds, collapse = ", "))))
	      	    abline(0, 1, lty = 2)
	      }
	      
	}
	
	results <- list(roc.tables = rocslist, best.subtraction = bestsubtraction, thresholds = thresholds, thres.tpr = tprs, thres.fpr = fprs)
	
	return(results)
}