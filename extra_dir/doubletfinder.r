RunDoubletFinder <- function(
	object,
	PCs = 1:30,
	doublet.rate = 0.06,
	db.ratio = 0.25,
	nb.size = 0.09,
	sct = FALSE,
	GT = FALSE,
	GT.calls = NULL,
	identity = NULL,
	...
){
	#Set pN-pK param sweep ranges
	pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01)) # 34
	pN <- seq(0.05,0.3,by=0.05)

	# Remove pK values with too few cells
	min.cells <- round(ncol(object)/(1-0.05) - ncol(object) )
	pK.test <- round(pK*min.cells)
	pK <- pK[which(pK.test >= 1)]

	# Down-sample cells to 10000 (when applicable) for computational effiency
	data <- SeuratObject::GetAssayData(object, slot = "counts")
	if ( ncol(object) > 10000 ) {
		real.cells <- colnames(object)[sample(1:ncol(object), 10000, replace=FALSE)]
		data <- data[ , real.cells]
	}else{
		real.cells <- colnames(object)
	}

	# Iterate through pN, computing pANN vectors at varying pK
	n.real.cells <- length(real.cells)
	output2 <- future.apply::future_lapply(pN, function(x){
		# Make merged real-artifical data
		print(paste("Creating artificial doublets for pN = ", x*100,"%",sep=""))
		n_doublets <- round(n.real.cells/(1 - x) - n.real.cells)
		real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
		real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
		doublets <- (data[, real.cells1] + data[, real.cells2])/2
		colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
		data_wdoublets <- cbind(data, doublets)

		## Pre-process Seurat object
		seu_wdoublets <- SeuratObject::CreateSeuratObject(counts = data_wdoublets)
		if (sct == FALSE){
			seu_wdoublets <- Seurat::NormalizeData(seu_wdoublets)
			seu_wdoublets <- Seurat::FindVariableFeatures(seu_wdoublets)
			seu_wdoublets <- Seurat::ScaleData(seu_wdoublets)
		}else{
			seu_wdoublets <- Seurat::SCTransform(seu_wdoublets)
		}

		seu_wdoublets <- Seurat::RunPCA(seu_wdoublets, verbose=FALSE)
		## Compute PC distance matrix
		nCells <- ncol(seu_wdoublets)
		pca.coord <- SeuratObject::Embeddings(seu_wdoublets, reduction = "pca")[,PCs]
		# pca.coord <- object_wdoublets@reductions$pca@cell.embeddings[ , PCs]
		rm(seu_wdoublets);gc()
		dist.mat <- fields::rdist(pca.coord)[,1:n.real.cells]

		# Pre-order PC distance matrix prior to iterating across pK for pANN computations
		for (i in 1:n.real.cells) {
			dist.mat[,i] <- order(dist.mat[,i])
		}

		# Trim PC distance matrix for faster manipulations
		ind <- round(nCells * max(pK))+5
		dist.mat <- dist.mat[1:ind, ]

		## Compute pANN across pK sweep
		print("Computing pANN across all pK...")
		sweep.res.list = list()
		list.ind = 0
		for (k in 1:length(pK)) {
			print(paste("pK = ", pK[k], "...", sep = ""))
			pk.temp <- round(nCells * pK[k])
			pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, ncol = 1))
			colnames(pANN) <- "pANN"
			rownames(pANN) <- real.cells
			list.ind <- list.ind + 1

			for (i in 1:n.real.cells) {
				neighbors <- dist.mat[2:(pk.temp + 1),i]
				pANN$pANN[i] <- length(which(neighbors > n.real.cells))/pk.temp
			}

			sweep.res.list[[list.ind]] <- pANN
		}
		sweep.res.list
	}, future.seed = 2020)

	## Write parallelized output into list
	sweep.list <- list()
	list.ind <- 0

	for(i in 1:length(output2)){
		for(j in 1:length(output2[[i]])){
			list.ind <- list.ind + 1
			sweep.list[[list.ind]] <- output2[[i]][[j]]
		}
	}

	## Assign names to list of results
	name.vec <- NULL
	for (j in 1:length(pN)) {
		name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
	}
	names(sweep.list) <- name.vec

	## Set pN-pK param sweep ranges
	# require(KernSmooth); require(ROCR)
	# name.vec <- names(sweep.list)
	name.vec <- unlist(strsplit(name.vec, split="pN_"))
	name.vec <- name.vec[seq(2, length(name.vec), by=2)]
	name.vec <- unlist(strsplit(name.vec, split="_pK_"))
	pN <- as.numeric(unique(name.vec[seq(1, length(name.vec), by=2)]))
	pK <- as.numeric(unique(name.vec[seq(2, length(name.vec), by=2)]))
	print(paste("PK ==+ ", pK, "..."))
	print(paste("PN ==+ ", pN, "..."))
	## Initialize data structure w/ or w/o AUC column, depending on whether ground-truth doublet classifications are available
	if ( GT == TRUE ) {
		sweep.stats <- as.data.frame(matrix(0L, nrow=length(sweep.list), ncol=4))
		colnames(sweep.stats) <- c("pN","pK","AUC","BCreal")
		sweep.stats$pN <- factor(rep(pN, each=length(pK), levels = pN))
		sweep.stats$pK <- factor(rep(pK, length(pN),levels = pK))
	}else{
		sweep.stats <- as.data.frame(matrix(0L, nrow=length(sweep.list), ncol=3))
		colnames(sweep.stats) <- c("pN","pK","BCreal")
		sweep.stats$pN <- factor(rep(pN, each=length(pK), levels = pN))
		sweep.stats$pK <- factor(rep(pK, length(pN),levels = pK))
	}

	## Perform pN-pK parameter sweep summary
	for (i in 1:length(sweep.list)) {
		res.temp <- sweep.list[[i]]

		## Use gaussian kernel density estimation of pANN vector to compute bimodality coefficient
		gkde <- approxfun(KernSmooth::bkde(res.temp$pANN, kernel="normal"))
		x <- seq(from=min(res.temp$pANN), to=max(res.temp$pANN), length.out=nrow(res.temp))
		gkdx = gkde(x)
		sweep.stats$BCreal[i] <- bimodality_coefficient(gkdx)

		if (GT == FALSE) { next }

		## If ground-truth doublet classifications are available, perform ROC analysis on logistic
		## regression model trained using pANN vector
		meta <- as.data.frame(matrix(0L, nrow=nrow(res.temp), ncol=2))
		meta[,1] <- GT.calls
		meta[,2] <- res.temp$pANN
		train.ind <- sample(1:nrow(meta), round(nrow(meta)/2), replace=FALSE)
		test.ind <- (1:nrow(meta))[-train.ind]
		colnames(meta) <- c("SinDub","pANN")
		meta$SinDub <- factor(meta$SinDub, levels = c("Doublet","Singlet"))
		model.lm <- glm(SinDub ~ pANN, family="binomial"(link='logit'), data=meta, subset=train.ind)
		prob <- predict(model.lm, newdata=meta[test.ind, ], type="response")
		ROCpred <- ROCR::prediction(predictions=prob, labels=meta$SinDub[test.ind])
		perf.auc <- ROCR::performance(ROCpred, measure="auc")
		sweep.stats$AUC[i] <- perf.auc@y.values[[1]]
	}

	## Implementation for data without ground-truth doublet classifications
	if ( !"AUC" %in% colnames(sweep.stats) ) {
		## Initialize data structure for results storage
		bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=5))
		colnames(bc.mvn) <- c("ParamID","pK","MeanBC","VarBC","BCmetric")
		bc.mvn$pK <- unique(sweep.stats$pK)
		bc.mvn$ParamID <- 1:nrow(bc.mvn)

		## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
		x <- 0
		for (i in unique(bc.mvn$pK)) {
			x <- x + 1
			ind <- which(sweep.stats$pK == i)
			bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
			bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
			bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
		}
	}else{
		# Implementation for data with ground-truth doublet classifications (e.g., MULTI-seq, CellHashing, Demuxlet, etc.)
		# Initialize data structure for results storage
		bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=6))
		colnames(bc.mvn) <- c("ParamID","pK","MeanAUC","MeanBC","VarBC","BCmetric")
		bc.mvn$pK <- unique(sweep.stats$pK)
		bc.mvn$ParamID <- 1:nrow(bc.mvn)

		## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
		x <- 0
		for (i in unique(bc.mvn$pK)) {
			x <- x + 1
			ind <- which(sweep.stats$pK == i)
			bc.mvn$MeanAUC[x] <- mean(sweep.stats[ind, "AUC"])
			bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
			bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
			bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
		}
	}

	# choose parameters
	maxBCmetric    <- max(bc.mvn$BCmetric, na.rm = TRUE)
	pK <- as.numeric(as.character(bc.mvn[bc.mvn$BCmetric==maxBCmetric, ]$pK))

	# compute doublet scores
	annotations <- object@meta.data$seurat_clusters  ## ex: annotations <- object_kidney@meta.data$ClusteringResults
	# homotypic.prop <- modelHomotypic(annotations)
	anno.freq <- table(annotations)/length(annotations)
	homotypic.prop <- sum(anno.freq^2)
	nExp_poi       <- round(doublet.rate*ncol(object))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
	nExp   <- round(nExp_poi*(1-homotypic.prop))
	# object.scored     <- doubletFinder_v3(object, PCs =PCs, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = use.SCT)
	print(paste( "pN =", pN, "pK=", pK, "nExp=", nExp))
	## Generate new list of doublet classificatons from existing pANN vector to save time
	## Make merged real-artifical data
	real.cells <- colnames(object)
	data <- SeuratObject::GetAssayData(object, slot = "counts")[, real.cells]
	n_real.cells <- length(real.cells)
	n_doublets <- round(n_real.cells/(1 - db.ratio) - n_real.cells)
	print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
	real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
	real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
	doublets <- (data[, real.cells1] + data[, real.cells2])/2
	colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
	data_wdoublets <- cbind(data, doublets)

	## Pre-process Seurat object
	print("Creating Seurat object...")
	seu_wdoublets <- SeuratObject::CreateSeuratObject(counts = data_wdoublets)
	if (sct == FALSE){
		print("Normalizing Seurat object...")
		seu_wdoublets <- Seurat::NormalizeData(seu_wdoublets)
		print("Finding variable genes...")
		seu_wdoublets <- Seurat::FindVariableFeatures(seu_wdoublets)
		print("Scaling data...")
		seu_wdoublets <- Seurat::ScaleData(seu_wdoublets)
	}else{
		print("Running SCTransform...")
		seu_wdoublets <- Seurat::SCTransform(seu_wdoublets)
	}
	print("Running PCA...")
	seu_wdoublets <- Seurat::RunPCA(seu_wdoublets, verbose=FALSE)
	## Compute PC distance matrix
	print("Calculating PC distance matrix...")
	nCells <- ncol(seu_wdoublets)
	pca.coord <- SeuratObject::Embeddings(seu_wdoublets, reduction = "pca")[,PCs]
	k <- round(ncol(seu_wdoublets) * nb.size)
	rm(seu_wdoublets);gc()

	## Compute PC distance matrix
	print("Calculating PC distance matrix...")
	dist.mat <- fields::rdist(pca.coord)

	## Compute pANN
	print("Computing pANN...")
	pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
	rownames(pANN) <- real.cells
	colnames(pANN) <- "pANN"
	for (i in 1:n_real.cells) {
		neighbors <- order(dist.mat[, i])
		neighbors <- neighbors[2:(k + 1)]
		neighbor.names <- rownames(dist.mat)[neighbors]
		pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
	}

	print("Classifying doublets..")
	classifications <- rep("FALSE",n_real.cells)
	classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "TRUE"
	object <- SeuratObject::AddMetaData(object, metadata = pANN[colnames(object), 1], col.name = paste("pANN", db.ratio, nb.size, nExp, sep="_"))
	object <- SeuratObject::AddMetaData(object, metadata = classifications, col.name = "predicted_doublets")
	return(object)
}

bimodality_coefficient <- function(x) {
      n <- length(x)
      S <- (1/n)*sum((x-mean(x))^3)/(((1/n)*sum((x-mean(x))^2))^1.5)
      G <- S*(sqrt(n*(n-1)))/(n-2)
      K <- (1/n)*sum((x-mean(x))^4)/(((1/n)*sum((x-mean(x))^2))^2)-3
      K <- ((n - 1)*((n+1)*K-3*(n-1))/((n-2)*(n-3)))+3
      B <- ((G^2)+1)/(K+((3*((n-1)^2))/((n-2)*(n-3))))
      return(B)
    }



