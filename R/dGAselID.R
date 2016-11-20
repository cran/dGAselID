####################################################################
#	InitialPopulation
####################################################################
#
InitialPopulation <- function(x, populationSize, startGenes, EveryGeneInInitialPopulation=TRUE) {
  cat(paste("\nGenerating the initial population...\n"))
  genomeLength=length(featureNames(x));
  cat("Generating random chromosomes...\n");
  populationInit = matrix(0, nrow=populationSize, ncol=genomeLength);

  if (EveryGeneInInitialPopulation==TRUE){
    noDGenes=ceiling(genomeLength/populationSize);
    index1<-sample(1:genomeLength, genomeLength, replace=FALSE)

    if (genomeLength%%populationSize==0){
      indexI<-matrix(index1, nrow=populationSize)
      cat(paste("\nNo rest...\n"))
    } else {
      cat(paste("\nWith rest...\n"))
      index2<-sample(setdiff(1:genomeLength,tail(index1, startGenes+1)), abs(populationSize*noDGenes-genomeLength), replace=FALSE)
      indexI<-matrix(c(index1,index2), nrow=populationSize, byrow = TRUE)
    }

    for (i in 1:populationSize) {
      for(j in 1:noDGenes){
        populationInit[i, indexI[i, j]] = 1;
      }
    }

    if (startGenes>noDGenes){
      for (i in 1:populationSize) {
        populationInit[i, sample(setdiff(1:genomeLength,indexI[i,]), startGenes-noDGenes, replace=FALSE)] = 1;
      }
    }
  } else {
    for (i in 1:populationSize) {
      populationInit[i, sample(1:genomeLength, startGenes, replace=FALSE)] = 1;
    }
  }

  colnames(populationInit)=featureNames(x);
  rownames(populationInit)=1:populationSize;
  return(populationInit);
}
####################################################################
#	Individuals
####################################################################
#
Individuals<-function(population){
  cat(paste("\tGenerating Individuals...\n"))
  if(nrow(population)%%2==1){
    population<-rbind(population, population[1,])
    rownames(population)<-c(1:nrow(population))
  }

  noIndividuals<-nrow(population)
  Id<-rep(1:(nrow(population)/2), each=2)
  population<-cbind(Id, population)
}

####################################################################
#	Split Chromosomes
####################################################################
#
splitChromosomes <- function(x, noChr=22) {

  noGenes<-matrix(c(1, 3000, 9.174312, 2, 2500, 7.645260, 3, 1900, 5.810398, 4, 1600, 4.892966, 5, 1700, 5.198777, 6, 1900, 5.810398, 7, 1800, 5.504587, 8, 1400, 4.281346, 9, 1400, 4.281346, 10, 1400, 4.281346, 11, 2000, 6.116208, 12, 1600, 4.892966, 13, 800, 2.446483, 14, 1200, 3.669725, 15, 1200, 3.669725, 16, 1300, 3.975535, 17, 1600, 4.892966, 18, 600, 1.834862, 19, 1700, 5.198777, 20, 900, 2.752294, 21, 400, 1.223242, 22, 800, 2.446483), nrow=22, ncol=3, byrow=TRUE)
  colnames(noGenes)<-c("Chromosome", "NoOfGenes",  "Percent")

  toSplit<-dim(exprs(x))[1];

  if (3*noChr>toSplit){
    cat(paste("\nToo many chromosomes for the given genome. Please specify a lower number of chromosomes...\n"));
    noChr<-toSplit%/%3
    cat(paste("\nAutomatically changing the number of chromosomes to "), noChr, paste(" ...\n"));
  }

  noGenes2<-noGenes[1:noChr,3];
  rm(noGenes);

  newPercent<-noGenes2/sum(noGenes2)*100;
  rm(noGenes2);

  chr<-as.integer(newPercent*toSplit/100);
  rm(newPercent);

  chr[noChr]<-chr[noChr]+(toSplit-sum(chr));  #Adjust for the length difference
  rm(toSplit);
  rm(noChr);

  chrNumber<-0;
  chrConfig<-c();
  for (i in chr) {
    chrNumber<-chrNumber+1
    chrConfig<-c(chrConfig, rep(chrNumber,i));
  }
  rm(chrNumber);
  rm(chr);
  rm(i);
  return(chrConfig);
}

####################################################################
#	RandomizePop
####################################################################
#
RandomizePop<-function(population){
  cat(paste("\tRandomizing the population...\n"))
  newIndex<-sample(1:dim(population)[1], replace=FALSE)
  newPopulation<-population[newIndex,]
  rownames(newPopulation)<-c(1:dim(newPopulation)[1])
  return(newPopulation)
}

####################################################################
#	Evaluate Population
####################################################################
#
EvaluationFunction <- function(x, individuals, response, method, trainTest, nnetSize=NA, nnetDecay=NA, rdaAlpha=NA, rdaDelta=NA, ...){
  cat(paste("\tEvaluating Fitnesses...\n"))
  if(toString(trainTest)=="LOO"){
    traintest<-xvalSpec("LOO")
  } else if (toString(trainTest)=="LOG"){
    traintest<-xvalSpec("LOG", 5, balKfold.xvspec(5))
  } else {
    traintest<-c(trainTest)
  }
  formula<-as.formula(paste(response, "~."))
  population<-individuals[,2:dim(individuals)[2]]
  populationSize<-nrow(population)
  results<-c()

  for (q in 1:populationSize) {
    if (!is.na(nnetSize)) {
      result = MLearn(formula, x[population[q,]==1,], .method = method, trainInd = traintest, size=nnetSize, decay=nnetDecay)
    } else if (!is.na(rdaAlpha)) {
      result = MLearn(formula, x[population[q,]==1,], .method = method, trainInd = traintest, alpha=rdaAlpha, delta=rdaDelta)
    } else {
      result = MLearn(formula, x[population[q,]==1,], .method = method, trainInd = traintest)
    }
    Accuracy<-(confuMat(result)[1]+confuMat(result)[4])/(confuMat(result)[1]+confuMat(result)[2]+confuMat(result)[3]+confuMat(result)[4])
    results<-rbind(results, Accuracy)
  }
  results<-cbind(individuals[,1], results)
  rownames(results)=1:populationSize

  avgs<-c()
  for (j in results[,1]){
    avgs<-rbind(avgs, mean(results[results[,1]==j,2]))
  }
  results<-cbind(results, avgs)
  colnames(results)=c("Id", "Accuracy", "Average Acc")
  return(results)
}

####################################################################
#	AnalyzeResults
####################################################################
#
AnalyzeResults<-function(individuals, results, randomAssortment=TRUE, chrConf){
  cat(paste("\tAnalyzing the results\n"))
  keep<-matrix(0, nrow=nrow(individuals), ncol=1)
  rownames(keep)<-rownames(results)
  colnames(keep)<-"Keep"
  crossOvers<-matrix(c(0), nrow=0, ncol=ncol(individuals))
  colnames(crossOvers)<-colnames(individuals)
  for (i in 1:nrow(results)){
    if (results[i, 2] >= results[i, 3])
      keep[i, 1]=1
  }

  cat(paste("\tApplying crossovers...\n"))

  if (randomAssortment==TRUE){
    for (j in 1:(nrow(individuals)/2)){
      selChr<-individuals[results[,1]==j, -1]
      newChr<-RandomAssortment(Crossover(selChr[1,], selChr[2,], chrConf), chrConf)
      newChromosomes<-cbind(c(j,j), newChr)
      crossOvers<-rbind(crossOvers, newChromosomes)
    }

    cat(paste("\tApplying Random Assortment...\n"))

    rownames(crossOvers)<-1:nrow(crossOvers)

  } else {
    for (j in 1:(nrow(individuals)/2)){
      selChr<-individuals[results[,1]==j, -1]
      newChr<-Crossover(selChr[1,], selChr[2,], chrConf)
      newChromosomes<-cbind(c(j,j), newChr)
      crossOvers<-rbind(crossOvers, newChromosomes)
    }

    rownames(crossOvers)<-1:nrow(crossOvers)
  }

  return(list(keep, crossOvers))
}

####################################################################
#	Crossover
####################################################################
#
Crossover<-function(c1, c2, chrConf){
  crossVector<-rep(0, length(c1))
  for (i in 1:max(chrConf)) {
    crossSlice<-crossVector[chrConf==i]
    crossIndexes<-sort(sample(1:(length(crossSlice)),2, replace=TRUE))
    crossSlice[crossIndexes[1]:crossIndexes[2]]=1
    rm(crossIndexes)
    crossVector[chrConf==i]<-crossSlice
  }

  c3<-rep(NA, length(c1))
  c4<-rep(NA, length(c1))
  maTemp<-rbind(chrConf, crossVector,c1, c2, c3, c4)
  rm(c3)
  rm(c4)
  rm(crossVector)

  maTemp[5,maTemp[2,]==0]<-maTemp[3,maTemp[2,]==0]
  maTemp[5,maTemp[2,]==1]<-maTemp[4,maTemp[2,]==1]
  maTemp[6,maTemp[2,]==0]<-maTemp[4,maTemp[2,]==0]
  maTemp[6,maTemp[2,]==1]<-maTemp[3,maTemp[2,]==1]
  newChrs<-maTemp[5:6,]
  rm(maTemp)

  return(newChrs)
}
####################################################################
#	RandomAssortment
####################################################################
#
RandomAssortment<-function(newChrs, chrConf){
  AssortIndex<-sample(c(0,1), size=max(chrConf), replace = TRUE)
  c3<-rep(NA, ncol(newChrs))
  c4<-rep(NA, ncol(newChrs))
  exchange<-c()
  for (i in 1:max(chrConf)) {
    exchange<-c(exchange, rep(AssortIndex[i], sum(as.numeric(chrConf==i))))
  }

  maTemp<-rbind(chrConf, exchange, newChrs[1,], newChrs[2,], c3, c4)
  rm(c3)
  rm(c4)
  rm(AssortIndex)
  rm(exchange)
  maTemp[5,maTemp[2,]==0]<-maTemp[3,maTemp[2,]==0]
  maTemp[5,maTemp[2,]==1]<-maTemp[4,maTemp[2,]==1]
  maTemp[6,maTemp[2,]==0]<-maTemp[4,maTemp[2,]==0]
  maTemp[6,maTemp[2,]==1]<-maTemp[3,maTemp[2,]==1]
  splittedChrs<-maTemp[5:6,]
  rm(maTemp)
  return(splittedChrs)
}

####################################################################
#	Mutation
####################################################################
#
Mutation<-function(individuals, mutationChance){
  noMutations<-floor(mutationChance/100*length(individuals[,-1]))
  cat(paste("\tApplying", noMutations, "mutations...\n"))
  indexes<-sample(0:length(individuals[,-1]),noMutations)
  individuals[,-1][indexes]=as.numeric(!as.logical(individuals[,-1][indexes]))
  return(individuals)
}

####################################################################
#	Elitism
####################################################################
#
Elitism<-function(results, elitism, ID){
  cat(paste("\tApplying Elitism...Keeping the Best ", elitism, "%\n"))
  if (ID=="ID2"){
    cat("Elitistic individuals...\n")
    elite<-3
  } else {
    cat("Elitistic genotypes...\n")
    elite<-2
  }

  keep<-matrix(0, nrow=nrow(results), ncol=1)
  rownames(keep)<-rownames(results)
  colnames(keep)<-"Keep"
  toKeep<-sort(results[,elite], decreasing = TRUE, index.return=TRUE)	#aici pune results[,2]
  toKeep<-toKeep$ix
  newIndex<-floor(length(toKeep)*elitism/100)
  rez<-results[toKeep,]
  toKeep<-toKeep[1:newIndex]
  keep[toKeep]<-1
  return(list(keep, toKeep))
}
####################################################################
#	EmbryonicSelection
####################################################################
#
EmbryonicSelection<-function(population, results, embryonicSelection){
  cat(paste("\tApplying Embryonic Selection for Fitness > ", embryonicSelection, "\n"))
  keep<-matrix(0, nrow=nrow(population), ncol=1)
  rownames(keep)<-rownames(results)
  colnames(keep)<-"Keep"

  keep[(results[,3] > embryonicSelection)==TRUE]=1

  return(keep)
}

####################################################################
#	PlotGenAlg
####################################################################
#
PlotGenAlg <- function(DGenes, dGenes, maxEval, meanEval){
  dev.off()
  dev.new()
  setFr <- layout(matrix(c(1,2,3,3),2,2,byrow = TRUE), TRUE)
  layout.show(setFr)
  par(las=2)
  index<-sort(DGenes, decreasing=TRUE, index.return=TRUE)
  DGenes<-DGenes[index$ix]
  plottedGenes<-length(DGenes)
  plot(maxEval, type="o", col="red", xlab="iteration no.", main="Maximum Accuracy")
  plot(meanEval, type="o", col="blue", xlab="iteration no.", main="Mean Accuracy")
  barplot(DGenes[1:plottedGenes], main="Genes inheritance", xlab="Gene", col=c("darkblue"), beside=FALSE, add=FALSE)
}

####################################################################
#	dGAselID
####################################################################
#
#'  dGAselID
#'
#' Initializes and starts the search with the genetic algorithm.
#' @aliases dGAselID
#' @param x Dataset in ExpressionSet format.
#' @param response Response variable
#' @param method Supervised classifier for fitness evaluation. Most of the supervised classifiers in MLInterfaces are acceptable. The default is knnI(k=3, l=2).
#' @param trainTest Cross-validation method. The default is "LOG".
#' @param startGenes Genes in the genotypes at initialization.
#' @param populationSize Number of genotypes in initial population.
#' @param iterations Number of iterations.
#' @param noChr Number of chromosomes. The default value is 22.
#' @param elitism Elite population in percentages.
#' @param ID Incomplete Dominance. The default value is "ID1". Use "ID2" for elitism applied to individuals.
#' @param pMutationChance Chance for a point mutation to occur.
#' @param randomAssortment Random Assortment of Chromosomes for recombinations. The default value is TRUE.
#' @param embryonicSelection Remove chromosomes with fitness < specified value. The default value is NA.
#' @param EveryGeneInInitialPopulation Request for every gene to be present in the initial population. The default value is TRUE.
#' @param nnetSize for nnetI. The default value is NA.
#' @param nnetDecay for nnetI. The default value is NA.
#' @param rdaAlpha for rdaI. The default value is NA.
#' @param rdaDelta for rdaI. The default value is NA.
#' @param ... Additional arguments.
#' @examples
#' \dontrun{
#'  library(genefilter)
#'  library(ALL)
#'  data(ALL)
#'  bALL = ALL[, substr(ALL$BT,1,1) == "B"]
#'  smallALL = bALL[, bALL$mol.biol %in% c("BCR/ABL", "NEG")]
#'  smallALL$mol.biol = factor(smallALL$mol.biol)
#'  smallALL$BT = factor(smallALL$BT)
#'  f1 <- pOverA(0.25, 9)
#'  f2 <- function(x) (IQR(x) > 0.75)
#'  f3 <- ttest(smallALL$mol.biol, p=0.1)
#'  ff <- filterfun(f1, f2, f3)
#'  selectedsmallALL <- genefilter(exprs(smallALL), ff)
#'  sum(selectedsmallALL)
#'
#'  set.seed(1357)
#'  res1<-dGAselID(smallALL, "mol.biol", method=knn.cvI(k=3, l=2), trainTest=1:79,
#'    startGenes=5, populationSize=50, iterations=4, noChr=5, mutationChance=0.05,
#'    elitism=10, ID="ID1")
#'  }
#' @export



dGAselID<-function(x, response, method=knnI(k=3, l=2), trainTest="LOG", startGenes, populationSize, iterations, noChr=22, elitism=NA, ID="ID1", pMutationChance=NA, randomAssortment=TRUE, embryonicSelection=NA, EveryGeneInInitialPopulation=TRUE, nnetSize=NA, nnetDecay=NA, rdaAlpha=NA, rdaDelta=NA, ...){

  if (typeof(x)!="S4") {
    stop("The supplied data is not an ExpressionSet.");
  }

  if(randomAssortment==TRUE){
    cat("The chromosomes will be randomly assigned...\n")
  }

  if (EveryGeneInInitialPopulation==TRUE){
    cat("Every gene will be present in the initial population...\n")
  }

  if (is.na(embryonicSelection)) {
    embryonicSelection = NA
  }

  if (is.na(pMutationChance)) {
    pMutationChance = 1/(2*length(featureNames(x))+1)
  }

  if (is.na(elitism)) {
    elitism = 0
  }

  cat("Elitism =", elitism, "%\n")
  cat("Mutations rate =", pMutationChance, "%\n")
  cat("Embryonic Selection for fitness > ", embryonicSelection, "\n")
  cat("Fitness evaluation function =", method@mlFunName, "\n")
  cat("Cross-validation =", trainTest, "\n")

  cat("\nInitial population...\n")
  initialPopulation<-InitialPopulation(x, populationSize, startGenes, EveryGeneInInitialPopulation)
  individuals<-Individuals(initialPopulation)

  cat(paste("Splitting the genotype in", noChr, "chromosomes\n"))
  chrConf<-splitChromosomes(x, noChr)

  kDGenes<-matrix(c(0), nrow=1, ncol=ncol(initialPopulation))
  colnames(kDGenes)<-colnames(initialPopulation)
  rownames(kDGenes)<-"DGenes"

  kdGenes<-matrix(c(0), nrow=1, ncol=ncol(initialPopulation))
  colnames(kdGenes)<-colnames(initialPopulation)
  rownames(kdGenes)<-"dGenes"

  MaxAcc<-c()
  MeanAcc<-c()
  MinAcc<-c()
  bestIndividual<-matrix(0, nrow=1, ncol=ncol(initialPopulation))

  iteration<-0
  dev.new()

  repeat{
    cat(paste("Starting iteration no.", iteration, "\n"))

    results<-EvaluationFunction(x, individuals, response, method, trainTest, nnetSize, nnetDecay, rdaAlpha, rdaDelta)

    iterMinAccuracy <- range(results[,2])[1]
    MinAcc<-c(MinAcc, iterMinAccuracy)
    cat(paste("\tMinimum Fitness in iteration no.", iteration, "equals", iterMinAccuracy*100, "%\n"))

    iterMeanAccuracy <- mean(results[,2])
    MeanAcc<-c(MeanAcc, iterMeanAccuracy)
    cat(paste("\tMean Fitness in iteration no.", iteration, "equals", iterMeanAccuracy*100, "%\n"))

    iterMaxAccuracy <- range(results[,2])[2]
    MaxAcc<-c(MaxAcc, iterMaxAccuracy)
    cat(paste("\tMaximum Fitness in iteration no.", iteration, "equals", iterMaxAccuracy*100, "%\n"))

    lostInEmbryonic<-0
    if (!is.na(embryonicSelection)) {
      keptEmbr<-EmbryonicSelection(individuals, results, embryonicSelection)
      keptIndividuals<-individuals[keptEmbr==1,]
      keptIndividuals[,1]<-rep(1:(nrow(keptIndividuals)/2), each=2)
      rownames(keptIndividuals)<-1:nrow(keptIndividuals)
      keptResults<-results[keptEmbr==1,]
      keptResults[,1]<-rep(1:(nrow(keptIndividuals)/2), each=2)
      rownames(keptResults)<-1:nrow(keptResults)
      discardedIndividuals<-individuals[keptEmbr==0,]
      lostInEmbryonic<-nrow(discardedIndividuals)
    } else {
      keptIndividuals<-individuals
      keptResults<-results
      discardedIndividuals<-individuals[0,]
    }

    iterRes<-AnalyzeResults(keptIndividuals, keptResults, randomAssortment, chrConf)

    keptIndividualsFromChildren<-iterRes[[2]]

    discardedIndividuals<-rbind(discardedIndividuals, keptIndividuals[!(iterRes[[1]]==1),])
    keptIndividualsFromParents<-keptIndividuals[(iterRes[[1]]==1),]
    rownames(keptIndividualsFromParents)<-1:nrow(keptIndividualsFromParents)
    keptResults<-keptResults[iterRes[[1]]==1,]
    rownames(keptResults)<-1:nrow(keptResults)

    keptInElitism<-0
    if (floor(nrow(keptResults)*elitism/100)==0) {
      tempResults<-keptResults
      tempIndividualsFromParents<-keptIndividualsFromParents
      keptIndividualsFromParents<-keptIndividualsFromParents[0,]
      keptResults<-keptResults[0,]
      keptElit<-Elitism(tempResults, 50, ID)
      forBest<-tempIndividualsFromParents[keptElit[[2]],]
      forBest<-forBest[1,]
      rm(tempResults)
      rm(tempIndividualsFromParents)
    } else {
      keptElit<-Elitism(keptResults, elitism, ID)
      keptIndividualsFromParents<-keptIndividualsFromParents[keptElit[[2]],]
      forBest<-keptIndividualsFromParents
      forBest<-forBest[1,]
      keptResults<-keptResults[keptElit[[2]],]
      keptInElitism<-nrow(keptResults)
    }

    best<-t(as.matrix(forBest))[,-1]
    bestIndividual<-rbind(bestIndividual, best)
    rm(forBest)

    if(lostInEmbryonic<keptInElitism){
      adjust<-keptInElitism-lostInEmbryonic
      toRemove<-sample(1:nrow(keptIndividualsFromChildren), adjust, replace=FALSE)
      keptIndividualsFromChildren<-keptIndividualsFromChildren[-toRemove,]
    }

    keptIndividuals<-rbind(keptIndividualsFromParents,keptIndividualsFromChildren)
    rownames(keptIndividuals)<-1:nrow(keptIndividuals)

    kDGenes[1,]<-kDGenes[1,] + colSums(keptIndividualsFromParents[,-1])

    if((nrow(discardedIndividuals)>0)&&(ncol(discardedIndividuals)>0)){
      kdGenes[1,]<-kdGenes[1,] + colSums(discardedIndividuals[,-1])
    }

    if(pMutationChance!=0){
      keptIndividuals<-Mutation(keptIndividuals, pMutationChance)
    }

    for (i in 1:dim(keptIndividuals)[1]) {
      if((sum(keptIndividuals[i,2:dim(keptIndividuals)[2]]))<3){
        index<-keptIndividuals[i,]==0
        interval<-1:dim(keptIndividuals)[2]
        interval<-interval[index]
        index<-sample(interval,1, replace=FALSE)
        keptIndividuals[i,index]<-1
        rm(index)
        rm(interval)
      }
    }

    keptIndividuals<-RandomizePop(keptIndividuals)
    keptIndividuals<-Individuals(keptIndividuals[,-1])
    cat(paste("\tPopulation Size in iteration no. ", iteration," = ", nrow(keptIndividuals), "\n"))

    PlotGenAlg(kDGenes, kdGenes, MaxAcc, MeanAcc)

    flush.console()

    individuals<-keptIndividuals

    iteration <- iteration+1
    if(iteration>=iterations) break()

  }
  rezultate<-list(DGenes=kDGenes, dGenes=kdGenes, MaximumAccuracy=MaxAcc, MeanAccuracy=MeanAcc, MinAccuracy=MinAcc, BestIndividuals=bestIndividual[-1,])
  return(rezultate)
}
####################################################################
