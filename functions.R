gnomeQuest <- function(x, type){
  #x = 5 column dataframe. $1=chrom, $2=start position, $3=stop position,
  #4 = variant, $5 = mutation type (SNV or INDEL)
  #exports a list of queries to submit to gnomad
  querymaker <- function(start, stop, chrom, ref_genome, dataset_id){
    paste0('{\n  region(start: ', start, ', stop: ', stop,
           ', chrom: "', chrom, '", reference_genome: ', ref_genome, ') {\n    clinvar_variants {\n      clinical_significance\n      clinvar_variation_id\n      gnomad {\n        exome {\n          ac\n          an\n          filters\n        }\n        genome {\n          ac\n          an\n          filters\n        }\n      }\n      gold_stars\n      hgvsc\n      hgvsp\n      in_gnomad\n      major_consequence\n      pos\n      review_status\n      transcript_id\n      variant_id\n    }\n    variants(dataset: ', dataset_id, ') {\n      consequence\n      flags\n      gene_id\n      gene_symbol\n      hgvs\n      hgvsc\n      hgvsp\n      lof\n      lof_filter\n      lof_flags\n      pos\n      rsids\n      transcript_id\n      transcript_version\n      variant_id\n      exome {\n        ac\n        ac_hemi\n        ac_hom\n        an\n        af\n        filters\n        populations {\n          id\n          ac\n          an\n          ac_hemi\n          ac_hom\n        }\n      }\n      genome {\n        ac\n        ac_hemi\n        ac_hom\n        an\n        af\n        filters\n        populations {\n          id\n          ac\n          an\n          ac_hemi\n          ac_hom\n        }\n      }\n      lof_curation {\n        verdict\n        flags\n      }\n    }\n  }\n}')
  }
  
  querymakerList <- list()
  
  if (type == "SNV") {
    for (i in seq_len(nrow(x))) {
      querymakerList[[i]] <- querymaker(start = x[i,2], 
                                        stop = x[i,2], 
                                        chrom = x[i,1], 
                                        ref_genome = "GRCh38", 
                                        dataset_id = "gnomad_r3")
    }
    names(querymakerList) <- paste0(x[,1], '_', x[,2])
  } else if(type == "INDEL"){
    for (i in seq_len(nrow(x))) {
      querymakerList[[i]] <- querymaker(start = x[i,2], 
                                      stop = x[i,3], 
                                      chrom = x[i,1], 
                                      ref_genome = "GRCh38", 
                                      dataset_id = "gnomad_r3")
    }
    names(querymakerList) <- paste0(x[,1], '_', x[,2], '-', x[,3])
  }
  return(querymakerList)
}

gnomeQuest2 <- function(x, type){
  #x = 5 column dataframe. $1=chrom, $2=start position, $3=stop position,
  #4 = variant, $5 = mutation type (SNV or INDEL)
  #exports a list of queries to submit to gnomad
  querymaker <- function(start, stop, chrom, ref_genome, dataset_id){
    paste0('{\n  region(start: ', start, ', stop: ', stop,
           ', chrom: "', chrom, '", reference_genome: ', ref_genome, ') {\n    clinvar_variants {\n      clinical_significance\n      clinvar_variation_id\n      gnomad {\n        exome {\n          ac\n          an\n          filters\n        }\n        genome {\n          ac\n          an\n          filters\n        }\n      }\n      gold_stars\n      hgvsc\n      hgvsp\n      in_gnomad\n      major_consequence\n      pos\n      review_status\n      transcript_id\n      variant_id\n    }\n    variants(dataset: ', dataset_id, ') {\n      consequence\n      flags\n      gene_id\n      gene_symbol\n      hgvs\n      hgvsc\n      hgvsp\n      lof\n      lof_filter\n      lof_flags\n      pos\n      rsids\n      transcript_id\n      transcript_version\n      variant_id\n      exome {\n        ac\n        ac_hemi\n        ac_hom\n        an\n        af\n        filters\n        populations {\n          id\n          ac\n          an\n          ac_hemi\n          ac_hom\n        }\n      }\n      genome {\n        ac\n        ac_hemi\n        ac_hom\n        an\n        af\n        filters\n        populations {\n          id\n          ac\n          an\n          ac_hemi\n          ac_hom\n        }\n      }\n      lof_curation {\n        verdict\n        flags\n      }\n    }\n  }\n}')
  }
  
  querymakerList <- list()
  
  if (type == "SNV") {
    for (i in seq_len(nrow(x))) {
      querymakerList[[i]] <- querymaker(start = x[i,2], 
                                        stop = x[i,2], 
                                        chrom = x[i,1], 
                                        ref_genome = "GRCh37", 
                                        dataset_id = "gnomad_r2_1")
    }
    names(querymakerList) <- paste0(x[,1], '_', x[,2])
  } else if(type == "INDEL"){
    for (i in seq_len(nrow(x))) {
      querymakerList[[i]] <- querymaker(start = x[i,2], 
                                        stop = x[i,3], 
                                        chrom = x[i,1], 
                                        ref_genome = "GRCh37", 
                                        dataset_id = "gnomad_r2_1")
    }
    names(querymakerList) <- paste0(x[,1], '_', x[,2], '-', x[,3])
  }
  return(querymakerList)
}

library(dplyr)
gnomeCatch <- function(x, gnomeQuestList, batch=5, sleep=20, type){
  jsonlist <- list()
  for (i in seq_len(nrow(x))) {
    if(i %% batch == 0) {Sys.sleep(sleep)}
    tryCatch(
      expr = {jsonlist[[i]] <- request("https://gnomad.broadinstitute.org/api/?") %>%
        req_body_json(list(query=gnomeQuestList[[i]], variables="null")) %>%
        req_perform() %>%
        resp_body_json()
      if (type == "SNV") {
        names(jsonlist[i]) <- paste0(x[,1][i], '_', x[,2][i])
      }else if (type == "INDEL") {
        names(jsonlist[i]) <- paste0(x[,1][i], '_', x[,2][i], '-', x[,3][i])
      }},
      error = function(e){
        message(paste0('Variant ', i, ' Caught an Error! - may need to re-run'))
        print(e)},
      warning = function(w){
        message(paste0('Variant ', i, ' Caught a Warning!'))
        print(w)},
      finally = {
        message(paste0('Variant ', i, ' Done!'))
      })}
  return(jsonlist)
}

gnomeCatch2 <- function(x, gnomeQuestList, batch=5, sleep=20, type){
  jsonlist <- list()
  for (i in seq_len(nrow(x))) {
    if(i %% batch == 0) {Sys.sleep(sleep)}
    tryCatch(
      expr = {jsonlist[[i]] <- request("https://gnomad.broadinstitute.org/api/?") %>%
        req_body_json(list(query=gnomeQuestList[[i]], variables="null")) %>%
        req_perform() %>%
        resp_body_json()
      if (type == "SNV") {
        names(jsonlist[[i]]) <- x$Name[i]
      }else if (type == "INDEL") {
        names(jsonlist[[i]]) <- x$Name[i]
      }
      message(paste0("Successfully caught gnomes at ", x$Name[i]))
      },
      error = function(e){
        message(paste0('Variant ', x$Name[i], ' was not caught - may need to re-run'))
        print(e)},
      warning = function(w){
        message(paste0('Variant ', i, ' Caught a Warning!'))
        print(w)}
    )}
  newjsonlist <- unlist(jsonlist,recursive=FALSE)
  return(newjsonlist)
}

gnomeMisses <- function(x, gnomeCatchList, type){
  dat <-  which(!lengths(gnomeCatchList))
  queryMissedList <- queryMakerList[dat]
  dfMissed <- df[dat,]
  
  catchList_2 <- gnomeCatch(x = dfMissed, gnomeQuestList = queryMissedList,
                                batch = 5, sleep = 20, type = type)
  
  for (i in seq_along(dat)) {
    this <- dat[[i]]
    gnomeCatchList[[this]] <- catchList_2[[i]]
  }
  return(gnomeCatchList)
}

gnomeVars <- function(x, querySentList, type){
  #Entry is the Jsonlist from gnomeCatch
  #exports variants which have data in gnomad, web-site info which have
  #data in gnomad, names of variants with data, number of variants found
  #in gnomad
  
  variantsList <- lapply(querySentList, function(X) {
    sapply(X$region$variants, function(X)X$variant_id)
  })
  
  if (type == "SNV") {
    names(variantsList) <- paste0(x[,1], '_', x[,2])
  }else if (type == "INDEL") {
    names(variantsList) <- paste0(x[,1], '_', x[,2], '-', x[,3])
  }
  
  return(variantsList)
}

gnomeAlleles <- function(variantList, siteList){
  
  genomicList <- vector("list", number+1)
  for (i in (seq_along(variantList))) {
    for(j in 1:length(siteList[[i]]$region$variants)){
      print(paste("This is iteration i =", i, "and j =", j))
      genomicList[[i]][[j]] <- siteList[[i]]$region$variants[[j]]$genome
      names(genomicList)[i] <- names(variantList)[i]
      #names(genomicList[[i]][j]) <- variantList[[i]][j]
    }
  }
  
  
  genomicList2 <- vector("list", number+1)
  for (i in (seq_along(variantList))) {
    for(j in 1:length(siteList[[i]]$region$variants)){
      print(paste("This is iteration i =", i, "and j =", j))
      genomicList2[[i]][[j]] <- genomicList[[i]][[j]][-c(6,7)]
      names(genomicList2)[i] <- names(variantList)[i]
      #names(genomicList2[[i]][j]) <- variantList[[i]][j]
    }
  }
  
  return(genomicList2)
}

gnomeDataFrame <- function(listOfAlleles, variants, 
                           siteData, varNames, allVariants){
  dfs <- unlist(listOfAlleles,recursive=FALSE)
  
  for (i in seq_len(length(dfs))) {
    if (length(dfs[[i]]) == 0) {
      dfs[[i]] <- list("NA", "NA", "NA", "NA", "NA")
    }
  }
  
  dfs <- do.call(rbind.data.frame, dfs)
  colnames(dfs) <- c("AlleleCount", "AlleleCountHemizygotes", "AlleleCountHomozygotes", 
                     "AlleleNumber", "AlleleFrequency")
  
  reparray <- vector()
  for (i in (seq_along(variants))) {
    reparray[i] <- length(listOfAlleles[[i]])
  }
  
  listrep <- list()
  for (i in (seq_along(variants))) {
    listrep[[i]] <- rep(varNames[i], reparray[i])
  }
  dfs$position <- unlist(listrep)
  
  dforder <- dfs[mixedorder(dfs$position),]
  dforder <- dforder[c(6, 1:5)]
  
  vec <- vector()
  for (i in (seq_along(variants))) {
    for(j in 1:length(siteList[[i]]$region$variants)){
      if (isFALSE(is.null(siteList[[i]]$region$variants[[j]]$genome))) {
        #print(paste("This is iteration i =", i, "and j =", j))
        element <- variants[[i]][[j]]
        vec <- c(vec, element)
      }
    }    
  }
  
  
  rownames(dforder) <- vec
  
  return(dforder)
}

gnomeRatify <- function(catchList, type){
  dat <-  which(!lengths(catchList))
  message("Currently ", print(length(dat)), "  samples are missing")
  if (length(dat) > 0) {
    querySentList <- gnomeMisses(x = df, gnomeCatchList = catchList,
                                 type = type)
    dat <-  which(!lengths(querySentList))
  } else (message("No need to re-run gnomeCatch - no failed variants found :)."))
  
  return(querySentList)
}

gnomeName <- function(x = df, type = type){
  
  if (type == "SNV") {
    x$Name <- paste0(paste0(x[,1], '_', x[,2]))
  }else if (type == "INDEL") {
    x$Name <- paste0(paste0(x[,1], '_', x[,2], '-', x[,3]))
  }
  return(x)
}


gnomeAlleles2 <- function(variantList, siteList){
  genomicList <- vector("list", number)
  genomicList2 <- vector("list", number)
  for (i in (seq_along(variants))) {
    for(j in 1:length(siteList[[i]]$region$variants)){
      if (isFALSE(is.null(siteList[[i]]$region$variants[[j]]$genome))) {
        print(paste("This is iteration i =", i, "and j =", j))
        genomicList[[i]][[j]] <- siteList[[i]]$region$variants[[j]]$genome
        names(genomicList)[i] <- names(variants)[i]
        genomicList2[[i]][[j]] <- genomicList[[i]][[j]][-c(6,7)]
        names(genomicList2)[i] <- names(variants)[i]
        
        
      }
    }
  }
  return(genomicList2)
}