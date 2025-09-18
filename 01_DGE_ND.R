# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Manuscript: Unveiling Common Transcriptomic Features between Melanoma Brain 
# Metastases and Neurodegenerative Diseases
# Authors: José Francisco Català Senent, Adolfo López Cerdanm Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# --- Load libraries ----------------------------------------------------------
library(Biobase)
library(limma)
library(edgeR)
library(dplyr)

# Local functions
limma.local = function(data, groups, cont, batch = NULL) {
  
  # Differential expression
  groups = as.factor(groups)
  
  if (!is.null(batch)) {
    design = model.matrix(~ 0 + groups + batch)
    colnames(design) = c(levels(groups), "batch")
  } else {
    design = model.matrix(~ 0 + groups)
    colnames(design) = levels(groups)
  }
  
  fit = lmFit(data, design)
  contrast.matrix = makeContrasts(contrasts = cont, levels = design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  
  SE.coef = sqrt(fit2$s2.post) * fit2$stdev.unscaled
  table(rownames(SE.coef) == rownames(fit2$coefficients))
  mat = cbind(fit2$coefficients, SE.coef)
  
  # Results
  results = lapply(1:length(cont), function(i) {
    tt = topTable(fit2, coef = i, adjust.method = "BH", number = "Inf", sort.by = sort, confint = sdev)
    tt$ID = rownames(tt)
    if (sdev == TRUE) {
      contrast = cont[[i]]
      submat = mat[, colnames(mat) == contrast]
      colnames(submat) = c("coef", "se.coef")
      
      tt[, "SE"] = (tt[, "CI.R"] - tt[, "CI.L"]) / 3.92
      
      head(submat)
      int = intersect(rownames(tt), rownames(submat))
      tt = cbind(tt, submat)
      head(tt)
      tt = tt[, c(9, 1:7, 10:12)]
    } else {
      tt = tt[, c(6, 1:5)]
    }
    tt
  })
  names(results) = names(cont)
  
  return(results)
}


limma.voom.local = function(data, groups, cont, batch = NULL) {
  
  ## Differential expression
  groups = as.factor(groups)
  if (!is.null(batch)) {
    design = model.matrix(~ 0 + groups + batch)
    colnames(design) = c(levels(groups), "batch")
  } else {
    design = model.matrix(~ 0 + groups)
    colnames(design) = levels(groups)
  }
  
  # Voom
  dgeL = DGEList(counts = data, group = groups)
  dgeN = calcNormFactors(dgeL, method = "TMM")
  v = voom(dgeN, design)
  
  fit = lmFit(v, design)
  contrast.matrix = makeContrasts(contrasts = cont, levels = design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2)
  
  SE.coef = sqrt(fit2$s2.post) * fit2$stdev.unscaled
  table(rownames(SE.coef) == rownames(fit2$coefficients))
  mat = cbind(fit2$coefficients, SE.coef)
  
  # Results
  results = lapply(1:length(cont), function(i) {
    tt = topTable(fit2, coef = i, adjust.method = "BH", number = "Inf", sort.by = sort, confint = sdev)
    tt$ID = rownames(tt)
    if (sdev == TRUE) {
      contrast = cont[[i]]
      submat = mat[, colnames(mat) == contrast]
      colnames(submat) = c("coef", "se.coef")
      
      tt[, "SE"] = (tt[, "CI.R"] - tt[, "CI.L"]) / 3.92
      
      head(submat)
      int = intersect(rownames(tt), rownames(submat))
      tt = cbind(tt, submat)
      head(tt)
      tt = tt[, c(9, 1:7, 10:12)]
    } else {
      tt = tt[, c(7, 1:5)]
    }
    tt
  })
  names(results) = names(cont)
  
  return(results)
}


################################################################################
#   DGE microarrays                                                            #
################################################################################
dataset = readRDS("TBD")
results = limma.local(data = dataset, groups = c("control", "ND"), cont = "disease", batch = "TBD")
write.table(x = results, file = "TBD")


################################################################################
#   DGE RNA-seq                                                                #
################################################################################
dataset = readRDS("TBD")
results = limma.voom.local(data = dataset, groups = c("control", "ND"), cont = "disease", batch = "TBD")
write.table(x = results, file = "TBD")



