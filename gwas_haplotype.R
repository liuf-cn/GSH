# Clear environment variables
rm(list=ls(all=TRUE))

# Load dependent packages
library(data.table)
library(BGLR)
library(MASS)
library(parallel)

# -------------------------- Configure Parameters (User-Modifiable) --------------------------
trait <- "PH"               # Trait name (consistent with the column name in the phenotype file)
chr <- "chr1"               # Chromosome number (e.g., chr1, chr2)
N_CPU_CORE <- 4             # Number of parallel cores
var.geno <- 771.992         # Genetic variance (can be estimated by BGLR)
var.err <- 43.76611         # Residual variance (can be estimated by BGLR)
windowsize <- 100000        # Inter-gene window size (bp)
path <- getwd()             # Working directory (current directory by default)
# ---------------------------------------------------------------------------

# Set output file path
outfile <- sprintf('GWAS/%s_sgene', trait)
if(!file.exists(outfile)){
  dir.create(outfile, recursive = TRUE)
}

# -------------------------- Read Input Files --------------------------
# Phenotype data
pheno <- read.table('phenotype.txt', header = TRUE, stringsAsFactors = FALSE)
# Genotype data (haplotype)
geno <- fread(sprintf('ld90/%s.hap.raw', chr), data.table = FALSE)
# Kinship matrix
K1 <- as.matrix(fread('kin_gemma.cXX.txt', data.table = FALSE))
rownames(K1) <- colnames(K1) <- geno$FID
# Genetic map (SNP position information)
map <- fread('hybrid_allchr_ld90.bim', data.table = FALSE)
colnames(map) <- c('chr', 'snp', 'mpos', 'pos', 'ref', 'alt')
map <- map[which(map$chr == gsub('chr', '', chr)), ]
rownames(map) <- map$snp
# Gene annotation file
gene <- fread('gene.gff3', data.table = FALSE)
gene <- gene[which(gene[, 1] == gsub('chr', '', chr)), ]
gene[1, 4] <- 0  # Set the start position of the first gene to 0
gene[nrow(gene), 5] <- max(map$pos) + 2  # Set the end position of the last gene to maximum SNP position + 2
# Haplotype allele files
allele3 <- read.table('hap_alleles.txt')
allele2 <- read.table('hap2_alleles.txt')

# -------------------------- Data Preprocessing --------------------------
# Match phenotype and genotype individuals (remove individuals with missing phenotypes)
geno$PHENOTYPE <- pheno[, grep(trait, colnames(pheno))]
ind <- which(is.na(geno$PHENOTYPE))
if(length(ind) > 0){
  geno <- geno[-ind, ]
  K1 <- K1[-ind, -ind]
}
ngen <- nrow(geno)

# Convert mixed linear model to simple linear model (variance component adjustment)
I <- diag(1, nrow(K1), ncol(K1))
decomposition <- eigen(K1)
D <- decomposition$values
TU <- t(decomposition$vectors)
W <- diag(1/sqrt(D*(var.geno/var.err) + 1))
coef <- W %*% TU

# Construct analysis dataset
onlyDATA <- data.frame(
  Y = coef %*% geno$PHENOTYPE,
  intercept = coef %*% rep(1, ngen)
)

# Initialize result file
start <- data.frame('Number', 'hap', 'chr', 'pos', 'pvalue_hap', 'pvalue_allele', 'r2')
write.table(
  start,
  file = sprintf('GWAS/%s_hap_step.txt', trait),
  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t'
)

# -------------------------- Haplotype Association Analysis Function --------------------------
HAP_LM <- function(m){
  # Define the SNP window of the current gene (based on gene annotation position)
  if(m == 1){
    start_pos <- gene[m, 4]
  }else{
    start_pos <- (gene[m-1, 5] + gene[m, 4])/2
  }
  if(m == nrow(gene)){
    end_pos <- gene[m, 5]
  }else{
    end_pos <- (gene[m, 5] + gene[m+1, 4])/2
  }
  chr_tmp <- map$chr[1]
  ind_p <- which(map$pos < end_pos & map$pos >= start_pos & map$chr == chr_tmp)
  allsnp <- map$snp[ind_p]
  
  # Perform analysis only if there are at least 2 SNPs
  if(length(allsnp) >= 2){
    tmpsnp <- allsnp
    nsplit <- ceiling(length(tmpsnp)/50)  # Split into SNP groups (each group ≤50 SNPs)
    for(sp in 1:nsplit){
      # Sample 50 SNPs (take all remaining for the last group)
      if(sp < nsplit){
        snp <- sample(tmpsnp, 50)
        tmpsnp <- setdiff(tmpsnp, snp)
      }else{
        snp <- tmpsnp
      }
      nsnp <- length(snp)
      
      # Stepwise regression to screen for optimal SNP combinations
      for(i in 1:nsnp){
        if(length(snp) > 1){
          snpg <- as.matrix(geno[, snp])
          tmpdata <- data.frame(Y = onlyDATA$Y, coef %*% snpg)
          lm1 <- lm(Y ~ . , data = tmpdata)
          lm2 <- stepAIC(lm1, direction = "both", trace = FALSE)
          snp1 <- rownames(summary(lm2)$coefficients)[-1]  # Screened SNPs
          
          # 3-SNP combination haplotype analysis
          if(length(snp1) >= 3){
            snp1 <- snp1[1:3]  # Take the top 3 optimal SNPs
            snp <- setdiff(snp, snp1)
            hpos <- round(mean(map[snp1, 'pos']))  # Average haplotype position
            snp_r <- paste(snp1, collapse = '_')  # SNP combination name
            ppg <- geno[, snp1]
            window <- apply(ppg, 1, paste, collapse = '')  # Construct haplotype
            hap_all <- as.matrix(allele3[window, ])
            rownames(hap_all) <- geno$FID
            hapDATA <- data.frame(onlyDATA, coef %*% hap_all)
            
            # Association test (ANOVA)
            hlm1 <- lm(Y ~ -1 + intercept, data = onlyDATA)
            hlm2 <- lm(Y ~ -1 + ., data = hapDATA)
            hap_anova <- anova(hlm1, hlm2)
            Pvalue_hap <- hap_anova$Pr[2]
            allele_min <- min(summary(hlm2)$coefficients[-1, 4])
            hadj_r <- summary(hlm2)$adj.r.squared
            
            # Save results
            result <- data.frame(m, snp_r, chr, hpos, Pvalue_hap, allele_min, hadj_r)
            write.table(
              result,
              file = sprintf('%s/%s_hap_%s', outfile, chr, m),
              quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE, sep = '\t'
            )
          }
          
          # 2-SNP combination haplotype analysis
          if(length(snp1) == 2){
            snp <- setdiff(snp, snp1)
            hpos <- round(mean(map[snp1, 'pos']))
            snp_r <- paste(snp1, collapse = '_')
            ppg <- geno[, snp1]
            window <- apply(ppg, 1, paste, collapse = '')
            hap_all <- as.matrix(allele2[window, ])
            rownames(hap_all) <- geno$FID
            hapDATA <- data.frame(onlyDATA, coef %*% hap_all)
            
            # Association test (ANOVA)
            hlm1 <- lm(Y ~ -1 + intercept, data = onlyDATA)
            hlm2 <- lm(Y ~ -1 + ., data = hapDATA)
            hap_anova <- anova(hlm1, hlm2)
            Pvalue_hap <- hap_anova$Pr[2]
            allele_min <- min(summary(hlm2)$coefficients[-1, 4])
            hadj_r <- summary(hlm2)$adj.r.squared
            
            # Save results
            result <- data.frame(m, snp_r, chr, hpos, Pvalue_hap, allele_min, hadj_r)
            write.table(
              result,
              file = sprintf('%s/%s_hap_%s', outfile, chr, m),
              quote = FALSE, append = TRUE, row.names = FALSE, col.names = FALSE, sep = '\t'
            )
          }
        }
      }
    }
  }
}

# -------------------------- Run Association Analysis in Parallel --------------------------
marker_pos <- 1:nrow(gene)
if(N_CPU_CORE > 1){
  res <- mclapply(marker_pos, mc.cores = N_CPU_CORE, FUN = HAP_LM, mc.preschedule = TRUE)
}else{
  res <- lapply(marker_pos, FUN = HAP_LM)
}

# -------------------------- Merge Results --------------------------
setwd(outfile)
system(sprintf("ls | grep '^%s_hap' | xargs cat > %s_%s_hap.txt", chr, trait, chr))
system(sprintf("cat %s_%s_hap.txt | sort -k1,1n >> ../%s_hap_step.txt", trait, chr, trait))

cat("Analysis completed! Result file: GWAS/", trait, "_hap_step.txt\n", sep = "")
