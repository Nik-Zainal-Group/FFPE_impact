library(signature.tools.lib)

#Load catalogues of interest

FFPE_impact <- function(SNV_path, ID_path, organ){

set.seed(1)
  
FFPE_SNV <- read.table(SNV_path)

FFPE_ID <- read.table(ID_path)

#Define SBS signature catalogues
SBS_sigs <- getSignaturesForFitting(organ)
SBS_sigs <- SBS_sigs$common

SBS_artefact <- readRDS("SBS_artefact.RDS")

SBS_sigs <- cbind(SBS_sigs, SBS_artefact)

#Define ID signature catalogues
ID_sigs <- readRDS("ID_sigs.RDS")
rownames(ID_sigs) <- ID_sigs$MutationType
ID_sigs$MutationType <- NULL
colnames(ID_sigs)[colnames(ID_sigs) == "S2"] <- "ID_FFPE"

#Fit SBS signatures

FFPE_Fit <- Fit(catalogues = FFPE_SNV,
                signatures = SBS_sigs,
                useBootstrap = T,
                nboot = 50,
                threshold_percent = 5,
                threshold_p.value = 0.05,
                nparallel=20)
FFPE_exposures <- FFPE_Fit$exposures
total_snv <- rowSums(FFPE_exposures)
FFPE_exposures <- cbind(FFPE_exposures, total_snv)


#Fit ID signatures

ID_FFPE_Fit <- Fit(catalogues = FFPE_ID,
                   signatures = ID_sigs,
                   useBootstrap = T,
                   nboot = 50,
                   threshold_percent = 5,
                   threshold_p.value = 0.05,
                   nparallel=20)
ID_FFPE_Fit_exposures <- ID_FFPE_Fit$exposures
total_ID <- rowSums(ID_FFPE_Fit_exposures)
ID_FFPE_Fit_exposures <- cbind(ID_FFPE_Fit_exposures, total_ID)

#Calculate FFPEimpact

Combined_mutations <- cbind(FFPE_exposures, ID_FFPE_Fit_exposures)
Combined_mutations <- Combined_mutations[,which(colnames(Combined_mutations) %in% c("SBSnew", "SBS57", "total_snv", "ID_FFPE", "total_ID"))]
Combined_mutations <- as.data.frame(Combined_mutations)
Combined_mutations$total_SNV_artefact <- Combined_mutations$SBS57 + Combined_mutations$SBSnew
Combined_mutations$FFPEimpact <- ((Combined_mutations$total_SNV_artefact/Combined_mutations$total_snv)+(Combined_mutations$ID_FFPE/Combined_mutations$total_ID))/2

#Generating new del.mh.prop by deleting artefact indels

catalogue <- FFPE_ID
exposures <- as.data.frame(ID_FFPE_Fit_exposures)
IDFFPE <- readRDS("IDFFPE.RDS")

artefactonly <- sapply(exposures$ID_FFPE, function(x){
  IDFFPE$S2*x
})

colnames(artefactonly) <- rownames(exposures)
rownames(artefactonly) <- rownames(catalogue)

cleancatalogue <- catalogue - artefactonly
cleancatalogue[cleancatalogue<0] <- 0 

all_deletions <- grepl(x = rownames(cleancatalogue), pattern = "Del", fixed = T)

a <- cleancatalogue[all_deletions,]
MHcount <- rbind(a,colSums(a))
MHcount <- rbind(MHcount, colSums(MHcount[36:47,]))
MHcount <- MHcount[-c(1:47),]
MHcount <- as.data.frame(t(MHcount))
MHcount$total_deletions <- MHcount[,1] 
MHcount$MH_number <- MHcount[,2]
MHcount <- MHcount[,-c(1,2)]
MHcount$del.mh.prop <- MHcount$MH_number/MHcount$total_deletions
MHcount$sampleID <- row.names(MHcount)

return(list("SBS_exposures"=FFPE_exposures,
            "ID_exposures"=ID_FFPE_Fit_exposures,
            "FFPE_impact"=Combined_mutations,
            "cleaned_catalogues"=cleancatalogue,
            "newMHcount"=MHcount
))

}

