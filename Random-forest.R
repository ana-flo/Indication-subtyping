rm(list=ls())
Sys.setenv(lang="en")
library(tidyverse)
library(data.table)
library(skimr)
library(tidymodels)
library(visdat)
library(modelStudio)
library(DALEX)
library(DALEXtra)
library(ggpubr)

###!!!!!!!!!!!!!!!!!
##########clean the  data for mutation so i can use it 

#=========================================================================================
#This 
source("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/Rscripts/Rlibraries/ClusteringFunctionLibraryRNAseq.R")
source("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/Rscripts/AnalysisFunctionLibrary.R")
#============================================
##-------------------

give.n <- function(x){
  return(data.frame(y = 0, label = paste0("n=", length(x))))
}

#------------
#read data 
#----------


df <- fread("C:/data/B_cell_Lymphoma_DataSets/remoldb/GSE117556-annotation.txt",sep="\t", header=TRUE, data.table=FALSE)
df.counts <- fread("C:/data/B_cell_Lymphoma_DataSets/remoldb/GSE117556-counts.txt",sep="\t", header=TRUE, data.table=FALSE)

#df <- read.table("C:/data/AML bulk/beatAML2/beataml_wv1to4_clinical_BM_GSVA.txt",sep="\t", header=TRUE)

rownames(df) <- make.unique(df$display_label)

df.annotation <- fread("C:/data/AML bulk/beatAML2/beataml_wv1to4_clinical_BM_GSVA.txt",sep="\t", header=TRUE, data.table = FALSE)
#df.annotation <-  fread("C:/data/AML bulk/beatAML2/beataml_wv1to4_clinical.txt",sep="\t", header=TRUE, data.table=FALSE)
df.t <- data.frame(dbgap_rnaseq_sample=colnames(df[,-c(1:5)]), t(df[,-c(1:5)]))
df.joined <- merge(df.annotation, df.t, by="dbgap_rnaseq_sample")

df.annotation <- df.annotation[,c(1:98,which(colnames(df.annotation)=="Flotetuzumab_response"))]
ListNonAML <- c("Acute erythroid leukaemia", "Acute megakaryoblastic leukaemia", "Acute promyelocytic leukaemia with t(15;17)(q22;q12); PML-RARA", "Acute undifferentiated leukaemia", 
                "Atypical chronic myeloid leukaemia, BCR-ABL1 negative", "Chronic myelomonocytic leukaemia", "Mixed phenotype acute leukaemia, B/myeloid, NOS", "Myelodysplastic syndrome associated with isolated del(5q)", 
                "Myelodysplastic/myeloproliferative neoplasm, unclassifiable", "Myeloid sarcoma", "Primary myelofibrosis", "Refractory cytopenia with multilineage dysplasia", 
                "unknown", "Acute megakaryoblastic leukaemia", "Blastic plasmacytoid dendritic cell neoplasm", "Essential thrombocythaemia", "Mixed phenotype acute leukaemia, T/myeloid, NOS", 
                "Myelodysplastic syndrome, unclassifiable", "Myeloid leukaemia associated with Down syndrome", "Plasma cell myeloma", "Refractory anaemia with excess blasts", "Unknown", "Mastocytosis", "Secondary myelofibrosis")

#first clean the data, remove irelevant columns and look at what i want to classify
skim(df.annotation)
#------------
#columns to remove : centerID,  overallSurvival, cohort, "used_manuscript_analyses" ,  "manuscript_dnaseq", "dbgap_rnaseq_sample",
#"manuscript_rnaseq", inferred_sex","inferred_ethnicity","rnaSeq","analysisRnaSeq","dbgap_dnaseq_sample"  ,"manuscript_inhibitor". "exomeSeq" ,"analysisExomeSeq" 
#columns that are relevant for survival but not for the sample
#"cumulativeChemo","cumulativeTreatmentTypeCount","cumulativeTreatmentRegimens" ,"mostRecentTreatmentType" ,"mostRecentTreatmentDuration"  , "causeOfDeath" ,"cumulativeTreatmentTypes", 
# "cumulativeTreatmentStageCount" ,"typeInductionTx" ,"currentRegimen"    , "vitalStatus",  "analysisDrug"  , "cumulativeTreatmentRegimenCount","cumulativeTreatmentStages" ,responseDurationToInductionTx"  ,
#"totalDrug"  ,"currentStage",  "overallSurvival"

# "surfaceAntigensImmunohistochemicalStains" ,"variantSummary" - require splitting to be usable for modeling 
#-------------------
df.clas <- df.annotation[which(!(df.annotation$specificDxAtInclusion %in% ListNonAML) & is.na(df.annotation$specificDxAtInclusion)==FALSE  & df.annotation$specimenType=="Bone Marrow Aspirate"),!(colnames(df.annotation) %in% c("dbgap_subject_id","centerID",  "overallSurvival", "cohort", "used_manuscript_analyses" ,  "manuscript_dnaseq", "dbgap_rnaseq_sample",
                                                           "manuscript_rnaseq", "inferred_sex","inferred_ethnicity","rnaSeq","analysisRnaSeq","dbgap_dnaseq_sample"  ,"manuscript_inhibitor", "exomeSeq" ,"analysisExomeSeq", 
                                                           "cumulativeChemo","cumulativeTreatmentTypeCount","cumulativeTreatmentRegimens" ,"mostRecentTreatmentType" ,"mostRecentTreatmentDuration"  , "causeOfDeath" ,"cumulativeTreatmentTypes", 
                                                            "cumulativeTreatmentStageCount" ,"typeInductionTx" ,"currentRegimen"    , "vitalStatus",  "analysisDrug"  , "cumulativeTreatmentRegimenCount","cumulativeTreatmentStages" ,"responseDurationToInductionTx",
                                                           "totalDrug"  ,"currentStage",  "overallSurvival", "surfaceAntigensImmunohistochemicalStains", "karyotype", "variantSummary", "otherCytogenetics", "RUNX1", "ASXL1", "TP53"))]
#----------------------------------
# test coding by blast count - first code 
df.clas$`Percent.Blasts.in.BM` <- as.numeric(df.clas$`Percent.Blasts.in.BM`)
df.clas$Blast.factor <- NA
for (i in 1: nrow(df.clas)){
  if (!is.na(df.clas$`Percent.Blasts.in.BM`[i]) & df.clas$`Percent.Blasts.in.BM`[i] < 25 ) df.clas$Blast.factor[i] <- "low"
  if (!is.na(df.clas$`Percent.Blasts.in.BM`[i]) & df.clas$`Percent.Blasts.in.BM`[i] > 25 & df.clas$`Percent.Blasts.in.BM`[i] < 50 ) df.clas$Blast.factor[i] <- "med low"
  if (!is.na(df.clas$`Percent.Blasts.in.BM`[i]) & df.clas$`Percent.Blasts.in.BM`[i] > 50 & df.clas$`Percent.Blasts.in.BM`[i] < 75) df.clas$Blast.factor[i] <- "med high"
  if (!is.na(df.clas$`Percent.Blasts.in.BM`[i]) & df.clas$`Percent.Blasts.in.BM`[i] > 75 ) df.clas$Blast.factor[i] <- "high"
}


# test coding by blast count - first code 

df.clas$Flotetuzumab_response.factor <- NA
qq <- quantile(df.clas$Flotetuzumab_response)
for (i in 1: nrow(df.clas)){
  if (!is.na(df.clas$Flotetuzumab_response[i]) & df.clas$Flotetuzumab_response[i] < qq[[2]] ) df.clas$Flotetuzumab_response.factor[i] <- "low"
  if (!is.na(df.clas$Flotetuzumab_response[i]) & df.clas$Flotetuzumab_response[i] > qq[[2]] & df.clas$Flotetuzumab_response[i] < qq[[2]] ) df.clas$Flotetuzumab_response.factor[i] <- "med low"
  if (!is.na(df.clas$Flotetuzumab_response[i]) & df.clas$Flotetuzumab_response[i] > qq[[3]] & df.clas$Flotetuzumab_response[i] < qq[[4]]) df.clas$Flotetuzumab_response.factor[i] <- "med high"
  if (!is.na(df.clas$Flotetuzumab_response[i]) & df.clas$Flotetuzumab_response[i] > qq[[4]] ) df.clas$Flotetuzumab_response.factor[i] <- "high"
}

#no2 remove samples that do not have blast.fact and remove balst count column
df.clas0 <- df.clas
df.clas$`Percent.Blasts.in.BM` <- NULL
df.clas$`Flotetuzumab_response` <- NULL
df.clas <- df.clas[which(!is.na(df.clas$Blast.factor)),]

vis_miss(df.clas, sort_miss = TRUE)

#keep age at specimen aquisiton 

#remove columns with too many missing values 
# "%.Immature.Granulocytes.in.PB","%.Nucleated.RBCs.in.PB" ,"LDH", "MCV","%.Basophils.in.PB"  ,"creatinine","totalProtein" , 
#"%.Neutrophils.in.PB" , "AST","ALT","%.Blasts.in.PB","%.Lymphocytes.in.PB" , "%.Monocytes.in.PB","hematocrit", "hemoglobin","%.Eosinophils.in.PB" ,"albumin" , "plateletCount"

#remove columns with more  than 10% missing values 
df.clas <- df.clas[which(is.na(df.clas$ageAtSpecimenAcquisition)==FALSE),which(!(colnames(df.clas) %in% c("ageAtDiagnosis","ageAtSpecimenAcquisition","wbcCount", "allelic_ratio", "Percent.Immature.Granulocytes.in.PB","Percent.Nucleated.RBCs.in.PB" ,"LDH", "MCV","Percent.Basophils.in.PB"  ,"creatinine","totalProtein" , 
"Percent.Neutrophils.in.PB" , "AST","albumin","ALT","Percent.Blasts.in.PB","albumin","Percent.Lymphocytes.in.PB" , "Percent.Monocytes.in.PB","hematocrit", "hemoglobin","Percent.Eosinophils.in.PB"  , "plateletCount",  "ELN2017")))]
vis_miss(df.clas, sort_miss = TRUE)

# split dataset into training and test

df.clas[sapply(df.clas, is.character)] <- lapply(df.clas[sapply(df.clas, is.character)], 
                                       as.factor)
str(df.clas)



# Put 3/4 of the data into the training set 
data_split <- initial_split(df.clas, 
                            prop = 3/4, 
                            strata = Flotetuzumab_response.factor)



# Create dataframes for the two sets:
train_data <- training(data_split) 

AML_recipe <- training(data_split) %>%
  recipe(Flotetuzumab_response.factor ~.) %>%
  step_corr(all_predictors()) %>%
  step_center(all_predictors(), -all_outcomes()) %>%
  step_scale(all_predictors(), -all_outcomes()) %>%
  prep()

test_data <- testing(data_split)

#run the model
clas.test <-  rand_forest(trees = 100, mode = "classification") %>%
  set_engine("ranger") %>% set_mode("classification") %>%
  fit(Flotetuzumab_response.factor ~ ., data = train_data)

#add the predicted values to the initial data frame as first column 

ml.output <-  predict(clas.test, test_data) %>%
  bind_cols(test_data) %>%
  glimpse()

ml.output %>% metrics(trut=Blast.factor, estimate=.pred_class)
#label is the name of the model - should macth the packahe used or let empty

explainer <- DALEX::explain(model = clas.test, data = train_data[,-ncol(train_data)], y=train_data$Flotetuzumab_response.factor)

vip_rf <- model_parts(explainer)
vip_rms <- model_parts(explainer, loss_function = loss_root_mean_square)
#-------------
#plot
ggplot_imp <- function(...) {
  obj <- list(...)
  metric_name <- attr(obj[[1]], "loss_name")
  metric_lab <- paste(metric_name, 
                      "after permutations\n(higher indicates more important)")
  
  full_vip <- bind_rows(obj) %>%
    filter(variable != "_baseline_")
  
  perm_vals <- full_vip %>% 
    filter(variable == "_full_model_") %>% 
    group_by(label) %>% 
    summarise(dropout_loss = mean(dropout_loss))
  
  p <- full_vip %>%
    filter(variable != "_full_model_") %>% 
    mutate(variable = fct_reorder(variable, dropout_loss)) %>%
    ggplot(aes(dropout_loss, variable)) 
  if(length(obj) > 1) {
    p <- p + 
      facet_wrap(vars(label)) +
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss, color = label),
                 size = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(aes(color = label, fill = label), alpha = 0.2)
  } else {
    p <- p + 
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss),
                 size = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(fill = "#91CBD765", alpha = 0.4)
    
  }
  p +
    theme(legend.position = "none") +
    labs(x = metric_lab, 
         y = NULL,  fill = NULL,  color = NULL)
}

ggplot_imp(vip_rf)+theme_light()


ggplot(df.clas0, aes(x=`specificDxAtAcquisition`, y=`Flotetuzumab_response`, fill=`specificDxAtAcquisition`))+geom_boxplot()+geom_point()+stat_compare_means()+theme_light()
