#Create Diagnosis Groups using ICD 9 and 10 codes

# Libraries ----
library(here)
library(tidyverse)

#Data
## Load ED Visit Data 
## Example dataset used below is from Nevada called nv_dxcleaner
## Functions work on a column of all different ICD codes from a given visit separated by commas


#Create Diagnosis Group Variables
dx_groups <- list(
  nasopharyngitis   = "J00|460",
  sinusitis         = "J01|461", 
  pharyngitis       = "J02|462",
  tonsilitis        = "J03|463",
  laryngitis        = "J04|J05|464",
  unsp_uri          = "J06|465",
  influenza         = "J09|J10|J11|487\\.0|487\\.1|487\\.2|487\\.3|487\\.4|487\\.5|487\\.6|487\\.7|487\\.8|487\\.9|488",
  pneumonia         = "J12|J13|J14|J15|J16|J17|J18|480|481|482|483|484|485|486",
  bronchitis        = "J20|J40|J41|J42|466|490|491",
  bronchiolitis     = "J21",
  lower_resp_inf    = "J22|519\\.8",
  allergy_rhinitis  = "J30|J31|477",
  emphysema         = "J43|492",
  copd              = "J44|496",
  asthma_attack     = "J45|493",
  bronchiectasis    = "J47|494",
  pneumoconiosis    = "J60|J61|J62|J63|J64|J65|500|501|502|503|504|505|506|507|508|515",
  inh_ext_agents    = "J66|J67|J68|J70|495",
  pneumonitis       = "J69|507",
  ards              = "J80",
  pul_edema         = "J81|J82|J83|J84|518\\.4",
  pleural_effusion  = "J90|J91|J92|J93|J94|511",
  resp_failure      = "J96|518\\.8",
  unsp_resp_disorder= "J98|J99|519\\.9")


# Build aggregate patterns from the atomic groups -- no copy-pasting
uad_groups <- c("unsp_uri", "laryngitis", "allergy_rhinitis",
                "nasopharyngitis", "sinusitis", "pharyngitis", "tonsilitis")

lad_groups <- c("asthma_attack", "copd", "pneumonitis", "lower_resp_inf",
                "bronchiolitis", "bronchitis", "emphysema", "pneumoconiosis",
                "inh_ext_agents", "pneumonia", "influenza", "bronchiectasis",
                "ards", "resp_failure", "pul_edema", "pleural_effusion",
                "unsp_resp_disorder")

cvd_groups <- c("acute_MI", "angina", "heart_failure", "cardiac_arrest",
                "cerebral_infarction", "ischemic_heart_disease", "dysrhythmia")

dx_groups$uad <- paste(dx_groups[uad_groups], collapse = "|")
dx_groups$lad <- paste(dx_groups[lad_groups], collapse = "|")
dx_groups$cvd <- paste(dx_groups[cvd_groups], collapse = "|")


#create matching function
match_diag <- function(pattern, diag_code) {
  # Wraps your combined pattern strings so that they must match from the start (^)
  # e.g., converts "J45|493" to "^(J45|493)"
  anchored_pattern <- paste0("^(", pattern, ")")
  
  ifelse(is.na(diag_code), NA_integer_, as.integer(grepl(anchored_pattern, diag_code)))
}



###apply matching function
# 

for (nm in names(dx_groups)) {
  nv_dxcleaner[[paste0(nm, "_prin")]] <- match_diag(dx_groups[[nm]], nv_dxcleaner$dx_prin_clean)
  nv_dxcleaner[[paste0(nm, "_sec")]]  <- match_diag(dx_groups[[nm]], nv_dxcleaner$dx_sec_clean)
  nv_dxcleaner[[paste0(nm, "_any")]]  <- match_diag(dx_groups[[nm]], nv_dxcleaner$dx_clean)
}



# Function For Matching all Respiratory ICD Codes 
match_respiratory <- function(diag_code) {
  pattern <- "J[0-9]{2}|4[6-9][0-9]|50[0-9]|51[0-9]"
  ifelse(is.na(diag_code), NA_integer_, as.integer(grepl(pattern, diag_code)))
}



# make grouped category for all respiratory outcomes
nv_dx_update <- nv_dxcleaner %>%
  mutate(
    resp_prin = match_respiratory(dx_prin_clean),
    resp_sec  = match_respiratory(dx_sec_clean),
    resp_any  = match_respiratory(dx_clean))


#make daily counts
nv_winter_daily<-nv_dx_update %>% 
  group_by(Date,city) %>%
  summarise(across(contains(c("prin","sec","any")), \(x) sum(x, na.rm = TRUE))) %>% ungroup()










