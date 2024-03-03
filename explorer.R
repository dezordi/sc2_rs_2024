library(dplyr)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(tidyr)


# This script performs a bunch of functions to get
# information of Brazilian samples.
# A lot of variables are hard-coded based on the
# necessary information to conduct the analysis
# of this paper: (in prep)

#The most recent genome in this study was collected on 2022-07-07
x_order_month <- c( "2020-01","2020-02", "2020-03", "2020-04",
                    "2020-05", "2020-06", "2020-07","2020-08",
                    "2020-09", "2020-10", "2020-11","2020-12",
                    "2021-01","2021-02", "2021-03", "2021-04", 
                    "2021-05", "2021-06", "2021-07","2021-08",
                    "2021-09", "2021-10", "2021-11", "2021-12",
                    "2022-01","2022-02", "2022-03", "2022-04",
                    "2022-05", "2022-06", "2022-07")

#The study focus on Gamma, Delta and Omicron variants, not including Omicron recombinants
target_lineages <- list("Gamma (P.1 e P.1*)",
                        "Delta (B.1.617.2 e AY.*)",
                        "Omicron (BA.*, BE.*, BQ.*)",
                        "P.2","B.1.1.28","B.1.1.33","B.1.1")


                             
                             

lineage_colors <- c("Others" = "#8B8B8B", 
                    "P.1" = "#FF6666",
                    "P.1.2" = "#FF9999",
                    "P.1.4" = "#CC0000",
                    "P.1.7" = "#FF0000",
                    "P.1.14" = "#990000",
                    "AY.34.1.1" = "#66CC00",
                    "AY.43" = "#4C9900",
                    "AY.46.3" = "#336600",
                    "AY.99.1" = "#99FF99",
                    "AY.99.2" = "#CCFF99",
                    "AY.100" = "#66FFB2",
                    "AY.101" = "#99FF33",
                    "BA.1*" = "#99CCFF",
                    "BA.1.1*" = "#3399FF",
                    "BA.2*" = "#0066CC",
                    "BA.4*" = "#004C99",
                    "BA.5*" = "#6666FF",
                    "BA.5.2*"= "#B266FF",
                    "BQ.1*" = "#7F00FF")

states_initials <- c("AC", "AL", "AP", "AM", "BA", "CE", 
                     "DF", "ES", "GO", "MA", "MT", "MS", 
                     "MG", "PA", "PB", "PR", "PE", "PI", 
                     "RJ", "RN", "RS", "RO", "RR", "SC", 
                     "SP", "SE", "TO")

state_abb <- c(
  "acre" = "ac",
  "alagoas" = "al",
  "amapa" = "ap",
  "amazonas" = "am",
  "bahia" = "ba",
  "ceara" = "ce",
  "federal district" = "df",
  "espirito santo" = "es",
  "goias" = "go",
  "maranhao" = "ma",
  "mato grosso do sul" = "ms",
  "minas gerais" = "mg",
  "paraiba" = "pb",
  "parana" = "pr",
  "pernambuco" = "pe",
  "piaui" = "pi",
  "rio de janeiro" = "rj",
  "rio grande do norte" = "rn",
  "rio grande do sul" = "rs",
  "rondonia" = "ro",
  "roraima" = "rr",
  "santa catarina" = "sc",
  "sao paulo" = "sp",
  "sergipe" = "se",
  "tocantins" = "to"
)

## data
# gisaid collection date 06/02/2023
# to avoid load all metadata gisaid.tsv:
# tar -xf metadata_tsv_2023_02_06.tar.xz
# cut -f 1,5,6,7,9,10,11,12,13,14,20,21,22,23 metadata.tsv > metadata_cut.tsv
# grep "Brazil" metadata_cut.tsv > metadata_cut_br.tsv
br_gisaid <- read.csv("metadata_cut_br.tsv", header=TRUE, sep="\t")

br_gisaid$Collection.date <- as.Date(br_gisaid$Collection.date)
br_gisaid$Collection.date <- as.character(format(br_gisaid$Collection.date, "%Y-%m"))
# covid19.gov.br 06/02/2023
brasil_cases <- read.csv("covid_19_br.csv", header=TRUE, sep=";")

## filters
# remove genomes without collection date, without location, and non human genomes
br_gisaid  <- br_gisaid %>% 
  filter(!is.na(Collection.date) & !is.na(Location), Host == "Human")

# remove registers without region
brasil_cases <- brasil_cases %>%
  filter(regiao != "Brasil") %>%
  select(c(regiao, estado, data, municipio, casosNovos))
brasil_cases$data <- as.Date(brasil_cases$data, "%Y-%m-%d")

## Functions
merge_by_period <- function(df1, df2) {
  
  sub_df1 <- df1 %>% select(period, genome_count)
  merged_df <- sub_df1 %>% inner_join(df2, by = "period")
  merged_df <- subset(merged_df, select = -lineage) %>%
    unique()  %>% 
    group_by(period) %>% 
    mutate(genome_cases = (genome_count / casosNovos)*100)
  
  return(merged_df)
}

format_data <- function(df) {
  df <- df %>%
    filter(period %in% x_order_month) %>%
    arrange(factor(period, levels = x_order_month))
  
  return(df)
}

get_genomes_cases_proportion <- function(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10,
                                df11, df12, df13, df14, df15, df16, df17, df18, df19, df20,
                                df21, df22, df23, df24, df25, df26, df27) {

  dfs <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10,
                           df11, df12, df13, df14, df15, df16, df17, df18, df19, df20,
                           df21, df22, df23, df24, df25, df26, df27)
  
  sum_genome_count <- sapply(dfs, function(df) sum(df$genome_count))
  sum_new_cases <- sapply(dfs, function(df) sum(df$casosNovos))
  
  states_initials <- c("AC", "AL", "AP", "AM", "BA", "CE", "DF", "ES", "GO", "MA", "MT", "MS", "MG", "PA", "PB", "PR", "PE", "PI", "RJ", "RN", "RS", "RO", "RR", "SC", "SP", "SE", "TO")
  
  
  new_df <- data.frame(regiao = states_initials, sum_genome_count, sum_new_cases)
  new_df[28, 1] <- "Brasil"
  new_df[28, 2] <- sum(sum_genome_count)
  new_df[28, 3] <- sum(sum_new_cases)
  
  frac_genome_count <- (new_df$sum_genome_count / new_df$sum_genome_count[28])*100
  frac_new_cases <- (new_df$sum_new_cases / new_df$sum_new_cases[28])*100
  
  new_df$frac_genome_count <- frac_genome_count
  new_df$frac_new_cases <- frac_new_cases
  
  return(new_df)
}

filter_loc <- function(target_location, df) {
  filtered_df <- df %>%
    filter(grepl(target_location, Location, ignore.case = TRUE))
  
  return(filtered_df)
}

filter_loc_excep <- function(target_location,exclude, df) {
  filtered_df <- df %>%
    filter(grepl(target_location, Location, ignore.case = TRUE) &
             !grepl(exclude, Location, ignore.case = TRUE))
  
  return(filtered_df)
}

get_total_cases_period_state <- function(df, estado_arg, min_date, max_date) {
  df2 <- df %>% 
    filter(estado == estado_arg & 
             data >= as.Date(min_date, "%Y-%m-%d") & 
             data <= as.Date(max_date, "%Y-%m-%d" ))
  df2 <- df2[df2$municipio != "", ]
  df2$period <- format(df2$data, "%Y-%m")
  df3 <- aggregate(casosNovos ~ period, data=df2, sum)
  return(df3)
}

target_lineages_br_period <- function(df, targets, min_genomes) {
  df <- df[grep(paste(targets, collapse = "|"), df$Pango.lineage),] %>%
    group_by(Pango.lineage) %>%
    mutate(total_genomes = n()) %>%
    ungroup() %>%
    mutate(lineage_new = if_else(total_genomes < min_genomes, "Others", Pango.lineage)) %>%
    select(Collection.date, lineage_new) %>%
    group_by(Collection.date) %>%
    mutate(genome_count = n()) %>%
    group_by(Collection.date, lineage_new) %>%
    mutate(lineage_freq = n()/genome_count, lineage_count = n())
  
  lineages_period <- unique(select(df, Collection.date, lineage_new, lineage_freq, lineage_count, genome_count))
  names(lineages_period)[names(lineages_period)=="Collection.date"] <- "period"
  names(lineages_period)[names(lineages_period)=="lineage_new"] <- "lineage"
  lineages_period <- arrange(lineages_period, period)
  
  return(lineages_period)
}

plot_lineages <- function(df, lineage_range, scale) {
  unique_lineages <- unique(df$lineage)
  plot_colors <- lineage_colors[unique_lineages]
  
  df <- df %>%
    group_by(period) %>%
    arrange(-lineage_freq) %>% 
    mutate(total_freq = sum(lineage_freq),
           ymin = cumsum(lineage_freq) - lineage_freq,
           ymax = cumsum(lineage_freq)) %>%
    ungroup()
  
  lineage_plot <- ggplot(df, aes(x=factor(period, levels = lineage_range), y = lineage_freq, fill = lineage)) +
    geom_rect(aes(ymin = ymin, ymax = ymax, xmin = as.numeric(factor(period, levels = lineage_range))-0.4, xmax = as.numeric(factor(period, levels = lineage_range))+0.4), color="black") +
    geom_line(aes(y = genome_count/scale, group = 1), color = "black", size = 2, alpha = 0.5) +
    labs(x = "", y = "") +
    scale_y_continuous(sec.axis = sec_axis(~.*scale, name = "")) +
    scale_x_discrete(limits = lineage_range) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
          axis.line = element_line(size = 0.5, colour = "black"),
          axis.title = element_text(size = 13),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill = "transparent"),
          legend.box.just = "left",
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.spacing = unit(0.1, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_blank(),
          plot.title = element_text(size = 14, face = "bold")) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_fill_manual(values = plot_colors, labels = unique_lineages)
  
  return(lineage_plot)
}

## Ranges
# Gamma:
#  Origin in Brazil November 2020: https://www.nature.com/articles/s41591-021-01378-7
#  Replaced by Delta in July 2021: https://academic.oup.com/ve/article/8/1/veac024/6550638?login=false
# Delta:
#  Spread in Brazil April 2021: https://academic.oup.com/ve/article/8/1/veac024/6550638?login=false
#  No data at national scale of replacement by Omicron,
# estimates of AM state: December 2021: https://www.nature.com/articles/s41467-023-37541-6
# Omicron: Global spread November 2021: https://www.nature.com/articles/s41588-022-01033-y
gamma_range <- c("2020-11","2020-12", "2021-01","2021-02",
                 "2021-03", "2021-04", "2021-05", "2021-06",
                 "2021-07")

delta_range <- c("2021-04","2021-05","2021-06","2021-07",
                 "2021-08","2021-09","2021-10","2021-11",
                 "2021-12")

omicron_range <- c("2021-11", "2021-12","2022-01","2022-02", 
                   "2022-03", "2022-04","2022-05", "2022-06",
                   "2022-07","2022-08" )

study_range <- c("2020-11","2020-12", "2021-01","2021-02",
                 "2021-03", "2021-04", "2021-05", "2021-06",
                 "2021-07","2021-08","2021-09","2021-10",
                 "2021-11", "2021-12","2022-01","2022-02", 
                 "2022-03", "2022-04","2022-05", "2022-06",
                 "2022-07","2022-08")

for (state in names(state_abb)) {
    assign(paste0("gisaid_", state_abb[state]), 
           filter_loc(paste0("Brazil / ", state), br_gisaid))
}

gisaid_mt <- filter_loc_excep("Brazil / mato grosso", "Brazil / mato grosso do sul", br_gisaid)
gisaid_pa <- filter_loc_excep("Brazil / para","Brazil / parana|Brazil / paraiba", br_gisaid)


#Gamma
targets <- c("P.1","P.1.*")

#rs
gisaid_rs_lineages <- target_lineages_br_period(gisaid_rs, targets, 10)
gisaid_rs_lineages <- format_data(gisaid_rs_lineages)
rs_casos <- get_total_cases_period_state(brasil_cases, "RS", "2020-11-01", "2021-07-31")
rs_casos_genomas <- merge_by_period(gisaid_rs_lineages, rs_casos)
#rs-SM
gisaid_rs_sm <- gisaid_rs[grep("Santa Maria", gisaid_rs$Location, ignore.case = TRUE),]
gisaid_rs_sm_lineages <- target_lineages_br_period(gisaid_rs_sm, targets, 0)
gisaid_rs_sm_lineages <- format_data(gisaid_rs_sm_lineages)
rs_sm_casos_genomas <- merge_by_period(gisaid_rs_sm_lineages, rs_casos)

#sc
gisaid_sc_lineages <- target_lineages_br_period(gisaid_sc, targets, 10)
gisaid_sc_lineages <- format_data(gisaid_sc_lineages)
sc_casos <- get_total_cases_period_state(brasil_cases, "SC", "2020-11-01", "2021-07-31")
sc_casos_genomas <- merge_by_period(gisaid_sc_lineages, sc_casos)

#pr
gisaid_pr_lineages <- target_lineages_br_period(gisaid_pr, targets, 10)
gisaid_pr_lineages <- format_data(gisaid_pr_lineages)
pr_casos <- get_total_cases_period_state(brasil_cases, "PR", "2020-11-01", "2021-07-31")
pr_casos_genomas <- merge_by_period(gisaid_pr_lineages, pr_casos)

#ac
gisaid_ac_lineages <- target_lineages_br_period(gisaid_ac, targets, 10)
gisaid_ac_lineages <- format_data(gisaid_ac_lineages)
ac_casos <- get_total_cases_period_state(brasil_cases, "AC", "2020-11-01", "2021-07-31")
ac_casos_genomas <- merge_by_period(gisaid_ac_lineages, ac_casos)

#al
gisaid_al_lineages <- target_lineages_br_period(gisaid_al, targets, 10)
gisaid_al_lineages <- format_data(gisaid_al_lineages)
al_casos <- get_total_cases_period_state(brasil_cases, "AL", "2020-11-01", "2021-07-31")
al_casos_genomas <- merge_by_period(gisaid_al_lineages, al_casos)

#am
gisaid_am_lineages <- target_lineages_br_period(gisaid_am, targets, 10)
gisaid_am_lineages <- format_data(gisaid_am_lineages)
am_casos <- get_total_cases_period_state(brasil_cases, "AM", "2020-11-01", "2021-07-31")
am_casos_genomas <- merge_by_period(gisaid_am_lineages, am_casos)

#ap
gisaid_ap_lineages <- target_lineages_br_period(gisaid_ap, targets, 10)
gisaid_ap_lineages <- format_data(gisaid_ap_lineages)
ap_casos <- get_total_cases_period_state(brasil_cases, "AP", "2020-11-01", "2021-07-31")
ap_casos_genomas <- merge_by_period(gisaid_ap_lineages, ap_casos)

#ba
gisaid_ba_lineages <- target_lineages_br_period(gisaid_ba, targets, 10)
gisaid_ba_lineages <- format_data(gisaid_ba_lineages)
ba_casos <- get_total_cases_period_state(brasil_cases, "BA", "2020-11-01", "2021-07-31")
ba_casos_genomas <- merge_by_period(gisaid_ba_lineages, ba_casos)

#ce
gisaid_ce_lineages <- target_lineages_br_period(gisaid_ce, targets, 10)
gisaid_ce_lineages <- format_data(gisaid_ce_lineages)
ce_casos <- get_total_cases_period_state(brasil_cases, "CE", "2020-11-01", "2021-07-31")
ce_casos_genomas <- merge_by_period(gisaid_ce_lineages, ce_casos)

#df
gisaid_df_lineages <- target_lineages_br_period(gisaid_df, targets, 10)
gisaid_df_lineages <- format_data(gisaid_df_lineages)
df_casos <- get_total_cases_period_state(brasil_cases, "DF", "2020-11-01", "2021-07-31")
df_casos_genomas <- merge_by_period(gisaid_df_lineages, df_casos)

#es
gisaid_es_lineages <- target_lineages_br_period(gisaid_es, targets, 10)
gisaid_es_lineages <- format_data(gisaid_es_lineages)
es_casos <- get_total_cases_period_state(brasil_cases, "ES", "2020-11-01", "2021-07-31")
es_casos_genomas <- merge_by_period(gisaid_es_lineages, es_casos)

#go
gisaid_go_lineages <- target_lineages_br_period(gisaid_go, targets, 10)
gisaid_go_lineages <- format_data(gisaid_go_lineages)
go_casos <- get_total_cases_period_state(brasil_cases, "GO", "2020-11-01", "2021-07-31")
go_casos_genomas <- merge_by_period(gisaid_go_lineages, go_casos)

#ma
gisaid_ma_lineages <- target_lineages_br_period(gisaid_ma, targets, 10)
gisaid_ma_lineages <- format_data(gisaid_ma_lineages)
ma_casos <- get_total_cases_period_state(brasil_cases, "MA", "2020-11-01", "2021-07-31")
ma_casos_genomas <- merge_by_period(gisaid_ma_lineages, ma_casos)

#mt
gisaid_mt_lineages <- target_lineages_br_period(gisaid_mt, targets, 10)
gisaid_mt_lineages <- format_data(gisaid_mt_lineages)
mt_casos <- get_total_cases_period_state(brasil_cases, "MT", "2020-11-01", "2021-07-31")
mt_casos_genomas <- merge_by_period(gisaid_mt_lineages, mt_casos)

#ms
gisaid_ms_lineages <- target_lineages_br_period(gisaid_ms, targets, 10)
gisaid_ms_lineages <- format_data(gisaid_ms_lineages)
ms_casos <- get_total_cases_period_state(brasil_cases, "MS", "2020-11-01", "2021-07-31")
ms_casos_genomas <- merge_by_period(gisaid_ms_lineages, ms_casos)

#mg
gisaid_mg_lineages <- target_lineages_br_period(gisaid_mg, targets, 10)
gisaid_mg_lineages <- format_data(gisaid_mg_lineages)
mg_casos <- get_total_cases_period_state(brasil_cases, "MG", "2020-11-01", "2021-07-31")
mg_casos_genomas <- merge_by_period(gisaid_mg_lineages, mg_casos)

#pa
gisaid_pa_lineages <- target_lineages_br_period(gisaid_pa, targets, 10)
gisaid_pa_lineages <- format_data(gisaid_pa_lineages)
pa_casos <- get_total_cases_period_state(brasil_cases, "PA", "2020-11-01", "2021-07-31")
pa_casos_genomas <- merge_by_period(gisaid_pa_lineages, pa_casos)

#pb
gisaid_pb_lineages <- target_lineages_br_period(gisaid_pb, targets, 10)
gisaid_pb_lineages <- format_data(gisaid_pb_lineages)
pb_casos <- get_total_cases_period_state(brasil_cases, "PB", "2020-11-01", "2021-07-31")
pb_casos_genomas <- merge_by_period(gisaid_pb_lineages, pb_casos)

#rj
gisaid_rj_lineages <- target_lineages_br_period(gisaid_rj, targets, 10)
gisaid_rj_lineages <- format_data(gisaid_rj_lineages)
rj_casos <- get_total_cases_period_state(brasil_cases, "RJ", "2020-11-01", "2021-07-31")
rj_casos_genomas <- merge_by_period(gisaid_rj_lineages, rj_casos)

#sp
gisaid_sp_lineages <- target_lineages_br_period(gisaid_sp, targets, 10)
gisaid_sp_lineages <- format_data(gisaid_sp_lineages)
sp_casos <- get_total_cases_period_state(brasil_cases, "SP", "2020-11-01", "2021-07-31")
sp_casos_genomas <- merge_by_period(gisaid_sp_lineages, sp_casos)

#rn
gisaid_rn_lineages <- target_lineages_br_period(gisaid_rn, targets, 10)
gisaid_rn_lineages <- format_data(gisaid_rn_lineages)
rn_casos <- get_total_cases_period_state(brasil_cases, "RN", "2020-11-01", "2021-07-31")
rn_casos_genomas <- merge_by_period(gisaid_rn_lineages, rn_casos)

#pi
gisaid_pi_lineages <- target_lineages_br_period(gisaid_pi, targets, 10)
gisaid_pi_lineages <- format_data(gisaid_pi_lineages)
pi_casos <- get_total_cases_period_state(brasil_cases, "PI", "2020-11-01", "2021-07-31")
pi_casos_genomas <- merge_by_period(gisaid_pi_lineages, pi_casos)

#pe
gisaid_pe_lineages <- target_lineages_br_period(gisaid_pe, targets, 10)
gisaid_pe_lineages <- format_data(gisaid_pe_lineages)
pe_casos <- get_total_cases_period_state(brasil_cases, "PE", "2020-11-01", "2021-07-31")
pe_casos_genomas <- merge_by_period(gisaid_pe_lineages, pe_casos)

#se
gisaid_se_lineages <- target_lineages_br_period(gisaid_se, targets, 10)
gisaid_se_lineages <- format_data(gisaid_se_lineages)
se_casos <- get_total_cases_period_state(brasil_cases, "SE", "2020-11-01", "2021-07-31")
se_casos_genomas <- merge_by_period(gisaid_se_lineages, se_casos)

#ro
gisaid_ro_lineages <- target_lineages_br_period(gisaid_ro, targets, 10)
gisaid_ro_lineages <- format_data(gisaid_ro_lineages)
ro_casos <- get_total_cases_period_state(brasil_cases, "RO", "2020-11-01", "2021-07-31")
ro_casos_genomas <- merge_by_period(gisaid_ro_lineages, ro_casos)

#rr
gisaid_rr_lineages <- target_lineages_br_period(gisaid_rr, targets, 10)
gisaid_rr_lineages <- format_data(gisaid_rr_lineages)
rr_casos <- get_total_cases_period_state(brasil_cases, "RR", "2020-11-01", "2021-07-31")
rr_casos_genomas <- merge_by_period(gisaid_rr_lineages, rr_casos)

#to
gisaid_to_lineages <- target_lineages_br_period(gisaid_to, targets, 10)
gisaid_to_lineages <- format_data(gisaid_to_lineages)
to_casos <- get_total_cases_period_state(brasil_cases, "TO", "2020-11-01", "2021-07-31")
to_casos_genomas <- merge_by_period(gisaid_to_lineages, to_casos)

#br
gisaid_br <- filter_loc("Brazil", br_gisaid)
gisaid_br_lineages <- target_lineages_br_period(gisaid_br, targets, 600)

sum_genomes_cases <- get_genomes_cases_proportion(ac_casos_genomas,al_casos_genomas,ap_casos_genomas,
                                           am_casos_genomas,ba_casos_genomas,ce_casos_genomas,
                                           df_casos_genomas,es_casos_genomas,go_casos_genomas,
                                           ma_casos_genomas,mt_casos_genomas,ms_casos_genomas,
                                           mg_casos_genomas,pa_casos_genomas,pb_casos_genomas,
                                           pr_casos_genomas,pe_casos_genomas,pi_casos_genomas,
                                           rj_casos_genomas,rn_casos_genomas,rs_casos_genomas,
                                           ro_casos_genomas,rr_casos_genomas,sc_casos_genomas,
                                           sp_casos_genomas,se_casos_genomas,to_casos_genomas)
sum_genomes_cases <- filter(sum_genomes_cases, regiao != "Brasil")

cor_gamma <- cor(sum_genomes_cases$sum_genome_count, sum_genomes_cases$sum_new_cases)

sum_genomes_cases$genomas_sampling <- round(sum_genomes_cases$frac_new_cases * 30)
write.table(sum_genomes_cases, file = "genomas_sampling_gamma.tsv", sep = "\t", row.names = FALSE)

# plots
rs_sm_gamma_lineage_plot <- plot_lineages(gisaid_rs_sm_lineages, gamma_range, 15)
rs_gamma_lineage_plot <- plot_lineages(gisaid_rs_lineages, gamma_range, 300)
br_gamma_lineage_plot <- plot_lineages(gisaid_br_lineages, gamma_range, 12500)

#Delta
targets <- c("B.1.617.2","AY.*")
#rs
gisaid_rs_lineages <- target_lineages_br_period(gisaid_rs, targets, 15)
gisaid_rs_lineages <- format_data(gisaid_rs_lineages)
rs_casos <- get_total_cases_period_state(brasil_cases, "RS", "2021-04-01", "2021-12-31")
rs_casos_genomas <- merge_by_period(gisaid_rs_lineages, rs_casos)
#rs-SM
gisaid_rs_sm <- gisaid_rs[grep("Santa Maria", gisaid_rs$Location, ignore.case = TRUE),]
gisaid_rs_sm_lineages <- target_lineages_br_period(gisaid_rs_sm, targets, 1)
gisaid_rs_sm_lineages <- format_data(gisaid_rs_sm_lineages)
rs_sm_casos_genomas <- merge_by_period(gisaid_rs_sm_lineages, rs_casos)

#sc
gisaid_sc_lineages <- target_lineages_br_period(gisaid_sc, targets, 15)
gisaid_sc_lineages <- format_data(gisaid_sc_lineages)
sc_casos <- get_total_cases_period_state(brasil_cases, "SC", "2021-04-01", "2021-12-31")
sc_casos_genomas <- merge_by_period(gisaid_sc_lineages, sc_casos)

#pr
gisaid_pr_lineages <- target_lineages_br_period(gisaid_pr, targets, 15)
gisaid_pr_lineages <- format_data(gisaid_pr_lineages)
pr_casos <- get_total_cases_period_state(brasil_cases, "PR", "2021-04-01", "2021-12-31")
pr_casos_genomas <- merge_by_period(gisaid_pr_lineages, pr_casos)

#ac
gisaid_ac_lineages <- target_lineages_br_period(gisaid_ac, targets, 15)
gisaid_ac_lineages <- format_data(gisaid_ac_lineages)
ac_casos <- get_total_cases_period_state(brasil_cases, "AC", "2021-04-01", "2021-12-31")
ac_casos_genomas <- merge_by_period(gisaid_ac_lineages, ac_casos)

#al
gisaid_al_lineages <- target_lineages_br_period(gisaid_al, targets, 15)
gisaid_al_lineages <- format_data(gisaid_al_lineages)
al_casos <- get_total_cases_period_state(brasil_cases, "AL", "2021-04-01", "2021-12-31")
al_casos_genomas <- merge_by_period(gisaid_al_lineages, al_casos)

#am
gisaid_am_lineages <- target_lineages_br_period(gisaid_am, targets, 15)
gisaid_am_lineages <- format_data(gisaid_am_lineages)
am_casos <- get_total_cases_period_state(brasil_cases, "AM", "2021-04-01", "2021-12-31")
am_casos_genomas <- merge_by_period(gisaid_am_lineages, am_casos)

#ap
gisaid_ap_lineages <- target_lineages_br_period(gisaid_ap, targets, 15)
gisaid_ap_lineages <- format_data(gisaid_ap_lineages)
ap_casos <- get_total_cases_period_state(brasil_cases, "AP", "2021-04-01", "2021-12-31")
ap_casos_genomas <- merge_by_period(gisaid_ap_lineages, ap_casos)

#ba
gisaid_ba_lineages <- target_lineages_br_period(gisaid_ba, targets, 15)
gisaid_ba_lineages <- format_data(gisaid_ba_lineages)
ba_casos <- get_total_cases_period_state(brasil_cases, "BA", "2021-04-01", "2021-12-31")
ba_casos_genomas <- merge_by_period(gisaid_ba_lineages, ba_casos)

#ce
gisaid_ce_lineages <- target_lineages_br_period(gisaid_ce, targets, 15)
gisaid_ce_lineages <- format_data(gisaid_ce_lineages)
ce_casos <- get_total_cases_period_state(brasil_cases, "CE", "2021-04-01", "2021-12-31")
ce_casos_genomas <- merge_by_period(gisaid_ce_lineages, ce_casos)

#df
gisaid_df_lineages <- target_lineages_br_period(gisaid_df, targets, 15)
gisaid_df_lineages <- format_data(gisaid_df_lineages)
df_casos <- get_total_cases_period_state(brasil_cases, "DF", "2021-04-01", "2021-12-31")
df_casos_genomas <- merge_by_period(gisaid_df_lineages, df_casos)

#es
gisaid_es_lineages <- target_lineages_br_period(gisaid_es, targets, 15)
gisaid_es_lineages <- format_data(gisaid_es_lineages)
es_casos <- get_total_cases_period_state(brasil_cases, "ES", "2021-04-01", "2021-12-31")
es_casos_genomas <- merge_by_period(gisaid_es_lineages, es_casos)

#go
gisaid_go_lineages <- target_lineages_br_period(gisaid_go, targets, 15)
gisaid_go_lineages <- format_data(gisaid_go_lineages)
go_casos <- get_total_cases_period_state(brasil_cases, "GO", "2021-04-01", "2021-12-31")
go_casos_genomas <- merge_by_period(gisaid_go_lineages, go_casos)

#ma
gisaid_ma_lineages <- target_lineages_br_period(gisaid_ma, targets, 15)
gisaid_ma_lineages <- format_data(gisaid_ma_lineages)
ma_casos <- get_total_cases_period_state(brasil_cases, "MA", "2021-04-01", "2021-12-31")
ma_casos_genomas <- merge_by_period(gisaid_ma_lineages, ma_casos)

#mt
gisaid_mt_lineages <- target_lineages_br_period(gisaid_mt, targets, 15)
gisaid_mt_lineages <- format_data(gisaid_mt_lineages)
mt_casos <- get_total_cases_period_state(brasil_cases, "MT", "2021-04-01", "2021-12-31")
mt_casos_genomas <- merge_by_period(gisaid_mt_lineages, mt_casos)

#ms
gisaid_ms_lineages <- target_lineages_br_period(gisaid_ms, targets, 15)
gisaid_ms_lineages <- format_data(gisaid_ms_lineages)
ms_casos <- get_total_cases_period_state(brasil_cases, "MS", "2021-04-01", "2021-12-31")
ms_casos_genomas <- merge_by_period(gisaid_ms_lineages, ms_casos)

#mg
gisaid_mg_lineages <- target_lineages_br_period(gisaid_mg, targets, 15)
gisaid_mg_lineages <- format_data(gisaid_mg_lineages)
mg_casos <- get_total_cases_period_state(brasil_cases, "MG", "2021-04-01", "2021-12-31")
mg_casos_genomas <- merge_by_period(gisaid_mg_lineages, mg_casos)

#pa
gisaid_pa_lineages <- target_lineages_br_period(gisaid_pa, targets, 15)
gisaid_pa_lineages <- format_data(gisaid_pa_lineages)
pa_casos <- get_total_cases_period_state(brasil_cases, "PA", "2021-04-01", "2021-12-31")
pa_casos_genomas <- merge_by_period(gisaid_pa_lineages, pa_casos)

#pb
gisaid_pb_lineages <- target_lineages_br_period(gisaid_pb, targets, 15)
gisaid_pb_lineages <- format_data(gisaid_pb_lineages)
pb_casos <- get_total_cases_period_state(brasil_cases, "PB", "2021-04-01", "2021-12-31")
pb_casos_genomas <- merge_by_period(gisaid_pb_lineages, pb_casos)

#rj
gisaid_rj_lineages <- target_lineages_br_period(gisaid_rj, targets, 15)
gisaid_rj_lineages <- format_data(gisaid_rj_lineages)
rj_casos <- get_total_cases_period_state(brasil_cases, "RJ", "2021-04-01", "2021-12-31")
rj_casos_genomas <- merge_by_period(gisaid_rj_lineages, rj_casos)

#sp
gisaid_sp_lineages <- target_lineages_br_period(gisaid_sp, targets, 15)
gisaid_sp_lineages <- format_data(gisaid_sp_lineages)
sp_casos <- get_total_cases_period_state(brasil_cases, "SP", "2021-04-01", "2021-12-31")
sp_casos_genomas <- merge_by_period(gisaid_sp_lineages, sp_casos)

#rn
gisaid_rn_lineages <- target_lineages_br_period(gisaid_rn, targets, 15)
gisaid_rn_lineages <- format_data(gisaid_rn_lineages)
rn_casos <- get_total_cases_period_state(brasil_cases, "RN", "2021-04-01", "2021-12-31")
rn_casos_genomas <- merge_by_period(gisaid_rn_lineages, rn_casos)

#pi
gisaid_pi_lineages <- target_lineages_br_period(gisaid_pi, targets, 15)
gisaid_pi_lineages <- format_data(gisaid_pi_lineages)
pi_casos <- get_total_cases_period_state(brasil_cases, "PI", "2021-04-01", "2021-12-31")
pi_casos_genomas <- merge_by_period(gisaid_pi_lineages, pi_casos)

#pe
gisaid_pe_lineages <- target_lineages_br_period(gisaid_pe, targets, 15)
gisaid_pe_lineages <- format_data(gisaid_pe_lineages)
pe_casos <- get_total_cases_period_state(brasil_cases, "PE", "2021-04-01", "2021-12-31")
pe_casos_genomas <- merge_by_period(gisaid_pe_lineages, pe_casos)

#se
gisaid_se_lineages <- target_lineages_br_period(gisaid_se, targets, 15)
gisaid_se_lineages <- format_data(gisaid_se_lineages)
se_casos <- get_total_cases_period_state(brasil_cases, "SE", "2021-04-01", "2021-12-31")
se_casos_genomas <- merge_by_period(gisaid_se_lineages, se_casos)

#ro
gisaid_ro_lineages <- target_lineages_br_period(gisaid_ro, targets, 15)
gisaid_ro_lineages <- format_data(gisaid_ro_lineages)
ro_casos <- get_total_cases_period_state(brasil_cases, "RO", "2021-04-01", "2021-12-31")
ro_casos_genomas <- merge_by_period(gisaid_ro_lineages, ro_casos)

#rr
gisaid_rr_lineages <- target_lineages_br_period(gisaid_rr, targets, 15)
gisaid_rr_lineages <- format_data(gisaid_rr_lineages)
rr_casos <- get_total_cases_period_state(brasil_cases, "RR", "2021-04-01", "2021-12-31")
rr_casos_genomas <- merge_by_period(gisaid_rr_lineages, rr_casos)

#to
gisaid_to_lineages <- target_lineages_br_period(gisaid_to, targets, 15)
gisaid_to_lineages <- format_data(gisaid_to_lineages)
to_casos <- get_total_cases_period_state(brasil_cases, "TO", "2021-04-01", "2021-12-31")
to_casos_genomas <- merge_by_period(gisaid_to_lineages, to_casos)

#br
gisaid_br <- filter_loc("Brazil", br_gisaid)
gisaid_br_lineages <- target_lineages_br_period(gisaid_br, targets, 1500)

sum_genomes_cases <- get_genomes_cases_proportion(ac_casos_genomas,al_casos_genomas,ap_casos_genomas,
                                            am_casos_genomas,ba_casos_genomas,ce_casos_genomas,
                                            df_casos_genomas,es_casos_genomas,go_casos_genomas,
                                            ma_casos_genomas,mt_casos_genomas,ms_casos_genomas,
                                            mg_casos_genomas,pa_casos_genomas,pb_casos_genomas,
                                            pr_casos_genomas,pe_casos_genomas,pi_casos_genomas,
                                            rj_casos_genomas,rn_casos_genomas,rs_casos_genomas,
                                            ro_casos_genomas,rr_casos_genomas,sc_casos_genomas,
                                            sp_casos_genomas,se_casos_genomas,to_casos_genomas)
sum_genomes_cases <- filter(sum_genomes_cases, regiao != "Brasil")
cor_delta <- cor(sum_genomes_cases$sum_genome_count, sum_genomes_cases$sum_new_cases)

sum_genomes_cases$genomas_sampling <- round(sum_genomes_cases$frac_new_cases * 30)
write.table(sum_genomes_cases, file = "genomas_sampling_delta.tsv", sep = "\t", row.names = FALSE)

# plots
rs_sm_gamma_lineage_plot <- plot_lineages(gisaid_rs_sm_lineages, delta_range, 2)
rs_delta_lineage_plot <- plot_lineages(gisaid_rs_lineages, delta_range, 200)
br_delta_lineage_plot <- plot_lineages(gisaid_br_lineages, delta_range, 12500)

#Omicron
omicron_lineages_br_period <- function(df, targets, min_genomes) {
  df <- df[grep(paste(targets, collapse = "|"), df$Pango.lineage),] %>%
    mutate(Pango.lineage = case_when(
      str_detect(Pango.lineage, "^BA\\.1\\.1(\\.\\S*)?$") ~ "BA.1.1*",
      str_detect(Pango.lineage, "^BA\\.1(\\.\\S*)?$") ~ "BA.1*",
      str_detect(Pango.lineage, "^BA\\.2(\\.\\S*)?$") ~ "BA.2*",
      str_detect(Pango.lineage, "^BA\\.4(\\.\\S*)?$") ~ "BA.4*",
      str_detect(Pango.lineage, "^BA\\.5\\.2(\\.\\S*)?$") ~ "BA.5.2*",
      str_detect(Pango.lineage, "^BA\\.5(\\.\\S*)?$") ~ "BA.5*",
      str_detect(Pango.lineage, "^BQ\\.1(\\.\\S*)?$") ~ "BQ.1*",
      str_detect(Pango.lineage, "^BE\\.1(\\.\\S*)?$") ~ "BE.1*",
      TRUE ~ Pango.lineage
    )) %>%
      group_by(Pango.lineage) %>%
      mutate(total_genomes = n()) %>%
      ungroup() %>%
      mutate(lineage_new = if_else(total_genomes < min_genomes, "Others", Pango.lineage)) %>%
      select(Collection.date, lineage_new) %>%
      group_by(Collection.date) %>%
      mutate(genome_count = n()) %>%
      group_by(Collection.date, lineage_new) %>%
      mutate(lineage_freq = n()/genome_count, lineage_count = n())
    
    lineages_period <- unique(select(df, Collection.date, lineage_new, lineage_freq, lineage_count, genome_count))
    names(lineages_period)[names(lineages_period)=="Collection.date"] <- "period"
    names(lineages_period)[names(lineages_period)=="lineage_new"] <- "lineage"
    lineages_period <- arrange(lineages_period, period)
    
    return(lineages_period)
}

targets <- c("BA.*","BE.*","BQ.*","DL.*")
#rs
gisaid_rs_lineages <- omicron_lineages_br_period(gisaid_rs, targets, 200)
gisaid_rs_lineages <- format_data(gisaid_rs_lineages)
rs_casos <- get_total_cases_period_state(brasil_cases, "RS", "2021-11-01", "2022-08-31")
rs_casos_genomas <- merge_by_period(gisaid_rs_lineages, rs_casos)
#rs-SM
gisaid_rs_sm <- gisaid_rs[grep("Santa Maria", gisaid_rs$Location, ignore.case = TRUE),]
gisaid_rs_sm_lineages <- omicron_lineages_br_period(gisaid_rs_sm, targets, 0)
gisaid_rs_sm_lineages <- format_data(gisaid_rs_sm_lineages)
rs_sm_casos_genomas <- merge_by_period(gisaid_rs_sm_lineages, rs_casos)
#sc
gisaid_sc_lineages <- omicron_lineages_br_period(gisaid_sc, targets, 200)
gisaid_sc_lineages <- format_data(gisaid_sc_lineages)
sc_casos <- get_total_cases_period_state(brasil_cases, "SC", "2021-11-01", "2022-08-31")
sc_casos_genomas <- merge_by_period(gisaid_sc_lineages, sc_casos)

#pr
gisaid_pr_lineages <- omicron_lineages_br_period(gisaid_pr, targets, 200)
gisaid_pr_lineages <- format_data(gisaid_pr_lineages)
pr_casos <- get_total_cases_period_state(brasil_cases, "PR", "2021-11-01", "2022-08-31")
pr_casos_genomas <- merge_by_period(gisaid_pr_lineages, pr_casos)

#ac
gisaid_ac_lineages <- omicron_lineages_br_period(gisaid_ac, targets, 200)
gisaid_ac_lineages <- format_data(gisaid_ac_lineages)
ac_casos <- get_total_cases_period_state(brasil_cases, "AC", "2021-11-01", "2022-08-31")
ac_casos_genomas <- merge_by_period(gisaid_ac_lineages, ac_casos)

#al
gisaid_al_lineages <- omicron_lineages_br_period(gisaid_al, targets, 200)
gisaid_al_lineages <- format_data(gisaid_al_lineages)
al_casos <- get_total_cases_period_state(brasil_cases, "AL", "2021-11-01", "2022-08-31")
al_casos_genomas <- merge_by_period(gisaid_al_lineages, al_casos)

#am
gisaid_am_lineages <- omicron_lineages_br_period(gisaid_am, targets, 200)
gisaid_am_lineages <- format_data(gisaid_am_lineages)
am_casos <- get_total_cases_period_state(brasil_cases, "AM", "2021-11-01", "2022-08-31")
am_casos_genomas <- merge_by_period(gisaid_am_lineages, am_casos)

#ap
gisaid_ap_lineages <- omicron_lineages_br_period(gisaid_ap, targets, 200)
gisaid_ap_lineages <- format_data(gisaid_ap_lineages)
ap_casos <- get_total_cases_period_state(brasil_cases, "AP", "2021-11-01", "2022-08-31")
ap_casos_genomas <- merge_by_period(gisaid_ap_lineages, ap_casos)

#ba
gisaid_ba_lineages <- omicron_lineages_br_period(gisaid_ba, targets, 200)
gisaid_ba_lineages <- format_data(gisaid_ba_lineages)
ba_casos <- get_total_cases_period_state(brasil_cases, "BA", "2021-11-01", "2022-08-31")
ba_casos_genomas <- merge_by_period(gisaid_ba_lineages, ba_casos)

#ce
gisaid_ce_lineages <- omicron_lineages_br_period(gisaid_ce, targets, 200)
gisaid_ce_lineages <- format_data(gisaid_ce_lineages)
ce_casos <- get_total_cases_period_state(brasil_cases, "CE", "2021-11-01", "2022-08-31")
ce_casos_genomas <- merge_by_period(gisaid_ce_lineages, ce_casos)

#df
gisaid_df_lineages <- omicron_lineages_br_period(gisaid_df, targets, 200)
gisaid_df_lineages <- format_data(gisaid_df_lineages)
df_casos <- get_total_cases_period_state(brasil_cases, "DF", "2021-11-01", "2022-08-31")
df_casos_genomas <- merge_by_period(gisaid_df_lineages, df_casos)

#es
gisaid_es_lineages <- omicron_lineages_br_period(gisaid_es, targets, 200)
gisaid_es_lineages <- format_data(gisaid_es_lineages)
es_casos <- get_total_cases_period_state(brasil_cases, "ES", "2021-11-01", "2022-08-31")
es_casos_genomas <- merge_by_period(gisaid_es_lineages, es_casos)

#go
gisaid_go_lineages <- omicron_lineages_br_period(gisaid_go, targets, 200)
gisaid_go_lineages <- format_data(gisaid_go_lineages)
go_casos <- get_total_cases_period_state(brasil_cases, "GO", "2021-11-01", "2022-08-31")
go_casos_genomas <- merge_by_period(gisaid_go_lineages, go_casos)

#ma
gisaid_ma_lineages <- omicron_lineages_br_period(gisaid_ma, targets, 200)
gisaid_ma_lineages <- format_data(gisaid_ma_lineages)
ma_casos <- get_total_cases_period_state(brasil_cases, "MA", "2021-11-01", "2022-08-31")
ma_casos_genomas <- merge_by_period(gisaid_ma_lineages, ma_casos)

#mt
gisaid_mt_lineages <- omicron_lineages_br_period(gisaid_mt, targets, 200)
gisaid_mt_lineages <- format_data(gisaid_mt_lineages)
mt_casos <- get_total_cases_period_state(brasil_cases, "MT", "2021-11-01", "2022-08-31")
mt_casos_genomas <- merge_by_period(gisaid_mt_lineages, mt_casos)

#ms
gisaid_ms_lineages <- omicron_lineages_br_period(gisaid_ms, targets, 200)
gisaid_ms_lineages <- format_data(gisaid_ms_lineages)
ms_casos <- get_total_cases_period_state(brasil_cases, "MS", "2021-11-01", "2022-08-31")
ms_casos_genomas <- merge_by_period(gisaid_ms_lineages, ms_casos)

#mg
gisaid_mg_lineages <- omicron_lineages_br_period(gisaid_mg, targets, 200)
gisaid_mg_lineages <- format_data(gisaid_mg_lineages)
mg_casos <- get_total_cases_period_state(brasil_cases, "MG", "2021-11-01", "2022-08-31")
mg_casos_genomas <- merge_by_period(gisaid_mg_lineages, mg_casos)

#pa
gisaid_pa_lineages <- omicron_lineages_br_period(gisaid_pa, targets, 200)
gisaid_pa_lineages <- format_data(gisaid_pa_lineages)
pa_casos <- get_total_cases_period_state(brasil_cases, "PA", "2021-11-01", "2022-08-31")
pa_casos_genomas <- merge_by_period(gisaid_pa_lineages, pa_casos)

#pb
gisaid_pb_lineages <- omicron_lineages_br_period(gisaid_pb, targets, 200)
gisaid_pb_lineages <- format_data(gisaid_pb_lineages)
pb_casos <- get_total_cases_period_state(brasil_cases, "PB", "2021-11-01", "2022-08-31")
pb_casos_genomas <- merge_by_period(gisaid_pb_lineages, pb_casos)

#rj
gisaid_rj_lineages <- omicron_lineages_br_period(gisaid_rj, targets, 200)
gisaid_rj_lineages <- format_data(gisaid_rj_lineages)
rj_casos <- get_total_cases_period_state(brasil_cases, "RJ", "2021-11-01", "2022-08-31")
rj_casos_genomas <- merge_by_period(gisaid_rj_lineages, rj_casos)

#sp
gisaid_sp_lineages <- omicron_lineages_br_period(gisaid_sp, targets, 200)
gisaid_sp_lineages <- format_data(gisaid_sp_lineages)
sp_casos <- get_total_cases_period_state(brasil_cases, "SP", "2021-11-01", "2022-08-31")
sp_casos_genomas <- merge_by_period(gisaid_sp_lineages, sp_casos)

#rn
gisaid_rn_lineages <- omicron_lineages_br_period(gisaid_rn, targets, 200)
gisaid_rn_lineages <- format_data(gisaid_rn_lineages)
rn_casos <- get_total_cases_period_state(brasil_cases, "RN", "2021-11-01", "2022-08-31")
rn_casos_genomas <- merge_by_period(gisaid_rn_lineages, rn_casos)

#pi
gisaid_pi_lineages <- omicron_lineages_br_period(gisaid_pi, targets, 200)
gisaid_pi_lineages <- format_data(gisaid_pi_lineages)
pi_casos <- get_total_cases_period_state(brasil_cases, "PI", "2021-11-01", "2022-08-31")
pi_casos_genomas <- merge_by_period(gisaid_pi_lineages, pi_casos)

#pe
gisaid_pe_lineages <- omicron_lineages_br_period(gisaid_pe, targets, 200)
gisaid_pe_lineages <- format_data(gisaid_pe_lineages)
pe_casos <- get_total_cases_period_state(brasil_cases, "PE", "2021-11-01", "2022-08-31")
pe_casos_genomas <- merge_by_period(gisaid_pe_lineages, pe_casos)

#se
gisaid_se_lineages <- omicron_lineages_br_period(gisaid_se, targets, 200)
gisaid_se_lineages <- format_data(gisaid_se_lineages)
se_casos <- get_total_cases_period_state(brasil_cases, "SE", "2021-11-01", "2022-08-31")
se_casos_genomas <- merge_by_period(gisaid_se_lineages, se_casos)

#ro
gisaid_ro_lineages <- omicron_lineages_br_period(gisaid_ro, targets, 200)
gisaid_ro_lineages <- format_data(gisaid_ro_lineages)
ro_casos <- get_total_cases_period_state(brasil_cases, "RO", "2021-11-01", "2022-08-31")
ro_casos_genomas <- merge_by_period(gisaid_ro_lineages, ro_casos)

#rr
gisaid_rr_lineages <- omicron_lineages_br_period(gisaid_rr, targets, 200)
gisaid_rr_lineages <- format_data(gisaid_rr_lineages)
rr_casos <- get_total_cases_period_state(brasil_cases, "RR", "2021-11-01", "2022-08-31")
rr_casos_genomas <- merge_by_period(gisaid_rr_lineages, rr_casos)

#tp
gisaid_to_lineages <- omicron_lineages_br_period(gisaid_to, targets, 200)
gisaid_to_lineages <- format_data(gisaid_to_lineages)
to_casos <- get_total_cases_period_state(brasil_cases, "TO", "2021-11-01", "2022-08-31")
to_casos_genomas <- merge_by_period(gisaid_to_lineages, to_casos)

gisaid_br <- filter_loc("Brazil", br_gisaid)
gisaid_br_lineages <- omicron_lineages_br_period(gisaid_br, targets, 1500)

sum_genomes_cases <- get_genomes_cases_proportion(ac_casos_genomas,al_casos_genomas,ap_casos_genomas,
                                            am_casos_genomas,ba_casos_genomas,ce_casos_genomas,
                                            df_casos_genomas,es_casos_genomas,go_casos_genomas,
                                            ma_casos_genomas,mt_casos_genomas,ms_casos_genomas,
                                            mg_casos_genomas,pa_casos_genomas,pb_casos_genomas,
                                            pr_casos_genomas,pe_casos_genomas,pi_casos_genomas,
                                            rj_casos_genomas,rn_casos_genomas,rs_casos_genomas,
                                            ro_casos_genomas,rr_casos_genomas,sc_casos_genomas,
                                            sp_casos_genomas,se_casos_genomas,to_casos_genomas)
sum_genomes_cases <- filter(sum_genomes_cases, regiao != "Brasil")
cor_omicron <- cor(sum_genomes_cases$sum_genome_count, sum_genomes_cases$sum_new_cases)


sum_genomes_cases$genomas_sampling <- round(sum_genomes_cases$frac_new_cases * 30)
write.table(sum_genomes_cases, file = "genomas_sampling_omicron.tsv", sep = "\t", row.names = FALSE)

rs_sm_omicron_lineage_plot <- plot_lineages(gisaid_rs_sm_lineages, omicron_range, 25)
rs_omicron_lineage_plot <- plot_lineages(gisaid_rs_lineages, omicron_range, 650)
br_omicron_lineage_plot <- plot_lineages(gisaid_br_lineages, omicron_range, 25000)

convert_others <- function(df, min_genomes) {
  df <- df %>%
    group_by(lineage) %>%
    mutate(total_genomes = n()) %>%
    ungroup() %>%
    mutate(lineage_new = if_else(total_genomes < min_genomes, "Others", lineage)) %>%
    select(Collection.date, lineage_new) %>%
    group_by(Collection.date) %>%
    mutate(genome_count = n()) %>%
    group_by(Collection.date, lineage_new) %>%
    mutate(lineage_freq = n()/genome_count, lineage_count = n())
  
  lineages_period <- unique(select(df, Collection.date, lineage_new, lineage_freq, lineage_count, genome_count))
  names(lineages_period)[names(lineages_period)=="Collection.date"] <- "period"
  names(lineages_period)[names(lineages_period)=="lineage_new"] <- "lineage"
  lineages_period <- arrange(lineages_period, period)
  
  return(lineages_period)
}

target_lineages <- list("Gamma (P.1 and P.1*)",
                        "Delta (B.1.617.2 and AY.*)",
                        "Omicron (BA.*, BE.*, BQ.* and DL.*)",
                        "P.2","B.1.1.28","B.1.1.33","B.1.1")

variant_colors <- c("B.1.1"='#FFFF99',"Others"='#C0C0C0', "P.2" = "994C00",
                    "B.1.1.28"='#FFCC99',"B.1.1.33"='#FF9933',
                    "Gamma (P.1 and P.1.*)"="#FF6666",
                    "Delta (B.1.617.2 and AY.*)"='#B2FF66',
                    "Omicron (BA.*, BE.*, BQ.* and DL.*)" = "#99CCFF")

convert_to_variants <- function(df) {
  df <- df %>%
    mutate(lineage = case_when(
      str_detect(Pango.lineage, "P\\.1.*") ~ "Gamma (P.1 and P.1.*)",
      str_detect(Pango.lineage, "AY.*") ~ "Delta (B.1.617.2 and AY.*)",
      Pango.lineage == "B.1.617.2" ~ "Delta (B.1.617.2 and AY.*)",
      str_detect(Pango.lineage, "BA\\..*|BE\\..*|BQ\\..*|DL\\..*") ~ "Omicron (BA.*, BE.*, BQ.* and DL.*)",
      Pango.lineage %in% unlist(target_lineages) ~ Pango.lineage,
      TRUE ~ "Others"
    ))
  
  return(df)
}

plot_variants <- function(df, lineage_range, scale) {
  unique_lineages <- unique(df$lineage)
  plot_colors <- variant_colors[unique_lineages]
  
  df <- df %>%
    group_by(period) %>%
    arrange(-lineage_freq) %>% 
    mutate(total_freq = sum(lineage_freq),
           ymin = cumsum(lineage_freq) - lineage_freq,
           ymax = cumsum(lineage_freq)) %>%
    ungroup()
  
  lineage_plot <- ggplot(df, aes(x=factor(period, levels = lineage_range), y = lineage_freq, fill = lineage)) +
    geom_rect(aes(ymin = ymin, ymax = ymax, xmin = as.numeric(factor(period, levels = lineage_range))-0.4, xmax = as.numeric(factor(period, levels = lineage_range))+0.4), color="black") +
    geom_line(aes(y = genome_count/scale, group = 1), color = "black", size = 2, alpha = 0.5) +
    labs(x = "Period", y = "Lineage Frequency") +
    scale_y_continuous(sec.axis = sec_axis(~.*scale, name = "Genomes (n)")) +
    scale_x_discrete(limits = lineage_range) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
          axis.line = element_line(size = 0.5, colour = "black"),
          axis.title = element_text(size = 13),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill = "transparent"),
          legend.box.just = "left",
          legend.box.margin = margin(0, 0, 0, 0),
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.spacing = unit(0.1, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_blank(),
          plot.title = element_text(size = 14, face = "bold")) +
    scale_fill_manual(values = plot_colors, labels = unique_lineages)
  
  return(lineage_plot)
}

br_gisaid_variants <- convert_to_variants(br_gisaid)
br_gisaid_variants <- convert_others(br_gisaid_variants,1)
rs_gisaid_variants <- convert_to_variants(gisaid_rs)
rs_gisaid_variants <- convert_others(rs_gisaid_variants,1)
gisaid_br_all_plot <- plot_variants(br_gisaid_variants, x_order_month, 26000)
gisaid_rs_all_plot <- plot_variants(rs_gisaid_variants, x_order_month, 700)


library(gridExtra)

grid_upper <- grid.arrange(gisaid_br_all_plot, gisaid_rs_all_plot, nrow = 2)
grid_lower <- grid.arrange(br_gamma_lineage_plot, rs_gamma_lineage_plot, br_delta_lineage_plot, rs_delta_lineage_plot, br_omicron_lineage_plot, rs_omicron_lineage_plot, nrow = 3, ncol = 2)
grid_final <- grid.arrange(grid_upper, grid_lower, nrow = 2)
grid_final
save.image(file = "2023-06-26_explorer.RData")



## Subclades timetree analysis
#ibge: https://censo2022.ibge.gov.br/panorama/ access: 07 Out 2023.
ibge_pop <- read.csv("~/pessoal_2/sars_rs/ibge.csv", header=TRUE)
br_gisaid <- read.csv("~/pessoal/sars_rs/2022-10-05_SARS_RS/data/metadata_cut_br.tsv", header=TRUE, sep="\t")
br_gisaid$Collection.date <- as.Date(br_gisaid$Collection.date)
br_gisaid  <- br_gisaid %>% 
  filter(!is.na(Collection.date) & !is.na(Location), Host == "Human")

for (state in names(state_abb)) {
  assign(paste0("gisaid_", state_abb[state]), 
         filter_loc(paste0("Brazil / ", state), br_gisaid))
}

gisaid_mt <- filter_loc_excep("Brazil / mato grosso", "Brazil / mato grosso do sul", br_gisaid)
gisaid_pa <- filter_loc_excep("Brazil / para","Brazil / parana|Brazil / paraiba", br_gisaid)

get_subclades_data <- function(df, data_minima, data_maxima, linhagem_alvo) {
  # Converter Collection.date para o formato Date
  df$Collection.date <- as.Date(df$Collection.date)
  
  # Filtrar o dataframe pelas datas mínima e máxima
  df_filtrado <- df %>%
    filter(Collection.date >= data_minima & Collection.date <= data_maxima)
  
  # Criar uma coluna de período no formato AAAA-MM
  df_filtrado <- df_filtrado %>%
    mutate(period = format(Collection.date, "%Y-%m"))
  
  # Contar o número de linhas da linhagem alvo por período
  linhas_linhagem_alvo <- df_filtrado %>%
    filter(Pango.lineage == linhagem_alvo) %>%
    group_by(period) %>%
    summarise(linhagem_alvo = n())
  
  # Contar o número total de linhas por período
  total_linhas <- df_filtrado %>%
    group_by(period) %>%
    summarise(total = n())
  
  # Juntar os resultados
  resultado <- right_join(linhas_linhagem_alvo, total_linhas, by = "period")
  
  # Calcular a frequência da linhagem alvo em relação ao total por período
  resultado <- resultado %>%
    mutate(frequencia = linhagem_alvo / total)
  resultado$frequencia <- ifelse(is.na(resultado$frequencia), 0, resultado$frequencia)
  
  return(resultado)
}

get_prop_cases_target_lineage <- function(df1, df2) {
  # Converter a coluna "periodo" para o mesmo tipo de dados em ambos os dataframes
  df1$period <- as.character(df1$period)
  df2$period <- as.character(df2$period)
  
  # Unir os dataframes pelo período
  novo_df <- merge(df1, df2, by = "period", all.x = TRUE)
  
  # Calcular o número de casos da linhagem alvo
  novo_df$numero_casos_linhagem_alvo <- novo_df$frequencia * novo_df$casosNovos
  novo_df$numero_casos_outras_linhagens <- novo_df$casosNovos - novo_df$numero_casos_linhagem_alvo
  soma_linhagem_alvo <- sum(novo_df$numero_casos_linhagem_alvo, na.rm = TRUE)
  soma_outras_linhagens <- sum(novo_df$numero_casos_outras_linhagens, na.rm = TRUE)
  soma_genomas_linhagem_alvo <- sum(novo_df$linhagem_alvo, na.rm = TRUE)
  soma_total_genomas <- sum(novo_df$total, na.rm = TRUE)
  
  # Retornar como uma lista
  resultado <- c(soma_linhagem_alvo, soma_outras_linhagens, soma_genomas_linhagem_alvo, soma_total_genomas)
  return(resultado)
}

get_genomes_cases_proportion <- function(list1, list2, list3, list4, list5, list6, list7, list8, list9, list10,
                                         list11, list12, list13, list14, list15, list16, list17, list18, list19, list20,
                                         list21, list22, list23, list24, list25, list26, list27) {
  
  lists <- list(list1, list2, list3, list4, list5, list6, list7, list8, list9, list10,
              list11, list12, list13, list14, list15, list16, list17, list18, list19, list20,
              list21, list22, list23, list24, list25, list26, list27)
  
  states_initials <- c("AC", "AL", "AP", "AM", "BA", "CE", "DF", "ES", "GO", "MA", "MT", "MS", "MG", "PA", "PB", "PR", "PE", "PI", "RJ", "RN", "RS", "RO", "RR", "SC", "SP", "SE", "TO")
  
  for (i in seq_along(lists)) {
    lists[[i]] <- as.integer(round(lists[[i]]))
  }
  
  # Criar um dataframe
  df <- data.frame(states_initials, 
                   target_lineage_cases= sapply(lists, function(x) x[1]),
                   other_lineage_cases = sapply(lists, function(x) x[2]),
                   target_lineage_genomes = sapply(lists, function(x) x[3]),
                   total_lineage_genomes = sapply(lists, function(x) x[4]))
  
  return(df)
}


## Gamma

## P1.14
min_date_p1.14 <- "2020-11-24"
max_date_p1.14 <- "2021-07-30"
target <- "P.1.14"
#rs
gisaid_rs_lineages <- get_subclades_data(gisaid_rs, min_date_p1.14,max_date_p1.14, target)
gisaid_rs_lineages <- format_data(gisaid_rs_lineages)
rs_casos <- get_total_cases_period_state(brasil_cases, "RS", min_date_p1.14, max_date_p1.14)
rs_casos_genomas <- get_prop_cases_target_lineage(gisaid_rs_lineages, rs_casos)

#sc
gisaid_sc_lineages <- get_subclades_data(gisaid_sc, min_date_p1.14,max_date_p1.14, target)
gisaid_sc_lineages <- format_data(gisaid_sc_lineages)
sc_casos <- get_total_cases_period_state(brasil_cases, "SC", min_date_p1.14, max_date_p1.14)
sc_casos_genomas <- get_prop_cases_target_lineage(gisaid_sc_lineages, sc_casos)

#pr
gisaid_pr_lineages <- get_subclades_data(gisaid_pr, min_date_p1.14,max_date_p1.14, target)
gisaid_pr_lineages <- format_data(gisaid_pr_lineages)
pr_casos <- get_total_cases_period_state(brasil_cases, "PR", min_date_p1.14, max_date_p1.14)
pr_casos_genomas <- get_prop_cases_target_lineage(gisaid_pr_lineages, pr_casos)

#ac
gisaid_ac_lineages <- get_subclades_data(gisaid_ac, min_date_p1.14,max_date_p1.14, target)
gisaid_ac_lineages <- format_data(gisaid_ac_lineages)
ac_casos <- get_total_cases_period_state(brasil_cases, "AC", min_date_p1.14, max_date_p1.14)
ac_casos_genomas <- get_prop_cases_target_lineage(gisaid_ac_lineages, ac_casos)

#al
gisaid_al_lineages <- get_subclades_data(gisaid_al, min_date_p1.14,max_date_p1.14, target)
gisaid_al_lineages <- format_data(gisaid_al_lineages)
al_casos <- get_total_cases_period_state(brasil_cases, "AL", min_date_p1.14, max_date_p1.14)
al_casos_genomas <- get_prop_cases_target_lineage(gisaid_al_lineages, al_casos)

#am
gisaid_am_lineages <- get_subclades_data(gisaid_am, min_date_p1.14,max_date_p1.14, target)
gisaid_am_lineages <- format_data(gisaid_am_lineages)
am_casos <- get_total_cases_period_state(brasil_cases, "AM", min_date_p1.14, max_date_p1.14)
am_casos_genomas <- get_prop_cases_target_lineage(gisaid_am_lineages, am_casos)

#ap
gisaid_ap_lineages <- get_subclades_data(gisaid_ap, min_date_p1.14,max_date_p1.14, target)
gisaid_ap_lineages <- format_data(gisaid_ap_lineages)
ap_casos <- get_total_cases_period_state(brasil_cases, "AP", min_date_p1.14, max_date_p1.14)
ap_casos_genomas <- get_prop_cases_target_lineage(gisaid_ap_lineages, ap_casos)

#ba
gisaid_ba_lineages <- get_subclades_data(gisaid_ba, min_date_p1.14,max_date_p1.14, target)
gisaid_ba_lineages <- format_data(gisaid_ba_lineages)
ba_casos <- get_total_cases_period_state(brasil_cases, "BA", min_date_p1.14, max_date_p1.14)
ba_casos_genomas <- get_prop_cases_target_lineage(gisaid_ba_lineages, ba_casos)

#ce
gisaid_ce_lineages <- get_subclades_data(gisaid_ce, min_date_p1.14,max_date_p1.14, target)
gisaid_ce_lineages <- format_data(gisaid_ce_lineages)
ce_casos <- get_total_cases_period_state(brasil_cases, "CE", min_date_p1.14, max_date_p1.14)
ce_casos_genomas <- get_prop_cases_target_lineage(gisaid_ce_lineages, ce_casos)

#df
gisaid_df_lineages <- get_subclades_data(gisaid_df, min_date_p1.14,max_date_p1.14, target)
gisaid_df_lineages <- format_data(gisaid_df_lineages)
df_casos <- get_total_cases_period_state(brasil_cases, "DF", min_date_p1.14, max_date_p1.14)
df_casos_genomas <- get_prop_cases_target_lineage(gisaid_df_lineages, df_casos)

#es
gisaid_es_lineages <- get_subclades_data(gisaid_es, min_date_p1.14,max_date_p1.14, target)
gisaid_es_lineages <- format_data(gisaid_es_lineages)
es_casos <- get_total_cases_period_state(brasil_cases, "ES", min_date_p1.14, max_date_p1.14)
es_casos_genomas <- get_prop_cases_target_lineage(gisaid_es_lineages, es_casos)

#go
gisaid_go_lineages <- get_subclades_data(gisaid_go, min_date_p1.14,max_date_p1.14, target)
gisaid_go_lineages <- format_data(gisaid_go_lineages)
go_casos <- get_total_cases_period_state(brasil_cases, "GO", min_date_p1.14, max_date_p1.14)
go_casos_genomas <- get_prop_cases_target_lineage(gisaid_go_lineages, go_casos)

#ma
gisaid_ma_lineages <- get_subclades_data(gisaid_ma, min_date_p1.14,max_date_p1.14, target)
gisaid_ma_lineages <- format_data(gisaid_ma_lineages)
ma_casos <- get_total_cases_period_state(brasil_cases, "MA", min_date_p1.14, max_date_p1.14)
ma_casos_genomas <- get_prop_cases_target_lineage(gisaid_ma_lineages, ma_casos)

#mt
gisaid_mt_lineages <- get_subclades_data(gisaid_mt, min_date_p1.14,max_date_p1.14, target)
gisaid_mt_lineages <- format_data(gisaid_mt_lineages)
mt_casos <- get_total_cases_period_state(brasil_cases, "MT", min_date_p1.14, max_date_p1.14)
mt_casos_genomas <- get_prop_cases_target_lineage(gisaid_mt_lineages, mt_casos)

#ms
gisaid_ms_lineages <- get_subclades_data(gisaid_ms, min_date_p1.14,max_date_p1.14, target)
gisaid_ms_lineages <- format_data(gisaid_ms_lineages)
ms_casos <- get_total_cases_period_state(brasil_cases, "MS", min_date_p1.14, max_date_p1.14)
ms_casos_genomas <- get_prop_cases_target_lineage(gisaid_ms_lineages, ms_casos)

#mg
gisaid_mg_lineages <- get_subclades_data(gisaid_mg, min_date_p1.14,max_date_p1.14, target)
gisaid_mg_lineages <- format_data(gisaid_mg_lineages)
mg_casos <- get_total_cases_period_state(brasil_cases, "MG", min_date_p1.14, max_date_p1.14)
mg_casos_genomas <- get_prop_cases_target_lineage(gisaid_mg_lineages, mg_casos)

#pa
gisaid_pa_lineages <- get_subclades_data(gisaid_pa, min_date_p1.14,max_date_p1.14, target)
gisaid_pa_lineages <- format_data(gisaid_pa_lineages)
pa_casos <- get_total_cases_period_state(brasil_cases, "PA", min_date_p1.14, max_date_p1.14)
pa_casos_genomas <- get_prop_cases_target_lineage(gisaid_pa_lineages, pa_casos)

#pb
gisaid_pb_lineages <- get_subclades_data(gisaid_pb, min_date_p1.14,max_date_p1.14, target)
gisaid_pb_lineages <- format_data(gisaid_pb_lineages)
pb_casos <- get_total_cases_period_state(brasil_cases, "PB", min_date_p1.14, max_date_p1.14)
pb_casos_genomas <- get_prop_cases_target_lineage(gisaid_pb_lineages, pb_casos)

#rj
gisaid_rj_lineages <- get_subclades_data(gisaid_rj, min_date_p1.14,max_date_p1.14, target)
gisaid_rj_lineages <- format_data(gisaid_rj_lineages)
rj_casos <- get_total_cases_period_state(brasil_cases, "RJ", min_date_p1.14, max_date_p1.14)
rj_casos_genomas <- get_prop_cases_target_lineage(gisaid_rj_lineages, rj_casos)

#sp
gisaid_sp_lineages <- get_subclades_data(gisaid_sp, min_date_p1.14,max_date_p1.14, target)
gisaid_sp_lineages <- format_data(gisaid_sp_lineages)
sp_casos <- get_total_cases_period_state(brasil_cases, "SP", min_date_p1.14, max_date_p1.14)
sp_casos_genomas <- get_prop_cases_target_lineage(gisaid_sp_lineages, sp_casos)

#rn
gisaid_rn_lineages <- get_subclades_data(gisaid_rn, min_date_p1.14,max_date_p1.14, target)
gisaid_rn_lineages <- format_data(gisaid_rn_lineages)
rn_casos <- get_total_cases_period_state(brasil_cases, "RN", min_date_p1.14, max_date_p1.14)
rn_casos_genomas <- get_prop_cases_target_lineage(gisaid_rn_lineages, rn_casos)

#pi
gisaid_pi_lineages <- get_subclades_data(gisaid_pi, min_date_p1.14,max_date_p1.14, target)
gisaid_pi_lineages <- format_data(gisaid_pi_lineages)
pi_casos <- get_total_cases_period_state(brasil_cases, "PI", min_date_p1.14, max_date_p1.14)
pi_casos_genomas <- get_prop_cases_target_lineage(gisaid_pi_lineages, pi_casos)

#pe
gisaid_pe_lineages <- get_subclades_data(gisaid_pe, min_date_p1.14,max_date_p1.14, target)
gisaid_pe_lineages <- format_data(gisaid_pe_lineages)
pe_casos <- get_total_cases_period_state(brasil_cases, "PE", min_date_p1.14, max_date_p1.14)
pe_casos_genomas <- get_prop_cases_target_lineage(gisaid_pe_lineages, pe_casos)

#se
gisaid_se_lineages <- get_subclades_data(gisaid_se, min_date_p1.14,max_date_p1.14, target)
gisaid_se_lineages <- format_data(gisaid_se_lineages)
se_casos <- get_total_cases_period_state(brasil_cases, "SE", min_date_p1.14, max_date_p1.14)
se_casos_genomas <- get_prop_cases_target_lineage(gisaid_se_lineages, se_casos)

#ro
gisaid_ro_lineages <- get_subclades_data(gisaid_ro, min_date_p1.14,max_date_p1.14, target)
gisaid_ro_lineages <- format_data(gisaid_ro_lineages)
ro_casos <- get_total_cases_period_state(brasil_cases, "RO", min_date_p1.14, max_date_p1.14)
ro_casos_genomas <- get_prop_cases_target_lineage(gisaid_ro_lineages, ro_casos)

#rr
gisaid_rr_lineages <- get_subclades_data(gisaid_rr, min_date_p1.14,max_date_p1.14, target)
gisaid_rr_lineages <- format_data(gisaid_rr_lineages)
rr_casos <- get_total_cases_period_state(brasil_cases, "RR", min_date_p1.14, max_date_p1.14)
rr_casos_genomas <- get_prop_cases_target_lineage(gisaid_rr_lineages, rr_casos)

#to
gisaid_to_lineages <- get_subclades_data(gisaid_to, min_date_p1.14,max_date_p1.14, target)
gisaid_to_lineages <- format_data(gisaid_to_lineages)
to_casos <- get_total_cases_period_state(brasil_cases, "TO", min_date_p1.14, max_date_p1.14)
to_casos_genomas <- get_prop_cases_target_lineage(gisaid_to_lineages, to_casos)

cases_by_state <- get_genomes_cases_proportion(ac_casos_genomas,al_casos_genomas,ap_casos_genomas,
                                               am_casos_genomas,ba_casos_genomas,ce_casos_genomas,
                                               df_casos_genomas,es_casos_genomas,go_casos_genomas,
                                               ma_casos_genomas,mt_casos_genomas,ms_casos_genomas,
                                               mg_casos_genomas,pa_casos_genomas,pb_casos_genomas,
                                               pr_casos_genomas,pe_casos_genomas,pi_casos_genomas,
                                               rj_casos_genomas,rn_casos_genomas,rs_casos_genomas,
                                               ro_casos_genomas,rr_casos_genomas,sc_casos_genomas,
                                               sp_casos_genomas,se_casos_genomas,to_casos_genomas)
cases_by_state
p1.14_data <- left_join(cases_by_state, ibge_pop, by = "states_initials")
p1.14_data <- setNames(p1.14_data, c("state", "cases_prop_p.1.14", "cases_prop_others","genomes_p.1.14","genomes_total","pop_2022"))
p1.14_data$cases_total <- p1.14_data$cases_prop_p.1.14 + p1.14_data$cases_prop_others
write.table(p1.14_data, file = "p1_14_traits.tsv", sep = "\t", row.names = FALSE)

## P1.2
min_date_p1.2 <- "2021-01-15"
max_date_p1.2 <- "2021-07-31"
target <- "P.1.2"
#rs
gisaid_rs_lineages <- get_subclades_data(gisaid_rs, min_date_p1.2,max_date_p1.2, target)
gisaid_rs_lineages <- format_data(gisaid_rs_lineages)
rs_casos <- get_total_cases_period_state(brasil_cases, "RS", min_date_p1.2, max_date_p1.2)
rs_casos_genomas <- get_prop_cases_target_lineage(gisaid_rs_lineages, rs_casos)

#sc
gisaid_sc_lineages <- get_subclades_data(gisaid_sc, min_date_p1.2,max_date_p1.2, target)
gisaid_sc_lineages <- format_data(gisaid_sc_lineages)
sc_casos <- get_total_cases_period_state(brasil_cases, "SC", min_date_p1.2, max_date_p1.2)
sc_casos_genomas <- get_prop_cases_target_lineage(gisaid_sc_lineages, sc_casos)

#pr
gisaid_pr_lineages <- get_subclades_data(gisaid_pr, min_date_p1.2,max_date_p1.2, target)
gisaid_pr_lineages <- format_data(gisaid_pr_lineages)
pr_casos <- get_total_cases_period_state(brasil_cases, "PR", min_date_p1.2, max_date_p1.2)
pr_casos_genomas <- get_prop_cases_target_lineage(gisaid_pr_lineages, pr_casos)

#ac
gisaid_ac_lineages <- get_subclades_data(gisaid_ac, min_date_p1.2,max_date_p1.2, target)
gisaid_ac_lineages <- format_data(gisaid_ac_lineages)
ac_casos <- get_total_cases_period_state(brasil_cases, "AC", min_date_p1.2, max_date_p1.2)
ac_casos_genomas <- get_prop_cases_target_lineage(gisaid_ac_lineages, ac_casos)

#al
gisaid_al_lineages <- get_subclades_data(gisaid_al, min_date_p1.2,max_date_p1.2, target)
gisaid_al_lineages <- format_data(gisaid_al_lineages)
al_casos <- get_total_cases_period_state(brasil_cases, "AL", min_date_p1.2, max_date_p1.2)
al_casos_genomas <- get_prop_cases_target_lineage(gisaid_al_lineages, al_casos)

#am
gisaid_am_lineages <- get_subclades_data(gisaid_am, min_date_p1.2,max_date_p1.2, target)
gisaid_am_lineages <- format_data(gisaid_am_lineages)
am_casos <- get_total_cases_period_state(brasil_cases, "AM", min_date_p1.2, max_date_p1.2)
am_casos_genomas <- get_prop_cases_target_lineage(gisaid_am_lineages, am_casos)

#ap
gisaid_ap_lineages <- get_subclades_data(gisaid_ap, min_date_p1.2,max_date_p1.2, target)
gisaid_ap_lineages <- format_data(gisaid_ap_lineages)
ap_casos <- get_total_cases_period_state(brasil_cases, "AP", min_date_p1.2, max_date_p1.2)
ap_casos_genomas <- get_prop_cases_target_lineage(gisaid_ap_lineages, ap_casos)

#ba
gisaid_ba_lineages <- get_subclades_data(gisaid_ba, min_date_p1.2,max_date_p1.2, target)
gisaid_ba_lineages <- format_data(gisaid_ba_lineages)
ba_casos <- get_total_cases_period_state(brasil_cases, "BA", min_date_p1.2, max_date_p1.2)
ba_casos_genomas <- get_prop_cases_target_lineage(gisaid_ba_lineages, ba_casos)

#ce
gisaid_ce_lineages <- get_subclades_data(gisaid_ce, min_date_p1.2,max_date_p1.2, target)
gisaid_ce_lineages <- format_data(gisaid_ce_lineages)
ce_casos <- get_total_cases_period_state(brasil_cases, "CE", min_date_p1.2, max_date_p1.2)
ce_casos_genomas <- get_prop_cases_target_lineage(gisaid_ce_lineages, ce_casos)

#df
gisaid_df_lineages <- get_subclades_data(gisaid_df, min_date_p1.2,max_date_p1.2, target)
gisaid_df_lineages <- format_data(gisaid_df_lineages)
df_casos <- get_total_cases_period_state(brasil_cases, "DF", min_date_p1.2, max_date_p1.2)
df_casos_genomas <- get_prop_cases_target_lineage(gisaid_df_lineages, df_casos)

#es
gisaid_es_lineages <- get_subclades_data(gisaid_es, min_date_p1.2,max_date_p1.2, target)
gisaid_es_lineages <- format_data(gisaid_es_lineages)
es_casos <- get_total_cases_period_state(brasil_cases, "ES", min_date_p1.2, max_date_p1.2)
es_casos_genomas <- get_prop_cases_target_lineage(gisaid_es_lineages, es_casos)

#go
gisaid_go_lineages <- get_subclades_data(gisaid_go, min_date_p1.2,max_date_p1.2, target)
gisaid_go_lineages <- format_data(gisaid_go_lineages)
go_casos <- get_total_cases_period_state(brasil_cases, "GO", min_date_p1.2, max_date_p1.2)
go_casos_genomas <- get_prop_cases_target_lineage(gisaid_go_lineages, go_casos)

#ma
gisaid_ma_lineages <- get_subclades_data(gisaid_ma, min_date_p1.2,max_date_p1.2, target)
gisaid_ma_lineages <- format_data(gisaid_ma_lineages)
ma_casos <- get_total_cases_period_state(brasil_cases, "MA", min_date_p1.2, max_date_p1.2)
ma_casos_genomas <- get_prop_cases_target_lineage(gisaid_ma_lineages, ma_casos)

#mt
gisaid_mt_lineages <- get_subclades_data(gisaid_mt, min_date_p1.2,max_date_p1.2, target)
gisaid_mt_lineages <- format_data(gisaid_mt_lineages)
mt_casos <- get_total_cases_period_state(brasil_cases, "MT", min_date_p1.2, max_date_p1.2)
mt_casos_genomas <- get_prop_cases_target_lineage(gisaid_mt_lineages, mt_casos)

#ms
gisaid_ms_lineages <- get_subclades_data(gisaid_ms, min_date_p1.2,max_date_p1.2, target)
gisaid_ms_lineages <- format_data(gisaid_ms_lineages)
ms_casos <- get_total_cases_period_state(brasil_cases, "MS", min_date_p1.2, max_date_p1.2)
ms_casos_genomas <- get_prop_cases_target_lineage(gisaid_ms_lineages, ms_casos)

#mg
gisaid_mg_lineages <- get_subclades_data(gisaid_mg, min_date_p1.2,max_date_p1.2, target)
gisaid_mg_lineages <- format_data(gisaid_mg_lineages)
mg_casos <- get_total_cases_period_state(brasil_cases, "MG", min_date_p1.2, max_date_p1.2)
mg_casos_genomas <- get_prop_cases_target_lineage(gisaid_mg_lineages, mg_casos)

#pa
gisaid_pa_lineages <- get_subclades_data(gisaid_pa, min_date_p1.2,max_date_p1.2, target)
gisaid_pa_lineages <- format_data(gisaid_pa_lineages)
pa_casos <- get_total_cases_period_state(brasil_cases, "PA", min_date_p1.2, max_date_p1.2)
pa_casos_genomas <- get_prop_cases_target_lineage(gisaid_pa_lineages, pa_casos)

#pb
gisaid_pb_lineages <- get_subclades_data(gisaid_pb, min_date_p1.2,max_date_p1.2, target)
gisaid_pb_lineages <- format_data(gisaid_pb_lineages)
pb_casos <- get_total_cases_period_state(brasil_cases, "PB", min_date_p1.2, max_date_p1.2)
pb_casos_genomas <- get_prop_cases_target_lineage(gisaid_pb_lineages, pb_casos)

#rj
gisaid_rj_lineages <- get_subclades_data(gisaid_rj, min_date_p1.2,max_date_p1.2, target)
gisaid_rj_lineages <- format_data(gisaid_rj_lineages)
rj_casos <- get_total_cases_period_state(brasil_cases, "RJ", min_date_p1.2, max_date_p1.2)
rj_casos_genomas <- get_prop_cases_target_lineage(gisaid_rj_lineages, rj_casos)

#sp
gisaid_sp_lineages <- get_subclades_data(gisaid_sp, min_date_p1.2,max_date_p1.2, target)
gisaid_sp_lineages <- format_data(gisaid_sp_lineages)
sp_casos <- get_total_cases_period_state(brasil_cases, "SP", min_date_p1.2, max_date_p1.2)
sp_casos_genomas <- get_prop_cases_target_lineage(gisaid_sp_lineages, sp_casos)

#rn
gisaid_rn_lineages <- get_subclades_data(gisaid_rn, min_date_p1.2,max_date_p1.2, target)
gisaid_rn_lineages <- format_data(gisaid_rn_lineages)
rn_casos <- get_total_cases_period_state(brasil_cases, "RN", min_date_p1.2, max_date_p1.2)
rn_casos_genomas <- get_prop_cases_target_lineage(gisaid_rn_lineages, rn_casos)

#pi
gisaid_pi_lineages <- get_subclades_data(gisaid_pi, min_date_p1.2,max_date_p1.2, target)
gisaid_pi_lineages <- format_data(gisaid_pi_lineages)
pi_casos <- get_total_cases_period_state(brasil_cases, "PI", min_date_p1.2, max_date_p1.2)
pi_casos_genomas <- get_prop_cases_target_lineage(gisaid_pi_lineages, pi_casos)

#pe
gisaid_pe_lineages <- get_subclades_data(gisaid_pe, min_date_p1.2,max_date_p1.2, target)
gisaid_pe_lineages <- format_data(gisaid_pe_lineages)
pe_casos <- get_total_cases_period_state(brasil_cases, "PE", min_date_p1.2, max_date_p1.2)
pe_casos_genomas <- get_prop_cases_target_lineage(gisaid_pe_lineages, pe_casos)

#se
gisaid_se_lineages <- get_subclades_data(gisaid_se, min_date_p1.2,max_date_p1.2, target)
gisaid_se_lineages <- format_data(gisaid_se_lineages)
se_casos <- get_total_cases_period_state(brasil_cases, "SE", min_date_p1.2, max_date_p1.2)
se_casos_genomas <- get_prop_cases_target_lineage(gisaid_se_lineages, se_casos)

#ro
gisaid_ro_lineages <- get_subclades_data(gisaid_ro, min_date_p1.2,max_date_p1.2, target)
gisaid_ro_lineages <- format_data(gisaid_ro_lineages)
ro_casos <- get_total_cases_period_state(brasil_cases, "RO", min_date_p1.2, max_date_p1.2)
ro_casos_genomas <- get_prop_cases_target_lineage(gisaid_ro_lineages, ro_casos)

#rr
gisaid_rr_lineages <- get_subclades_data(gisaid_rr, min_date_p1.2,max_date_p1.2, target)
gisaid_rr_lineages <- format_data(gisaid_rr_lineages)
rr_casos <- get_total_cases_period_state(brasil_cases, "RR", min_date_p1.2, max_date_p1.2)
rr_casos_genomas <- get_prop_cases_target_lineage(gisaid_rr_lineages, rr_casos)

#to
gisaid_to_lineages <- get_subclades_data(gisaid_to, min_date_p1.2,max_date_p1.2, target)
gisaid_to_lineages <- format_data(gisaid_to_lineages)
to_casos <- get_total_cases_period_state(brasil_cases, "TO", min_date_p1.2, max_date_p1.2)
to_casos_genomas <- get_prop_cases_target_lineage(gisaid_to_lineages, to_casos)

cases_by_state <- get_genomes_cases_proportion(ac_casos_genomas,al_casos_genomas,ap_casos_genomas,
                                                  am_casos_genomas,ba_casos_genomas,ce_casos_genomas,
                                                  df_casos_genomas,es_casos_genomas,go_casos_genomas,
                                                  ma_casos_genomas,mt_casos_genomas,ms_casos_genomas,
                                                  mg_casos_genomas,pa_casos_genomas,pb_casos_genomas,
                                                  pr_casos_genomas,pe_casos_genomas,pi_casos_genomas,
                                                  rj_casos_genomas,rn_casos_genomas,rs_casos_genomas,
                                                  ro_casos_genomas,rr_casos_genomas,sc_casos_genomas,
                                                  sp_casos_genomas,se_casos_genomas,to_casos_genomas)
cases_by_state
p1.2_data <- left_join(cases_by_state, ibge_pop, by = "states_initials")
p1.2_data <- setNames(p1.2_data, c("state", "cases_prop_p.1.2", "cases_prop_others","genomes_p.1.2","genomes_total","pop_2022"))
p1.2_data$cases_total <- p1.2_data$cases_prop_p.1.2 + p1.2_data$cases_prop_others
write.table(p1.2_data, file = "p1_2_traits.tsv", sep = "\t", row.names = FALSE)

## P1.7
min_date_p1.7 <- "2021-01-15"
max_date_p1.7 <- "2021-07-31"
target <- "P.1.7"

#rs
gisaid_rs_lineages <- get_subclades_data(gisaid_rs, min_date_p1.7,max_date_p1.7, target)
gisaid_rs_lineages <- format_data(gisaid_rs_lineages)
rs_casos <- get_total_cases_period_state(brasil_cases, "RS", min_date_p1.7, max_date_p1.7)
rs_casos_genomas <- get_prop_cases_target_lineage(gisaid_rs_lineages, rs_casos)

#sc
gisaid_sc_lineages <- get_subclades_data(gisaid_sc, min_date_p1.7,max_date_p1.7, target)
gisaid_sc_lineages <- format_data(gisaid_sc_lineages)
sc_casos <- get_total_cases_period_state(brasil_cases, "SC", min_date_p1.7, max_date_p1.7)
sc_casos_genomas <- get_prop_cases_target_lineage(gisaid_sc_lineages, sc_casos)

#pr
gisaid_pr_lineages <- get_subclades_data(gisaid_pr, min_date_p1.7,max_date_p1.7, target)
gisaid_pr_lineages <- format_data(gisaid_pr_lineages)
pr_casos <- get_total_cases_period_state(brasil_cases, "PR", min_date_p1.7, max_date_p1.7)
pr_casos_genomas <- get_prop_cases_target_lineage(gisaid_pr_lineages, pr_casos)

#ac
gisaid_ac_lineages <- get_subclades_data(gisaid_ac, min_date_p1.7,max_date_p1.7, target)
gisaid_ac_lineages <- format_data(gisaid_ac_lineages)
ac_casos <- get_total_cases_period_state(brasil_cases, "AC", min_date_p1.7, max_date_p1.7)
ac_casos_genomas <- get_prop_cases_target_lineage(gisaid_ac_lineages, ac_casos)

#al
gisaid_al_lineages <- get_subclades_data(gisaid_al, min_date_p1.7,max_date_p1.7, target)
gisaid_al_lineages <- format_data(gisaid_al_lineages)
al_casos <- get_total_cases_period_state(brasil_cases, "AL", min_date_p1.7, max_date_p1.7)
al_casos_genomas <- get_prop_cases_target_lineage(gisaid_al_lineages, al_casos)

#am
gisaid_am_lineages <- get_subclades_data(gisaid_am, min_date_p1.7,max_date_p1.7, target)
gisaid_am_lineages <- format_data(gisaid_am_lineages)
am_casos <- get_total_cases_period_state(brasil_cases, "AM", min_date_p1.7, max_date_p1.7)
am_casos_genomas <- get_prop_cases_target_lineage(gisaid_am_lineages, am_casos)

#ap
gisaid_ap_lineages <- get_subclades_data(gisaid_ap, min_date_p1.7,max_date_p1.7, target)
gisaid_ap_lineages <- format_data(gisaid_ap_lineages)
ap_casos <- get_total_cases_period_state(brasil_cases, "AP", min_date_p1.7, max_date_p1.7)
ap_casos_genomas <- get_prop_cases_target_lineage(gisaid_ap_lineages, ap_casos)

#ba
gisaid_ba_lineages <- get_subclades_data(gisaid_ba, min_date_p1.7,max_date_p1.7, target)
gisaid_ba_lineages <- format_data(gisaid_ba_lineages)
ba_casos <- get_total_cases_period_state(brasil_cases, "BA", min_date_p1.7, max_date_p1.7)
ba_casos_genomas <- get_prop_cases_target_lineage(gisaid_ba_lineages, ba_casos)

#ce
gisaid_ce_lineages <- get_subclades_data(gisaid_ce, min_date_p1.7,max_date_p1.7, target)
gisaid_ce_lineages <- format_data(gisaid_ce_lineages)
ce_casos <- get_total_cases_period_state(brasil_cases, "CE", min_date_p1.7, max_date_p1.7)
ce_casos_genomas <- get_prop_cases_target_lineage(gisaid_ce_lineages, ce_casos)

#df
gisaid_df_lineages <- get_subclades_data(gisaid_df, min_date_p1.7,max_date_p1.7, target)
gisaid_df_lineages <- format_data(gisaid_df_lineages)
df_casos <- get_total_cases_period_state(brasil_cases, "DF", min_date_p1.7, max_date_p1.7)
df_casos_genomas <- get_prop_cases_target_lineage(gisaid_df_lineages, df_casos)

#es
gisaid_es_lineages <- get_subclades_data(gisaid_es, min_date_p1.7,max_date_p1.7, target)
gisaid_es_lineages <- format_data(gisaid_es_lineages)
es_casos <- get_total_cases_period_state(brasil_cases, "ES", min_date_p1.7, max_date_p1.7)
es_casos_genomas <- get_prop_cases_target_lineage(gisaid_es_lineages, es_casos)

#go
gisaid_go_lineages <- get_subclades_data(gisaid_go, min_date_p1.7,max_date_p1.7, target)
gisaid_go_lineages <- format_data(gisaid_go_lineages)
go_casos <- get_total_cases_period_state(brasil_cases, "GO", min_date_p1.7, max_date_p1.7)
go_casos_genomas <- get_prop_cases_target_lineage(gisaid_go_lineages, go_casos)

#ma
gisaid_ma_lineages <- get_subclades_data(gisaid_ma, min_date_p1.7,max_date_p1.7, target)
gisaid_ma_lineages <- format_data(gisaid_ma_lineages)
ma_casos <- get_total_cases_period_state(brasil_cases, "MA", min_date_p1.7, max_date_p1.7)
ma_casos_genomas <- get_prop_cases_target_lineage(gisaid_ma_lineages, ma_casos)

#mt
gisaid_mt_lineages <- get_subclades_data(gisaid_mt, min_date_p1.7,max_date_p1.7, target)
gisaid_mt_lineages <- format_data(gisaid_mt_lineages)
mt_casos <- get_total_cases_period_state(brasil_cases, "MT", min_date_p1.7, max_date_p1.7)
mt_casos_genomas <- get_prop_cases_target_lineage(gisaid_mt_lineages, mt_casos)

#ms
gisaid_ms_lineages <- get_subclades_data(gisaid_ms, min_date_p1.7,max_date_p1.7, target)
gisaid_ms_lineages <- format_data(gisaid_ms_lineages)
ms_casos <- get_total_cases_period_state(brasil_cases, "MS", min_date_p1.7, max_date_p1.7)
ms_casos_genomas <- get_prop_cases_target_lineage(gisaid_ms_lineages, ms_casos)

#mg
gisaid_mg_lineages <- get_subclades_data(gisaid_mg, min_date_p1.7,max_date_p1.7, target)
gisaid_mg_lineages <- format_data(gisaid_mg_lineages)
mg_casos <- get_total_cases_period_state(brasil_cases, "MG", min_date_p1.7, max_date_p1.7)
mg_casos_genomas <- get_prop_cases_target_lineage(gisaid_mg_lineages, mg_casos)

#pa
gisaid_pa_lineages <- get_subclades_data(gisaid_pa, min_date_p1.7,max_date_p1.7, target)
gisaid_pa_lineages <- format_data(gisaid_pa_lineages)
pa_casos <- get_total_cases_period_state(brasil_cases, "PA", min_date_p1.7, max_date_p1.7)
pa_casos_genomas <- get_prop_cases_target_lineage(gisaid_pa_lineages, pa_casos)

#pb
gisaid_pb_lineages <- get_subclades_data(gisaid_pb, min_date_p1.7,max_date_p1.7, target)
gisaid_pb_lineages <- format_data(gisaid_pb_lineages)
pb_casos <- get_total_cases_period_state(brasil_cases, "PB", min_date_p1.7, max_date_p1.7)
pb_casos_genomas <- get_prop_cases_target_lineage(gisaid_pb_lineages, pb_casos)

#rj
gisaid_rj_lineages <- get_subclades_data(gisaid_rj, min_date_p1.7,max_date_p1.7, target)
gisaid_rj_lineages <- format_data(gisaid_rj_lineages)
rj_casos <- get_total_cases_period_state(brasil_cases, "RJ", min_date_p1.7, max_date_p1.7)
rj_casos_genomas <- get_prop_cases_target_lineage(gisaid_rj_lineages, rj_casos)

#sp
gisaid_sp_lineages <- get_subclades_data(gisaid_sp, min_date_p1.7,max_date_p1.7, target)
gisaid_sp_lineages <- format_data(gisaid_sp_lineages)
sp_casos <- get_total_cases_period_state(brasil_cases, "SP", min_date_p1.7, max_date_p1.7)
sp_casos_genomas <- get_prop_cases_target_lineage(gisaid_sp_lineages, sp_casos)

#rn
gisaid_rn_lineages <- get_subclades_data(gisaid_rn, min_date_p1.7,max_date_p1.7, target)
gisaid_rn_lineages <- format_data(gisaid_rn_lineages)
rn_casos <- get_total_cases_period_state(brasil_cases, "RN", min_date_p1.7, max_date_p1.7)
rn_casos_genomas <- get_prop_cases_target_lineage(gisaid_rn_lineages, rn_casos)

#pi
gisaid_pi_lineages <- get_subclades_data(gisaid_pi, min_date_p1.7,max_date_p1.7, target)
gisaid_pi_lineages <- format_data(gisaid_pi_lineages)
pi_casos <- get_total_cases_period_state(brasil_cases, "PI", min_date_p1.7, max_date_p1.7)
pi_casos_genomas <- get_prop_cases_target_lineage(gisaid_pi_lineages, pi_casos)

#pe
gisaid_pe_lineages <- get_subclades_data(gisaid_pe, min_date_p1.7,max_date_p1.7, target)
gisaid_pe_lineages <- format_data(gisaid_pe_lineages)
pe_casos <- get_total_cases_period_state(brasil_cases, "PE", min_date_p1.7, max_date_p1.7)
pe_casos_genomas <- get_prop_cases_target_lineage(gisaid_pe_lineages, pe_casos)

#se
gisaid_se_lineages <- get_subclades_data(gisaid_se, min_date_p1.7,max_date_p1.7, target)
gisaid_se_lineages <- format_data(gisaid_se_lineages)
se_casos <- get_total_cases_period_state(brasil_cases, "SE", min_date_p1.7, max_date_p1.7)
se_casos_genomas <- get_prop_cases_target_lineage(gisaid_se_lineages, se_casos)

#ro
gisaid_ro_lineages <- get_subclades_data(gisaid_ro, min_date_p1.7,max_date_p1.7, target)
gisaid_ro_lineages <- format_data(gisaid_ro_lineages)
ro_casos <- get_total_cases_period_state(brasil_cases, "RO", min_date_p1.7, max_date_p1.7)
ro_casos_genomas <- get_prop_cases_target_lineage(gisaid_ro_lineages, ro_casos)

#rr
gisaid_rr_lineages <- get_subclades_data(gisaid_rr, min_date_p1.7,max_date_p1.7, target)
gisaid_rr_lineages <- format_data(gisaid_rr_lineages)
rr_casos <- get_total_cases_period_state(brasil_cases, "RR", min_date_p1.7, max_date_p1.7)
rr_casos_genomas <- get_prop_cases_target_lineage(gisaid_rr_lineages, rr_casos)

#to
gisaid_to_lineages <- get_subclades_data(gisaid_to, min_date_p1.7,max_date_p1.7, target)
gisaid_to_lineages <- format_data(gisaid_to_lineages)
to_casos <- get_total_cases_period_state(brasil_cases, "TO", min_date_p1.7, max_date_p1.7)
to_casos_genomas <- get_prop_cases_target_lineage(gisaid_to_lineages, to_casos)

cases_by_state <- get_genomes_cases_proportion(ac_casos_genomas,al_casos_genomas,ap_casos_genomas,
                                               am_casos_genomas,ba_casos_genomas,ce_casos_genomas,
                                               df_casos_genomas,es_casos_genomas,go_casos_genomas,
                                               ma_casos_genomas,mt_casos_genomas,ms_casos_genomas,
                                               mg_casos_genomas,pa_casos_genomas,pb_casos_genomas,
                                               pr_casos_genomas,pe_casos_genomas,pi_casos_genomas,
                                               rj_casos_genomas,rn_casos_genomas,rs_casos_genomas,
                                               ro_casos_genomas,rr_casos_genomas,sc_casos_genomas,
                                               sp_casos_genomas,se_casos_genomas,to_casos_genomas)
cases_by_state
p1.7_data <- left_join(cases_by_state, ibge_pop, by = "states_initials")
p1.7_data <- setNames(p1.7_data, c("state", "cases_prop_p.1.7", "cases_prop_others","genomes_p.1.7","genomes_total","pop_2022"))
p1.7_data$cases_total <- p1.7_data$cases_prop_p.1.7 + p1.7_data$cases_prop_others
write.table(p1.7_data, file = "p1_7_traits.tsv", sep = "\t", row.names = FALSE)




### UFSM PLOTS
sm_palmeiras <- read.csv("/Users/drt64727/pessoal_2/sars_rs/data/seq_info.tsv", header=TRUE, sep="\t")

sm_palmeiras$collection_date <- as.Date(sm_palmeiras$collection_date )
sm_palmeiras$collection_date <- as.character(format(sm_palmeiras$collection_date, "%Y-%m"))

## filters
# remove genomes without collection date, without location, and non human genomes
sm_palmeiras  <- sm_palmeiras %>% 
  filter(!is.na(collection_date))

sm_palmeiras <- sm_palmeiras %>%
  mutate(lineage = case_when(
    grepl("P.1", lineage) ~ "Gamma (P.1 e P.1*)",
    grepl("B.1.617.2", lineage) ~ "Delta (B.1.617.2 e AY.*)",
    grepl("AY.", lineage) ~ "Delta (B.1.617.2 e AY.*)",
    grepl("BA", lineage) ~ "Omicron (BA.*, BE.*, BQ.*)",
    grepl("BE.", lineage) ~ "Omicron (BA.*, BE.*, BQ.*)",
    grepl("BQ.", lineage) ~ "Omicron (BA.*, BE.*, BQ.*)",
    TRUE ~ lineage
  ))

sm_palmeiras$lineage[!(sm_palmeiras$lineage %in% target_lineages)] <- "Others"


br_gisaid <- br_gisaid %>%
  mutate(Pango.lineage = case_when(
    grepl("P.1", Pango.lineage) ~ "Gamma (P.1 e P.1*)",
    grepl("B.1.617.2", Pango.lineage) ~ "Delta (B.1.617.2 e AY.*)",
    grepl("AY.", Pango.lineage) ~ "Delta (B.1.617.2 e AY.*)",
    grepl("BA", Pango.lineage) ~ "Omicron (BA.*, BE.*, BQ.*)",
    grepl("BE.", Pango.lineage) ~ "Omicron (BA.*, BE.*, BQ.*)",
    grepl("BQ.", Pango.lineage) ~ "Omicron (BA.*, BE.*, BQ.*)",
    TRUE ~ Pango.lineage
  ))

br_gisaid$Pango.lineage[!(br_gisaid$Pango.lineage %in% target_lineages)] <- "Others"

getwd()
## plot ufsm-palmeiras
sm_palmeiras$collection_date <- factor(sm_palmeiras$collection_date, levels = x_order_month)

sm_palmeiras_summary <- sm_palmeiras %>%
  count(collection_date, lineage) %>%
  group_by(collection_date) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

sm_palmeiras_summary$collection_date <- factor(sm_palmeiras_summary$collection_date, levels = x_order_month)

target_lineages_colors <-  c("Others" = "#8B8B8B", 
                             "Gamma (P.1 e P.1*)" = "#FF6666",
                             "Delta (B.1.617.2 e AY.*)" = "#99FF33",
                             "Omicron (BA.*, BE.*, BQ.*)" = "#99CCFF",
                             "P.2" = "#FFB266",
                             "B.1.1.28" = "#CCFFFF",
                             "B.1.1.33" = "#FF99FF",
                             "B.1.1" = "#FFFF99")


sm_palmeiras_summary <- sm_palmeiras_summary %>%
  group_by(collection_date) %>%
  arrange(-freq) %>% 
  mutate(total_freq = sum(freq),
         ymin = cumsum(freq) - freq,
         ymax = cumsum(freq)) %>%
  ungroup()

ggplot(sm_palmeiras_summary, aes(x=factor(collection_date, levels = x_order_month), y = freq, fill = lineage)) +
  geom_rect(aes(ymin = ymin, ymax = ymax, xmin = as.numeric(factor(collection_date, levels = x_order_month))-0.4, xmax = as.numeric(factor(collection_date, levels = x_order_month))+0.4), color="black") +
  geom_line(aes(y = n/200, group = 1), color = "black", size = 2, alpha = 0.5) +
  labs(x = "Período (ano-mês)", y = "Frequência") +
  scale_y_continuous(sec.axis = sec_axis(~.*200, name = "Genomas Sequenciados (n)")) +
  scale_x_discrete(limits = x_order_month) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.line = element_line(size = 0.5, colour = "black"),
        axis.title = element_text(size = 13),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent"),
        legend.box.just = "left",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.spacing = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = target_lineages_colors)


gisaid_rs_summary <- gisaid_rs %>%
  count(Collection.date, Pango.lineage) %>%
  group_by(Collection.date) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

gisaid_rs_summary <- gisaid_rs_summary %>%
  group_by(Collection.date) %>%
  arrange(-freq) %>% 
  mutate(total_freq = sum(freq),
         ymin = cumsum(freq) - freq,
         ymax = cumsum(freq)) %>%
  ungroup()


ggplot(gisaid_rs_summary, aes(x=factor(Collection.date, levels = x_order_month), y = freq, fill = Pango.lineage)) +
  geom_rect(aes(ymin = ymin, ymax = ymax, xmin = as.numeric(factor(Collection.date, levels = x_order_month))-0.4, xmax = as.numeric(factor(Collection.date, levels = x_order_month))+0.4), color="black") +
  geom_line(aes(y = n/650, group = 1), color = "black", size = 2, alpha = 0.5) +
  labs(x = "Período (ano-mês)", y = "Frequência") +
  scale_y_continuous(sec.axis = sec_axis(~.*650, name = "Genomas Sequenciados (n)")) +
  scale_x_discrete(limits = x_order_month) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.line = element_line(size = 0.5, colour = "black"),
        axis.title = element_text(size = 13),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent"),
        legend.box.just = "left",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.spacing = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = target_lineages_colors)


br_gisaid_summary <- br_gisaid %>%
  count(Collection.date, Pango.lineage) %>%
  group_by(Collection.date) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

br_gisaid_summary <- br_gisaid_summary %>%
  group_by(Collection.date) %>%
  arrange(-freq) %>% 
  mutate(total_freq = sum(freq),
         ymin = cumsum(freq) - freq,
         ymax = cumsum(freq)) %>%
  ungroup()


ggplot(br_gisaid_summary, aes(x=factor(Collection.date, levels = x_order_month), y = freq, fill = Pango.lineage)) +
  geom_rect(aes(ymin = ymin, ymax = ymax, xmin = as.numeric(factor(Collection.date, levels = x_order_month))-0.4, xmax = as.numeric(factor(Collection.date, levels = x_order_month))+0.4), color="black") +
  geom_line(aes(y = n/25000, group = 1), color = "black", size = 2, alpha = 0.5) +
  labs(x = "Período (ano-mês)", y = "Frequência") +
  scale_y_continuous(sec.axis = sec_axis(~.*25500, name = "Genomas Sequenciados (n)")) +
  scale_x_discrete(limits = x_order_month) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.line = element_line(size = 0.5, colour = "black"),
        axis.title = element_text(size = 13),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "transparent"),
        legend.box.just = "left",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.spacing = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        plot.title = element_text(size = 14, face = "bold")) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = target_lineages_colors)
