load_data_bam_vs_cram <- function(input_folder, input_file, repetition_name, format) {
  
  location <- paste0(input_folder, input_file, "")
  storage <- read.csv(file = Sys.glob(file.path(location, '*tsv')), sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  trace <- read.csv(file = Sys.glob(file.path(location, 'execution*')), sep = "\t", stringsAsFactors = FALSE)
  
  names(storage)[names(storage) == "V2"] <- "workdir"
  names(storage)[names(storage) == "V1"] <- "workdir_bytes"
  trace <- trace %>% mutate('repetition_name' = repetition_name)
  trace <- trace %>% mutate('format' = format)
  
  merged <- merge(trace, storage, by = "workdir")
  
  return(merged)
}

format_aws_costs <- function(cost_table, date) {
  
  # transpose data frame & use first column from original dataframe as column headers in transposed dataframe
  cost_table_formatted <- setNames(data.frame(t(cost_table[,-1])), cost_table[,1])
  
  # Rename the `date` column to costs, since it carries the costs for that date for better handling
  names(cost_table_formatted)[2] <- "costs"

  # create a new column named 'type' that has the current row names (needed for ggplot2)
  cost_table_formatted <- tibble::rownames_to_column(cost_table_formatted, "type")

  # Round too last 3 digits & remove costs below 0.00$
  cost_table_formatted <- cost_table_formatted %>% mutate_if(is.numeric, round, digits=3)
  cost_table_formatted <- filter(cost_table_formatted, costs >= 0.001 , .preserve = FALSE)
  return (cost_table_formatted)
}

bytesto <- function(bytes, to, bsize=1024){
  a <- data.frame ( unit = c('k', 'm', 'g', 't', 'p', 'e' ), factor = c(1,2,3,4,5,6), stringsAsFactors=FALSE)

  get_factor <- a$factor
  names(get_factor) <- a$unit
  
  return (bytes / (bsize ^ get_factor[to]))
}

format_bam_vs_cram <- function(df, format){
  
  df <- df %>%   
    # Get last element (by :) from process name
    mutate('simple_name' = sub("^.*:", "", process)) %>%
    # Remove unrelevant processes
    filter(simple_name != 'ASSESS_SIGNIFICANCE') %>%
    filter(simple_name != 'BCFTOOLS_STATS') %>%
    filter(simple_name != 'BCFTOOLS_SORT') %>%
    filter(simple_name != 'CALCULATECONTAMINATION') %>%
    filter(simple_name != 'CAT_MPILEUP') %>%
    filter(simple_name != 'CNNSCOREVARIANTS') %>%
    filter(simple_name != 'CREATE_INTERVALS_BED') %>%
    filter(simple_name != 'BWAMEM1_MEM') %>%
    filter(simple_name != 'FASTP') %>%
    filter(simple_name != 'FASTQC') %>%
    filter(simple_name != 'FILTERMUTECTCALLS') %>%
    filter(simple_name != 'FILTERVARIANTTRANCHES') %>%
    filter(simple_name != 'FREEC2CIRCOS') %>%
    filter(simple_name != 'FREEC2BED') %>%
    filter(simple_name != 'GATHERPILEUPSUMMARIES_NORMAL') %>%
    filter(simple_name != 'GATHERPILEUPSUMMARIES_TUMOR') %>%
    filter(simple_name != 'GATHERBQSRREPORTS') %>%
    filter(simple_name != 'SAMTOOLS_BAMTOCRAM') %>%
    #filter(simple_name != 'GATK4_MARKDUPLICATES') %>%
    filter(simple_name != 'LEARNREADORIENTATIONMODEL') %>%
    filter(simple_name != 'MAKEGRAPH') %>%
    filter(simple_name != 'MERGE_DEEPVARIANT_GVCF') %>%
    filter(simple_name != 'MERGE_DEEPVARIANT_VCF') %>%
    filter(simple_name != 'MERGE_HAPLOTYPECALLER') %>%
    filter(simple_name != 'MERGE_FREEBAYES') %>%
    filter(simple_name != 'MERGE_MANTA_DIPLOID') %>%
    filter(simple_name != 'MERGE_MANTA_SMALL_INDELS') %>%
    filter(simple_name != 'MERGE_MANTA_SOMATIC') %>%
    filter(simple_name != 'MERGE_MANTA_SV') %>%
    filter(simple_name != 'MERGE_MUTECT2') %>%
    filter(simple_name != 'MERGE_STRELKA') %>%
    filter(simple_name != 'MERGE_STRELKA_GENOME') %>%
    filter(simple_name != 'MERGE_STRELKA_INDELS') %>%
    filter(simple_name != 'MERGE_STRELKA_SNVS') %>%
    filter(simple_name != 'MERGEMUTECTSTATS') %>%
    filter(simple_name != 'SVDB_MERGE') %>%
    filter(simple_name != 'TABIX_BGZIP_TIDDIT_SV') %>%
    filter(simple_name != 'VCFTOOLS_SUMMARY') %>%
    filter(simple_name != 'VCFTOOLS_TSTV_COUNT') %>%
    filter(simple_name != 'VCFTOOLS_TSTV_QUAL') %>%
    filter(simple_name != 'UNZIP_RT') %>%
    # Rename process to something less convoluted
    mutate(simple_name = recode(simple_name, 'SAMTOOLS_STATS_CRAM' = 'SAMTOOLS_STATS')) %>% 
    mutate(simple_name = recode(simple_name, 'GATK4_MARKDUPLICATES' = 'MARKDUPLICATES')) %>%
    mutate(simple_name = recode(simple_name, 'GATK4_BASERECALIBRATOR' = 'BASERECALIBRATOR')) %>%
    mutate(simple_name = recode(simple_name, 'GATK4_APPLYBQSR' = 'APPLYBQSR')) %>%
    mutate(simple_name = recode(simple_name, 'GATK4_GATHERBQSRREPORTS' = 'GATHERBQSRREPORTS')) %>%
    mutate(simple_name = recode(simple_name, 'INDEX_CRAM' = 'INDEX_RECALIBRATED')) %>%
    mutate(simple_name = recode(simple_name, 'MERGE_CRAM' = 'MERGE_RECALIBRATED')) %>%
    mutate(simple_name = recode(simple_name, 'GATK4_HAPLOTYPECALLER' = 'HAPLOTYPECALLER')) %>%
    mutate(simple_name = recode(simple_name, 'GETPILEUPSUMMARIES_NORMAL' = 'GETPILEUPSUMMARIES')) %>%
    mutate(simple_name = recode(simple_name, 'GETPILEUPSUMMARIES_TUMOR' = 'GETPILEUPSUMMARIES'))
  
  
  return (df %>% 
    # Remove rows with column value '-' and cast columns to numeric
    filter(vmem !='-') %>%
    #filter(vmem !='0') %>%
    mutate_all(type.convert, as.is=TRUE) %>%
    # Convert vmem to megabytes and gigabytes
    mutate('vmem_MB' = bytesto(vmem, 'm')) %>%
    mutate('vmem_GB' = bytesto(vmem, 'g')) %>%
    # Convert peak_vmem to megabytes and gigabytes
    mutate('peak_vmem_MB' = bytesto(peak_vmem, 'm')) %>%
    mutate('peak_vmem_GB' = bytesto(peak_vmem, 'g')) %>%
    # Convert requested memory to megabytes and gigabytes
    mutate('memory_MB' = bytesto(memory, 'm')) %>%
    mutate('memory_GB' = bytesto(memory, 'g')) %>%
    # Get an actual CPU number by dividing the used cpu number by 100
    mutate('num_cpu' = X.cpu/100) %>%
    # Convert workdir storage to megabytes and gigabytes
    mutate('workdir_MB' = bytesto(workdir_bytes, 'm')) %>%
    mutate('workdir_GB' = bytesto(workdir_bytes, 'g')) %>%
    mutate('workdir_TB' = bytesto(workdir_bytes, 't')) %>%
    # Convert seconds to minutes
    mutate('realtime_min' = realtime/(60*1000)) %>%
    # Get tumor/normal
    mutate('status' = case_when(
      grepl("vs", tag) ~ "paired",
      grepl("N", tag) ~ "normal", 
      grepl("T", tag) ~ "tumor",
      TRUE ~ "other")) %>%
    filter(status !='other') %>%
    mutate('sample' = case_when(
      grepl("CHC892", tag) ~ "CHC892", 
      grepl("CHC912", tag) ~ "CHC912", 
      grepl("CHC2111", tag) ~ "CHC2111", 
      grepl("CHC2113", tag) ~ "CHC2113",
      grepl("CHC2115", tag) ~ "CHC2115",
      grepl("HCC1395", tag) ~ "HCC1395",
      TRUE ~ "other")) %>%
    # Classify processes into subcategories
    mutate('type' = case_when(
      grepl("APPLYBQSR", simple_name) ~ "Pre-processing", 
      grepl("ASCAT", simple_name) ~ "CNV",
      grepl("BASERECALIBRATOR", simple_name) ~ "Pre-processing", 
      grepl("CALCULATECONTAMINATION", simple_name) ~ "SNP",
      grepl("CNVKIT_BATCH", simple_name) ~ "CNV",
      grepl("FREEBAYES", simple_name) ~ "SNP",
      grepl("GATHERBQSRREPORT", simple_name) ~ "Pre-processing", 
      grepl("MARKDUPLICATES", simple_name) ~ "Pre-processing", 
      grepl("GETPILEUPSUMMARIES_NORMAL", simple_name) ~ "SNP",
      grepl("GETPILEUPSUMMARIES_TUMOR", simple_name) ~ "SNP",
      grepl("HAPLOTYPECALLER", simple_name) ~ "SNP",
      grepl("INDEX_CRAM", simple_name) ~ "Pre-processing", 
      grepl("MANTA_GERMLINE", simple_name) ~ "SV",
      grepl("MANTA_SOMATIC", simple_name) ~ "SV",
      grepl("MERGE_CRAM", simple_name) ~ "Pre-processing", 
      grepl("MANTA_DIPLOID", simple_name) ~ "SV",
      grepl("MOSDEPTH", simple_name) ~ "Pre-processing", 
      grepl("MSISENSORPRO_MSI_SOMATIC", simple_name) ~ "SNP",
      grepl("SAMTOOLS_BAMTOCRAM", simple_name) ~ "Pre-processing", 
      grepl("SAMTOOLS_STATS", simple_name) ~ "Pre-processing", 
      grepl("SAMTOOLS_STATS_CRAM", simple_name) ~ "Pre-processing", 
      grepl("STRELKA_SINGLE", simple_name) ~ "SNP",
      grepl("STRELKA_SOMATIC", simple_name) ~ "SNP",
      grepl("MUTECT2", simple_name) ~ "SNP",
      grepl("TIDDIT_SV", simple_name) ~ "SV",
      TRUE ~ "other"))
  )
}

bam_vs_cram_categories <- list(
  "CNV" = "CNV",
  "Pre-processing" = "Preprocessing",
  "SNP" = "SNP/Indels",
  "SV"  = "SV"
  
)

resource_names <- list(
  "peak_vmem_GB" = "Mem(GB)",
  "num_cpu" = "CPUs",
  "workdir_GB" = "storage(GB)",
  "realtime_min" = "Time(min)"
)

variable_labeller_y <- function(variable,value){
  if (variable=='measure') {
    return(resource_names[value])
  } else {
    return(bam_vs_cram_categories)
  }
}


plot_bam_vs_cram_boxplot <- function(df, xaxis, yaxis, boxplot_color, outputname, results_folder, yaxisname) {
  my_plot <- 
    ggboxplot(df, x=xaxis, y=yaxis, fill = boxplot_color) + 
    # Increase point size in the legend
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) +
    theme(axis.text.x = element_text(angle=90, hjust=1, size=10),
          axis.text.y = element_text(size=10),
          axis.title.x = element_blank(),
          legend.position="top",  legend.box="vertical", legend.margin=margin()) +
    labs(y = yaxisname) #+
    #stat_compare_means(aes_string(group = boxplot_color), label = "p.signif")
  
  ### plots each run individually and adds dots per plot
  # my_plot <- 
  #   ggboxplot(df, x=xaxis, y=yaxis, fill = boxplot_color, outlier.colour="black", outlier.shape=16, outlier.size=.5, notch=FALSE, lwd =.1, width = 0.5) +
  #   # Add dots for each (sub)sample
  #   geom_point(shape=16, aes_string(fill=boxplot_color, colour = jitter_color), size = 0.1, alpha = 1, position = position_jitterdodge(jitter.width = 0.05)) +
  #   # Plot requested resources as black *
  #   { if(!is.null(reference)) geom_point(shape=8, aes_string(fill=boxplot_color, y=reference), position = position_jitterdodge(jitter.width = 0)) }+
  #   # Increase point size in the legend
  #   guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) +
  #   theme(axis.text.x = element_text(angle=90, hjust=1, size=10),
  #        axis.text.y = element_text(size=10),
  #        axis.title.x = element_blank(),
  #        legend.position="top",  legend.box="vertical", legend.margin=margin()) +
  #   labs(y = yaxisname)
  
  # ggsave(plot=my_plot, filename = paste0(results_folder,outputname, ".png"), device="png",
  #        width = 20, height = 10, units="cm")
  # ggsave(plot=my_plot, filename = paste0(results_folder,outputname, ".pdf"), device="pdf", 
  #        width=20, height=10, units="cm")
  return(my_plot)
}

plot_cumulative_storage <- function(merged, outputfile) {
  merged_cumulative_storage <- merged %>%
    filter(repetition_name != 'cram_2') %>%
    filter(repetition_name != 'cram_3') %>%
    filter(repetition_name != 'bam_2') %>%
    filter(repetition_name != 'bam_3') %>%
    group_by(simple_name, format) %>%
    summarise(storage = sum(workdir_GB), .groups = 'drop')
  
  
  storage_cumulative <-
    ggbarplot(merged_cumulative_storage,
              x = "simple_name", y = "storage", fill = "format",
              position = position_dodge(0.8), orientation = "horiz") +
    geom_text(aes(label = round(storage,3),  group = format), colour = "black", size = 3,
              position = position_dodge(.9),
              hjust = 1.2) +
    rremove("legend.title") +
    theme(axis.text.x = element_text(angle=45, hjust = 1, size=15),
          axis.text.y = element_text(size=15),legend.position="top",
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=15),
          legend.text = element_text(size = 15)) +
    labs(y = "work dir (GB)") +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma=0.001),
      breaks = c(0, 10^(-2:6)))#+
    #stat_compare_means(aes(group = format), label = "p.signif")
  
  
  ggsave(plot=storage_cumulative, filename = paste0(results_folder,outputfile), device="png", dpi = 600)
  
  return(storage_cumulative)
}

plot_summary <- function(merged, file_name) {
  summary_plot <- ggboxplot(merged, x="simple_name", 
                               y="value", fill = "format", outlier.colour="black", outlier.shape=16, 
                               outlier.size=.5, notch=FALSE, lwd =.1, width = 0.5) + 
    rremove("legend.title") +
    #geom_point(shape=16, aes_string(fill="format", colour = "sample"), 
    #           size = 0.1, alpha = 1, position = position_jitterdodge(jitter.width = 0.05)) + 
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) + 
    facet_grid(measure~., scales="free", space="free_x", 
               labeller=variable_labeller_y, switch="y") +
    theme(axis.text.x = element_text(angle=45, hjust = 1, size=10), 
          axis.text.y = element_text(size=10),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          strip.text.y = element_text(size = 10) ,
          legend.position="top", legend.box="vertical", legend.margin=margin(),
          legend.text = element_text(size = 10)) #+ 
  # stat_compare_means(aes(group = format), label = "p.signif")
  #
  
  ggsave(plot=summary_plot, filename = paste0(results_folder, file_name), device="png", 
         dpi = 600)
  
  return(summary_plot)
}

format_splitting <- function(df){
  return (df %>% 
            # Remove rows with column value '-' and cast columns to numeric
            #filter(vmem != '-') %>%
            #filter(vmem != '0') %>%
            # Only keep relevant processes
            #filter(grepl('BWAMEM1_MEM|FASTP|GATK4_MARKDUPLICATES', process))
            filter(grepl('BWAMEM1_MEM|FASTP|GATK4_MARKDUPLICATES|APPLYBQSR|STRELKA|MUTECT2', process)) %>%
            mutate_all(type.convert, as.is=TRUE) %>%
            # Convert vmem to megabytes and gigabytes
            #mutate('vmem_MB' = bytesto(vmem, 'm')) %>%
            #mutate('vmem_GB' = bytesto(vmem, 'g')) %>%
            # Convert peak_vmem to megabytes and gigabytes
            #mutate('peak_vmem_MB' = bytesto(peak_vmem, 'm')) %>%
            #mutate('peak_vmem_GB' = bytesto(peak_vmem, 'g')) %>%
            # Convert requested memory to megabytes and gigabytes
            #mutate('memory_MB' = bytesto(memory, 'm')) %>%
            #mutate('memory_GB' = bytesto(memory, 'g')) %>%
            # Get an actual CPU number by dividing the used cpu number by 100
            #mutate('num_cpu' = X.cpu/100) %>%
            # Get last element (by :) from process name
            mutate('simple_name' = sub("^.*:", "", process)) %>%
            # Convert workdir storage to megabytes and gigabytes
            mutate('workdir_MB' = bytesto(workdir_bytes, 'm')) %>%
            mutate('workdir_GB' = bytesto(workdir_bytes, 'g')) %>%
            # Convert seconds to minutes
            mutate('realtime_min' = realtime/(60*1000)) %>%
            # Get tumor/normal
            mutate('status' = case_when(
              grepl("vs", tag) ~ "paired",
              grepl("N", tag) ~ "sample1 (12.5GB)", 
              grepl("T", tag) ~ "sample2 (14.8GB)",
              TRUE ~ "other")) %>%
            filter(status !='other') %>%
            mutate('sample' = case_when(
              grepl("CHC892", tag) ~ "CHC892", 
              grepl("CHC912", tag) ~ "CHC912", 
              grepl("CHC2111", tag) ~ "CHC2111", 
              grepl("CHC2113", tag) ~ "CHC2113",
              grepl("HCC1395", tag) ~ "HCC1395",
              TRUE ~ "other"))
 
  )
}