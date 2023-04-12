#####################
## BAM vs CRAM ######
#####################
resource_names <- list("max_peak_vmem_GB" = "Mem(GB)","max_num_cpu" = "CPUs","avg_realtime_min" = "Time(min)")

load_data_bam_vs_cram <- function(input_folder, input_file, repetition_name, format) {
  
  location <- paste0(input_folder, input_file, "")
  storage <- read.csv(file = Sys.glob(file.path(location, '*tsv')), sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  trace <- read.csv(file = Sys.glob(file.path(location, 'execution*')), sep = "\t", stringsAsFactors = FALSE)
  
  storage <- rename_storage_columns(storage)
  trace <- trace %>% mutate('repetition_name' = repetition_name)
  trace <- trace %>% mutate('format' = format)
  
  merged <- merge(trace, storage, by = "workdir")
  
  return(merged)
}

format_bam_vs_cram <- function(df, format){
  
  df_filter <- df %>%   
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
    filter(simple_name != 'FREEC_SOMATIC') %>%
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
    mutate(simple_name = recode(simple_name, 'GATK4_MARKDUPLICATES_SPARK' = 'MARKDUPLICATES_SPARK')) %>%
    mutate(simple_name = recode(simple_name, 'GATK4_ESTIMATELIBRARYCOMPLEXITY' = 'ESTIMATELIBRARYCOMPLEXITY')) %>%
    mutate(simple_name = recode(simple_name, 'GATK4_BASERECALIBRATOR' = 'BASERECALIBRATOR')) %>%
    mutate(simple_name = recode(simple_name, 'GATK4_APPLYBQSR' = 'APPLYBQSR')) %>%
    mutate(simple_name = recode(simple_name, 'GATK4_GATHERBQSRREPORTS' = 'GATHERBQSRREPORTS')) %>%
    mutate(simple_name = recode(simple_name, 'INDEX_CRAM' = 'INDEX_RECALIBRATED')) %>%
    mutate(simple_name = recode(simple_name, 'MERGE_CRAM' = 'MERGE_RECALIBRATED')) %>%
    mutate(simple_name = recode(simple_name, 'GATK4_HAPLOTYPECALLER' = 'HAPLOTYPECALLER')) %>%
    mutate(simple_name = recode(simple_name, 'GETPILEUPSUMMARIES_NORMAL' = 'GETPILEUPSUMMARIES')) %>%
    mutate(simple_name = recode(simple_name, 'GETPILEUPSUMMARIES_TUMOR' = 'GETPILEUPSUMMARIES')) %>%
    mutate(simple_name = recode(simple_name, 'STRELKA_SINGLE' = 'STRELKA')) %>%
    mutate(simple_name = recode(simple_name, 'STRELKA_SOMATIC' = 'STRELKA')) %>%
    mutate(simple_name = recode(simple_name, 'MANTA_GERMLINE' = 'MANTA')) %>%
    mutate(simple_name = recode(simple_name, 'MANTA_SOMATIC' = 'MANTA'))
    
  return (df_filter %>% 
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
           TRUE ~ "other")) 
    )
  
}

collect_compute_values <- function(df) {
  return (df %>%     
          # Of processes that are split across many sub-samples, i.e. ApplyBQSR only keep the one with the max time per sample & repetition -> 30 data points for each process
          group_by(tag, repetition_name, simple_name) %>%
          filter(realtime_min == max(realtime_min)) %>%
          ungroup() %>%
          
          ## Use average time, max cpus, and max memory over the 3 repetitions per sample
          group_by(tag, format, simple_name) %>%
          mutate(avg_realtime_min = mean(realtime_min)) %>%
          mutate(max_num_cpu = max(num_cpu)) %>%
          mutate(max_vmem_GB = max(vmem_GB)) %>%
          mutate(max_peak_vmem_GB = max(peak_vmem_GB)) %>%
          ungroup() %>%
          
          ## Collapse repetitions -> 10 data points left per preprocessing process
          select(simple_name, sample, tag, format, cpus, memory_GB, avg_realtime_min,
                 max_num_cpu, max_peak_vmem_GB, max_vmem_GB) %>%
          group_by(simple_name, tag, format) %>%
          distinct(avg_realtime_min, max_num_cpu, max_peak_vmem_GB, max_vmem_GB,  .keep_all = TRUE))
}

plot_bam_vs_cram_boxplot <- function(df, xaxis, yaxis, boxplot_color, outputname, results_folder, yaxisname, reference) {
  
  my_plot <- 
    ggboxplot(df, x=xaxis, y=yaxis, width = 0.5, outlier.size=.5, color = boxplot_color) + 
    # Plot requested resources as black triangles
    { if(!is.null(reference)) geom_point(shape=6, mapping = aes(!!!ensyms(fill=boxplot_color, y=reference)), position = position_jitterdodge(jitter.width = 0)) }+
    
    # Increase point size in the legend
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) +
    theme(axis.text.x = element_text(angle=45, hjust=1, size=9),
          axis.text.y = element_text(size=9),
          axis.title.y = element_text(size=10),
          axis.title.x = element_blank(),
          legend.position="top",  legend.box="vertical", legend.margin=margin()) +
    labs(y = yaxisname) + 
    # Add sign marker for paired wolcoxon test
    stat_compare_means(mapping = aes(!!!ensyms(group = boxplot_color)), label = "p.signif",  method = "wilcox.test", paired = TRUE, hide.ns = TRUE, vjust = 1)
  
  if(!is.null(outputname)) {
    ggsave(plot=my_plot, filename = paste0(results_folder,outputname, ".png"), device="png",
          width = 20, height = 10, units="cm")
    ggsave(plot=my_plot, filename = paste0(results_folder,outputname, ".pdf"), device="pdf", 
          width=20, height=10, units="cm")
  }
  
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
              x = "simple_name", y = "storage", fill = "format", color = "format",
              position = position_dodge(0.8), orientation = "horiz") +
    geom_text(aes(label = round(storage,3),  group = format, fontface = "bold"), colour = "white", size = 5,
              position = position_dodge(.85),
              hjust = 1.1) +
    rremove("legend.title") +
    theme(axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15), legend.position="top",
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=15),
          legend.text = element_text(size = 15)) +
    labs(y = "work dir (GB)") +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma=0.001),
      breaks = c(0, 10^(-2:6))) +
    stat_compare_means(aes(group = format), label = "p.signif", vjust = 1, hide.ns = TRUE, paired=TRUE)

  ggsave(plot=storage_cumulative, filename = paste0(results_folder,outputfile), device="png", width =10, height = 10)
  
  return(storage_cumulative)
}

plot_summary <- function(df, file_name, results_folder) {
  cpu <- plot_bam_vs_cram_boxplot(df=df,
                                  xaxis="simple_name",
                                  yaxis="max_num_cpu",
                                  boxplot_color="format",
                                  yaxisname = "CPUs",
                                  reference = NULL,
                                  outputname = NULL,
                                  results_folder = NULL) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    rremove("legend.title")
  
  mem <- plot_bam_vs_cram_boxplot(df=df, 
                                  xaxis="simple_name", 
                                  yaxis="max_peak_vmem_GB", 
                                  boxplot_color="format",
                                  yaxisname = "Memory (GB)",
                                  reference = NULL,
                                  outputname = NULL,
                                  results_folder = NULL)  + 
    theme(axis.text.x = element_text(angle=20, size=8)) +
    rremove("legend.title")
  
  time <- plot_bam_vs_cram_boxplot(df=df, 
                                   xaxis="simple_name", 
                                   yaxis="avg_realtime_min", 
                                   boxplot_color="format",
                                   yaxisname = "Time (min)",
                                   reference = NULL,
                                   outputname = NULL,
                                   results_folder = NULL) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    rremove("legend.title")
  
  summary_plot <-ggpubr::ggarrange(time, cpu, mem, ncol=1, nrow=3, common.legend = TRUE, align = "v", heights = c(0.7,0.7,1)) 
  ggsave(plot=summary_plot, filename = paste0(results_folder, file_name), device="png", dpi = 600)

  return(summary_plot)
}

#####################
##### Dataflow ######
#####################
format_splitting <- function(df){
  return (df %>% 
            # Remove rows with status FAILED
            filter(status != 'FAILED') %>%
            # Only keep relevant processes
            filter(grepl('BWAMEM1_MEM|FASTP|GATK4_MARKDUPLICATES|GATK4_BASERECALIBRATOR|GATK4_GATHERBQSRREPORTS|GATK4_APPLYBQSR|MERGE_CRAM|DEEPVARIANT|HAPLOTYPECALLER|STRELKA|MUTECT2|FREEBAYES', process)) %>%
            mutate_all(type.convert, as.is=TRUE) %>%
            # Get an actual CPU number by dividing the used cpu number by 100
            mutate('num_cpu' = X.cpu/100) %>%
            # Get last element (by :) from process name
            mutate('simple_name' = sub("^.*:", "", process)) %>%
            # Convert workdir storage to megabytes and gigabytes
            mutate('workdir_MB' = bytesto(workdir_bytes, 'm')) %>%
            mutate('workdir_GB' = bytesto(workdir_bytes, 'g')) %>%
            # Convert milli-seconds to minutes
            mutate('realtime_min' = realtime/(60*1000)) %>%
            mutate('realtime_h' = realtime/(3600*1000)) %>%
            mutate('cpu_h' = realtime_h * num_cpu) %>%
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
              TRUE ~ "other"))
  )
}

rename_storage_columns <- function(df) {
  names(df)[names(df) == "V2"] <- "workdir"
  names(df)[names(df) == "V1"] <- "workdir_bytes"
  return (df)
}

format_process_mapping <- function(df, process_name){
  df_process <- df %>% filter(grepl(process_name, simple_name))
  
  df_time <- df_process %>% group_by(fastp, sample, status, sample_status) %>% summarise(
    across(realtime_min, max, .names = 'max_sample_realtime'),
    across(cpu_h, sum, .names = 'sum_sample_cpuh'))
  
  df_storage <- df_process %>% group_by(fastp, sample, status, sample_status) %>% summarise(sum = sum(workdir_GB))
  df_storage$fastp <- as.factor(df_storage$fastp)
  
  return(list(process = df_process, time= df_time, storage=df_storage))
}

format_process_splitting <- function(df, process_name){
  df_process <- df %>% filter(grepl(process_name, simple_name))
  
  df_time <- df_process %>% group_by(intervals, sample, status, sample_status) %>% summarise(across(realtime_min, max, .names = 'max_sample_realtime'),
                                                                                             across(cpu_h, sum, .names = 'sum_sample_cpuh'))
  
  df_storage <- df_process %>% group_by(intervals, sample, status, sample_status) %>% summarise(sum = sum(workdir_GB))
  df_storage$intervals <- as.factor(df_storage$intervals)
  
  return(list(process = df_process, time= df_time, storage=df_storage))
}

x_axis_dataflow <- function(type) {
  if(identical('fastp', type)) {
    scale_x_discrete(breaks = c(1, 4, 8, 12, 16), labels = c(1, 4, 8, 12, 16)) }
  else {
    scale_x_discrete(breaks = c(1, 21, 40, 78, 124), labels = c(1, 21, 40, 78, 124))
  }
}

dataflow_theme = theme(axis.title.x = element_text(size=15),
              axis.title.y = element_text(size=15),
              plot.title = element_text(size = 20, hjust = 0.))

plot_dataflow_single_process <- function(df_max_time, df, df_storage, group, title, xaxis, outputname, results_folder){

  line <-  ggviolin(data=df_max_time, x=group, y="max_sample_realtime", draw_quantile = 0.5, color=group, width = 0.9) +
            x_axis_dataflow(group) +
            labs(y = "time (min)", x = xaxis) +
            dataflow_theme

  line_cpuh <- ggviolin(data=df_max_time, x=group, y="sum_sample_cpuh", draw_quantile = 0.5, color = group,  width = 0.9) +
                x_axis_dataflow(group) +
                labs(y = "CPUh", x = xaxis) +
                dataflow_theme

  bar <- ggviolin(df_storage, x=group, y="sum", draw_quantile = 0.5, color = group, width = 0.9) +
          x_axis_dataflow(group) +
          labs(y = "work dir (GB)", x = xaxis) +
          dataflow_theme +
          theme(legend.text = element_text(size = 8))

  plot <- ggpubr::ggarrange(line, bar, line_cpuh, legend= "none", ncol=3)
  ann_plot <- annotate_figure(plot, top = text_grob(title, face = "bold", size = 14))

  ggsave(plot=ann_plot, filename = paste0(results_folder,outputname, ".png"), device="png", width=20, height=20, units="cm")
  #ggsave(plot=ann_plot, filename = paste0(results_folder,outputname, ".pdf"), device="pdf",
  #       width=20, height=20, units="cm")

  return(ann_plot)
}

plot_splitting_summary <- function(df_time, group, xaxis, df_storage, title, outputname){
  
  line <- ggviolin(data=df_time, x=group, y="sum_combined_per_sample_realtime", draw_quantile = 0.5, color = group,  width = 0.9) +
          x_axis_dataflow(group) +
          labs(y = "time (min)", x = xaxis) +
          dataflow_theme
              
  line_cpuh <-  ggviolin(data=df_time, x=group, y="sum_combined_per_sample_cpuh", draw_quantile = 0.5, color = group, width = 0.9) +
                x_axis_dataflow(group) +
                labs(y = "CPUh", x = xaxis) +
                dataflow_theme
  
  bar <- ggviolin(data=df_storage, x=group, y="sum",  draw_quantile = 0.5, color = group, width = 0.9) +
            x_axis_dataflow(group) +
            labs(y = "work dir (GB)", x = xaxis) +
            dataflow_theme +
            theme(legend.text = element_text(size = 8)) +
            rremove("legend.title")
  
  plot <- ggpubr::ggarrange(line, bar, line_cpuh, legend = "none", ncol = 3)
  ann_plot <- annotate_figure(plot, top = text_grob(title, face = "bold", size = 14))
  
  ggsave(plot=ann_plot, filename = paste0(results_folder,outputname, ".png"), device="png",
         width=20, height=10, units="cm")
  #ggsave(plot=ann_plot, filename = paste0(results_folder,outputname, ".pdf"), device="pdf", 
  #      width=20, height=10, units="cm")
  
  return(ann_plot)
}

plot_variantcaller_summary <-  function(df_time, group, df_storage, title, outputname) {
  
  line <- facet(ggboxplot(data=df_time, x=group, y="sum_combined_per_sample_realtime", color = 'caller') +
                  labs(y = "time (min)") +
                  theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
                  facet.by = 'caller',
                  panel.labs.font = list(face = "bold"), nrow=1)
  
  line_cpuh <-  facet(ggboxplot(data=df_time, x=group, y="sum_combined_per_sample_cpuh", color = 'caller') +
                        labs(y = "CPUh", x = "#interval groups") +
                        theme(strip.text = element_blank()),
                        facet.by  = 'caller', nrow=1)
  
  bar <- facet(ggboxplot(data=df_storage, x=group, y="sum",color = 'caller') +
                 labs(y = "work dir (GB)") +
                 theme(legend.text = element_text(size = 8)) +
                 rremove("legend.title") + 
                 theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.text = element_blank()),
                facet.by  = 'caller', nrow=1)
  
  plot <- ggpubr::ggarrange(line, bar, line_cpuh, legend='none', ncol=1, align = "v")
  ann_plot <- annotate_figure(plot, top = text_grob(title, face = "bold", size = 14))
  ggsave(plot=ann_plot, filename = paste0(results_folder,outputname, ".png"), device="png", dpi=600 )
  ann_plot
}
#####################
######## AWS ########
#####################
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

#####################
###### Helper #######
#####################
bytesto <- function(bytes, to, bsize=1024){
  a <- data.frame ( unit = c('k', 'm', 'g', 't', 'p', 'e' ), factor = c(1,2,3,4,5,6), stringsAsFactors=FALSE)
  
  get_factor <- a$factor
  names(get_factor) <- a$unit
  
  return (bytes / (bsize ^ get_factor[to]))
}