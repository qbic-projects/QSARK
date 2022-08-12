format_input <- function(cost_table) {
  
  # transpose data frame & use first column from original dataframe as column headers in transposed dataframe
  cost_table_formatted <- setNames(data.frame(t(cost_table[,-1])), cost_table[,1])

  # create a new column named 'type' that has the current row names (needed for ggplot2)
  cost_table_formatted <- tibble::rownames_to_column(cost_table_formatted, "type")

  # Round too last 3 digits & remove costs below 0.00$
  cost_table_formatted <- cost_table_formatted %>% mutate_if(is.numeric, round, digits=3)
  cost_table_formatted <- filter(cost_table_formatted, `2022-08-07` > 0.00 , .preserve = FALSE)
  return (cost_table_formatted)
}