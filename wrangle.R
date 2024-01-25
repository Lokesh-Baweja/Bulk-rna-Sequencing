# Step 1: Set the path to your directory containing the files
directory_path <- getwd()

# Step 2: Get the list of file names with a specific extension (e.g., .txt)
file_names <- list.files(path = directory_path, pattern = "\\.txt", full.names = TRUE)

# Step 3: Initialize an empty data frame to store the results
result_df <- data.frame(FirstColumn = character(), stringsAsFactors = FALSE)

# Step 4: Loop through each file
for (file in file_names) {
  # Read the data from the file
  current_data <- read.table(file, header = TRUE, stringsAsFactors = FALSE)

  # Extract the second column
  second_column <- current_data[, 2]

  # Extract the name after stripping the first string and the first underscore
  col_name <- sub("^[^_]*_", "", basename(file))
  col_name <- gsub("\\.txt$", "", col_name)  # Remove the ".txt" extension if present

  # Create a data frame with the first column and the extracted second column
  current_df <- data.frame(FirstColumn = current_data[, 1], 
                           SecondColumn = second_column,
                           stringsAsFactors = FALSE)

  # Set the column name
  names(current_df)[2] <- col_name

  # Merge the data frames based on the first column
  result_df <- merge(result_df, current_df, by = "FirstColumn", all = TRUE)
}

# Step 5: Print the resulting data frame
print(rownames(result_df))

print(colnames(result_df))

write.csv(result_df, file = "resulting_data.csv", row.names = FALSE)
