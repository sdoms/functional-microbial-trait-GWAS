import pandas as pd
import os

os.chdir("/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/")

# Load the data (assuming it's a CSV format; if not, modify accordingly)
file_path = "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/brite_abundance.csv"  # Replace with the correct file path
data = pd.read_csv(file_path, delimiter=",")


# Define function to extract rows based on prefix
def extract_rows_by_prefix(data, prefix):
    return data[data["Identifier"].str.startswith(prefix)]


# Extract different hierarchical levels
top_level = extract_rows_by_prefix(data, "091")
subcategories = extract_rows_by_prefix(data, "00")
# specific_genes = extract_rows_by_prefix(data, "K")

categories_all = extract_rows_by_prefix(data, "0")

# Save the extracted data to separate files
top_level.to_csv("top_level_091_BRITE.csv", index=False)
subcategories.to_csv("subcategories_00_BRITE.csv", index=False)
categories_all.to_csv("all_categories_BRITE.csv", index=False)


# Define function to create core subset
def create_core_subset(df, min_count=10, min_percentage=0.25):
    # Create a boolean DataFrame where values are True if greater than or equal to min_count
    presence_df = (
        df.iloc[:, 1:] >= min_count
    )  # Skip the 'Identifier' column for this calculation

    # Calculate the number of individuals where the function is present
    presence_count = presence_df.sum(axis=1)

    # Calculate the total number of individuals (columns)
    total_individuals = presence_df.shape[1]

    # Calculate the percentage of individuals where the function is present
    percentage_present = presence_count / total_individuals

    # Filter rows based on the percentage threshold
    core_subset = df[percentage_present >= min_percentage]

    return core_subset


top_level_core = create_core_subset(top_level)
subcategories_core = create_core_subset(subcategories)
categories_all_core = create_core_subset(categories_all)


# Define function to create core subset with new criteria
def create_core_subset_median(df, min_median=10, min_percentage=0.25):
    # Skip the 'Identifier' column for this calculation
    abundance_df = df.iloc[:, 1:]

    # Create a boolean DataFrame where values are True if greater than zero
    presence_df = abundance_df > 0

    # Calculate the number of individuals where the function is present
    presence_count = presence_df.sum(axis=1)

    # Calculate the total number of individuals (columns)
    total_individuals = presence_df.shape[1]

    # Calculate the percentage of individuals where the function is present
    percentage_present = presence_count / total_individuals

    # Calculate the median abundance of non-zero values for each row
    median_abundance = abundance_df.apply(lambda row: row[row > 0].median(), axis=1)

    # Filter rows based on the percentage threshold and median abundance
    core_subset = df[
        (percentage_present >= min_percentage) & (median_abundance >= min_median)
    ]

    return core_subset


top_level_core_median = create_core_subset_median(top_level)
subcategories_core_median = create_core_subset_median(subcategories)
categories_all_core_median = create_core_subset_median(categories_all)

top_level_core_median = create_core_subset_median(top_level, min_median=50)
subcategories_core_median = create_core_subset_median(subcategories, min_median=50)
categories_all_core_median = create_core_subset_median(categories_all, min_median=50)

# Save core subsets to CSV files (if needed)
top_level_core_median.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/top_level_core_median_50.csv",
    index=False,
)
subcategories_core_median.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/subcategories_core_median_50.csv",
    index=False,
)
categories_all_core_median.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/categories_all_core_median_50.csv",
    index=False,
)
