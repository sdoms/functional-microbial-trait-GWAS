import pandas as pd

# Load the data
file_path = "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/brite_enzymes_abundance.csv"  # Replace with the correct file path
df = pd.read_csv(file_path)

# Define regular expressions for hierarchical levels
patterns = {
    "level_1": r"^\d+\.\s",  # Matches 1. followed by space
    "level_2": r"^\d+\.\d+\s",  # Matches 1.1 followed by space
    "level_3": r"^\d+\.\d+\.\d+\s",  # Matches 1.1.1 followed by space
    "level_4": r"^\d+\.\d+\.\d+\.\d+\s",  # Matches 1.1.1.1 followed by space
}

# Filter rows based on patterns
df_level_1 = df[df["Identifier"].str.match(patterns["level_1"], na=False)]
df_level_2 = df[df["Identifier"].str.match(patterns["level_2"], na=False)]
df_level_3 = df[df["Identifier"].str.match(patterns["level_3"], na=False)]
df_level_4 = df[df["Identifier"].str.match(patterns["level_4"], na=False)]

# Save to CSV files (if needed)
df_level_1.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/level_1.csv",
    index=False,
)
df_level_2.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/level_2.csv",
    index=False,
)
df_level_3.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/level_3.csv",
    index=False,
)
df_level_4.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/level_4.csv",
    index=False,
)


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


# Apply the function to each DataFrame
df_level_1_core = create_core_subset(df_level_1)
df_level_2_core = create_core_subset(df_level_2)
df_level_3_core = create_core_subset(df_level_3)
df_level_4_core = create_core_subset(df_level_4)

# Save core subsets to CSV files (if needed)
df_level_1_core.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/level_1_core.csv",
    index=False,
)
df_level_2_core.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/level_2_core.csv",
    index=False,
)
df_level_3_core.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/level_3_core.csv",
    index=False,
)
df_level_4_core.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/level_4_core.csv",
    index=False,
)


# Define function to create core subset with new criteria
def create_core_subset_median(df, min_median=50, min_percentage=0.25):
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


# Apply the function to each DataFrame
df_level_1_core = create_core_subset_median(df_level_1)
df_level_2_core = create_core_subset_median(df_level_2)
df_level_3_core = create_core_subset_median(df_level_3)
df_level_4_core = create_core_subset_median(df_level_4)

# Save core subsets to CSV files (if needed)
df_level_1_core.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/EC_level_1_core.csv",
    index=False,
)
df_level_2_core.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/EC_level_2_core.csv",
    index=False,
)
df_level_3_core.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/EC_level_3_core.csv",
    index=False,
)
df_level_4_core.to_csv(
    "/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_mapping/EC_level_4_core.csv",
    index=False,
)
