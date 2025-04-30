import requests
from collections import defaultdict
import json
import pandas as pd
import os

# Set working directory to where KO parsing will occur
os.chdir("/home/doms/shotgun_mapping/KO_parsing")


def remove_redundant(x):
    """Remove redundant entries from a list while preserving order."""
    return list(dict.fromkeys(x))


def get_kegg_data(ko_ids):
    """
    Fetch KEGG annotation data for a list of KO IDs using the KEGG REST API.
    Returns a nested dict with keys: BRITE_KO, BRITE_Enzymes, Pathways, Reactions, GO Terms, Brite Hierarchy.
    """
    base_url = "http://rest.kegg.jp/get/"
    kegg_data = defaultdict(dict)

    for ko_id in ko_ids:
        # Build URL to request full KO entry details
        details_url = f"{base_url}{ko_id}"
        details_response = requests.get(details_url)

        if details_response.status_code == 200:
            # Split response by line for parsing
            details_data = details_response.text.split("\n")

            # Initialize placeholders for extracted fields
            extracted_data = {
                "symbol": "",
                "name": "",
                "pathways": [],
                "reaction": "",
                "brite": {
                    "KO": [],
                    "enzymes": [],
                    "hierarchy": defaultdict(list),  # Map parent Brite terms to sub-level terms
                },
                "cog": "",
                "go": "",
            }
            line_name = "NA"
            in_ko_section = True
            current_brite = None  # Track which BRITE term we're in

            # Parse each line for the fields of interest
            for line in details_data:
                if line.startswith("SYMBOL"):
                    extracted_data["symbol"] = line.split(maxsplit=1)[1].strip()
                    line_name = "SYMBOL"
                elif line.startswith("NAME"):
                    extracted_data["name"] = line.split(maxsplit=1)[1].strip()
                    line_name = "NAME"
                elif line.startswith("PATHWAY"):
                    line_name = "PATHWAY"
                    pathway = line.split(maxsplit=1)[1].strip()
                    extracted_data["pathways"].append(pathway)
                elif line.startswith("REACTION"):
                    line_name = "REACTION"
                    extracted_data["reaction"] = line.split(maxsplit=1)[1].strip()
                elif line.startswith("BRITE"):
                    # Enter BRITE section: reset tracking variables
                    line_name = "BRITE"
                    in_ko_section = True
                    current_brite = None
                    continue
                elif "COG" in line and ":" in line:
                    extracted_data["cog"] = line.split(":", maxsplit=1)[1].strip()
                elif "GO" in line and ":" in line:
                    extracted_data["go"] = line.split(":", maxsplit=1)[1].strip()
                elif line.startswith("GENES"):
                    # Stop before GENES section
                    break
                elif line.startswith("  "):
                    # Handle indented lines for PATHWAY continuation or BRITE hierarchy
                    if line_name == "BRITE" and in_ko_section:
                        # First level under BRITE: KO categories
                        if "Enzymes" in line:
                            # Switch to enzymes subsection
                            in_ko_section = False
                            continue
                        else:
                            current_brite = line.strip()
                            extracted_data["brite"]["KO"].append(current_brite)
                    elif line_name == "BRITE" and not in_ko_section:
                        # Under BRITE enzymes subsection
                        extracted_data["brite"]["enzymes"].append(line.strip())
                    elif line_name == "PATHWAY":
                        # Multi-line pathway entries
                        extracted_data["pathways"].append(line.strip())
                    else:
                        # Sub-level under a BRITE parent category
                        if current_brite:
                            extracted_data["brite"]["hierarchy"][current_brite].append(line.strip())

            # After parsing, store unique entries in the output dict
            kegg_data[ko_id]["BRITE_KO"] = remove_redundant(extracted_data["brite"]["KO"])
            kegg_data[ko_id]["BRITE_Enzymes"] = remove_redundant(extracted_data["brite"]["enzymes"])
            kegg_data[ko_id]["Pathways"] = remove_redundant(extracted_data["pathways"])
            kegg_data[ko_id]["Reactions"] = remove_redundant([extracted_data["reaction"]])
            kegg_data[ko_id]["GO Terms"] = remove_redundant([extracted_data["go"]])
            kegg_data[ko_id]["Brite Hierarchy"] = dict(extracted_data["brite"]["hierarchy"])

        else:
            # Log failures to fetch
            print(f"Failed to fetch data for {ko_id}: {details_response.status_code}")

    return kegg_data


def save_kegg_data(kegg_data, filename):
    """Write the KEGG data dict to a JSON file."""
    with open(filename, "w") as json_file:
        json.dump(kegg_data, json_file, indent=4)


def load_kegg_data(filename):
    """Load KEGG data dict from a JSON file."""
    with open(filename, "r") as json_file:
        return json.load(json_file)


def create_abundance_tables(kegg_results, abundance_df):
    """
    Using the KEGG annotation, aggregate KO abundances into BRITE, enzymes, pathways,
    reactions, and GO terms for each sample.
    """
    # Initialize nested dicts: {category: {sample: cumulative_abundance}}
    brite_KO_abundance = defaultdict(lambda: defaultdict(float))
    brite_Enzymes_abundance = defaultdict(lambda: defaultdict(float))
    pathway_abundance = defaultdict(lambda: defaultdict(float))
    reactions_abundance = defaultdict(lambda: defaultdict(float))
    go_abundance = defaultdict(lambda: defaultdict(float))

    # Iterate through each KO in the results
    for ko_id, data in kegg_results.items():
        if ko_id in abundance_df.index:
            # For each sample column, add the KO TPM to relevant categories
            for sample in abundance_df.columns:
                abundance_value = abundance_df.at[ko_id, sample]
                for brite_ko in data.get("BRITE_KO", []):
                    brite_KO_abundance[brite_ko][sample] += abundance_value
                for brite_enzymes in data.get("BRITE_Enzymes", []):
                    brite_Enzymes_abundance[brite_enzymes][sample] += abundance_value
                for pathway in data.get("Pathways", []):
                    pathway_abundance[pathway][sample] += abundance_value
                for reaction in data.get("Reactions", []):
                    reactions_abundance[reaction][sample] += abundance_value
                for go in data.get("GO Terms", []):
                    go_abundance[go][sample] += abundance_value

    return (
        brite_KO_abundance,
        brite_Enzymes_abundance,
        pathway_abundance,
        reactions_abundance,
        go_abundance,
    )


def save_abundance_to_csv(abundance_dict, filename):
    """Convert an abundance dict to a DataFrame and save it as CSV."""
    df = pd.DataFrame.from_dict(abundance_dict, orient="index").fillna(0)
    df.to_csv(filename)


# -------------------- Main workflow --------------------

# Load KEGG abundance data (KO TPMs per sample)
abundance_df = pd.read_csv("KO_tpm_F2.csv", index_col=0)

# Fetch KEGG annotation for all KO IDs in the abundance table
kegg_results = get_kegg_data(abundance_df.index)

# Save the raw KEGG annotation as JSON for reproducibility
save_kegg_data(kegg_results, "kegg_data.json")

# Build aggregated abundance tables by functional categories
(
    brite_abundance,
    brite_Enzymes_abundance,
    pathway_abundance,
    reactions_abundance,
    go_abundance,
) = create_abundance_tables(kegg_results, abundance_df)

# Write each aggregated table to CSV
save_abundance_to_csv(brite_abundance, "brite_abundance.csv")
save_abundance_to_csv(brite_Enzymes_abundance, "brite_enzymes_abundance.csv")
save_abundance_to_csv(pathway_abundance, "pathway_abundance.csv")
save_abundance_to_csv(reactions_abundance, "reactions_abundance.csv")
save_abundance_to_csv(go_abundance, "go_abundance.csv")
