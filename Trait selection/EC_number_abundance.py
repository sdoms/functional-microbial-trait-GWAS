#!/usr/bin/env python3

import requests
import pandas as pd
import re


# Load KO abundance table
abundance_df = pd.read_csv("/Users/doms/Documents/Research/Experiments/shotgun_hybrids/KO_tpm_F2.csv", index_col=0)



def get_ec_numbers(ko_id):
    ec_numbers = []
    try:
        # Getting the full information for the KO
        ec_response = requests.get(f"http://rest.kegg.jp/get/{ko_id}")
        if ec_response.status_code == 200:
            # Use a regular expression to find EC numbers with 1 to 4 parts
            ec_numbers = re.findall(r"EC\s*:\s*(\d+\.\d+\.\d*\.?\d*)", ec_response.text)
    except Exception as e:
        print(f"Error processing KO {ko_id}: {e}")
    
    return ec_numbers

# Initialize a dictionary to hold the aggregated EC abundance data
ec_abundance = {}

# Iterate through each KO in the abundance table
for ko_id in abundance_df.index:
    ec_numbers = get_ec_numbers(ko_id)
    
    # Aggregate abundance by EC numbers
    for ec in ec_numbers:
        if ec not in ec_abundance:
            ec_abundance[ec] = abundance_df.loc[ko_id]
        else:
            ec_abundance[ec] += abundance_df.loc[ko_id]

# Convert the EC abundance dictionary to a DataFrame
ec_df = pd.DataFrame(ec_abundance).T

# Save the EC abundance table
ec_df.to_csv("/Users/doms/Documents/Research/Experiments/shotgun_hybrids/ec_abundance.csv")