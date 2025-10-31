# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 16:15:33 2025

@author: L-F-S
"""

from conf import CHEMBL_INPUT_DATA_DIR, cell_lines
import requests
import pandas as pd
from tqdm import tqdm
from collections import defaultdict

activity_type = "IC50"

# API endpoints
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
CELL_URL = f"{BASE_URL}/cell_line.json"
ACTIVITY_URL = f"{BASE_URL}/activity.json"
ASSAY_URL = f"{BASE_URL}/assay.json"

#%% Functions

def get_total(json_page):
    return json_page.get("page_meta", {}).get("total_count", 0)

def download_cell_line_map(limit=1000):
    cell_line_name_2_chembl = {  }
    cell_line_chembl_2_name = {  }
    offset = 0
    # First request just to get total count
    initial_params = {"limit": limit, "offset": offset}
    initial_response = requests.get(CELL_URL, params=initial_params)
    initial_response.raise_for_status()
    initial_page = initial_response.json()
    initial_cells = initial_page['cell_lines']
    for cell in initial_cells:
        cell_line_name_2_chembl[cell["cell_name"]] = cell["cell_chembl_id"]
        cell_line_chembl_2_name[cell["cell_chembl_id"]] = cell["cell_name"]
    
    total = get_total(initial_page)
    print(f"Total cells available: {total}")
    offset += limit

    with tqdm(total=total, desc="Fetching cell lines") as pbar:
        pbar.update(len(initial_page["cell_lines"]))
        while len(cell_line_name_2_chembl) < total:
            params = {"limit": limit, "offset": offset}
            response = requests.get(CELL_URL, params=params)
            page = response.json()
            cells = page["cell_lines"]
            
            for cell in cells:
                cell_line_name_2_chembl[cell["cell_name"]] = cell["cell_chembl_id"]
                cell_line_chembl_2_name[cell["cell_chembl_id"]] = cell["cell_name"]
            if len(cells) < limit:
                print(len(cells))
                break
            pbar.update(len(cells))
            offset += limit
        
    return cell_line_name_2_chembl, cell_line_chembl_2_name


def fetch_assays_for_cell(cell_chembl_id, limit=1000):
    """Fetch assay_chembl_ids associated with a specific cell_chembl_id"""
    assays = []
    offset = 0

    # First request to get total count
    initial_params = {
        "cell_chembl_id": cell_chembl_id,
        "limit": limit,
        "offset": offset
    }
    initial_response = requests.get(ASSAY_URL, params=initial_params)
    initial_response.raise_for_status()
    initial_page = initial_response.json()
    total = get_total(initial_page)
    print(f"total assays available: {total}")
    assays.extend(initial_page.get("assays", []))
    offset += limit

    with tqdm(total=total, desc=f"Fetching assays for {cell_chembl_id}") as pbar:
        pbar.update(len(assays))
        while len(assays) < total:
            params = {
                "cell_chembl_id": cell_chembl_id,
                "limit": limit,
                "offset": offset
            }
            response = requests.get(ASSAY_URL, params=params)
            response.raise_for_status()
            page = response.json()
            data = page.get("assays", [])
            assays.extend(data)
            pbar.update(len(data))
            if len(data) < limit:
                break
            offset += limit

    # Extract only the assay_chembl_ids
    assay_ids = [a["assay_chembl_id"] for a in assays]
    return assay_ids


def fetch_ic50_data_for_assay(assay_id, limit=1000):
    """Fetch IC50 activity data for a given cell line"""
    activities = []
    offset = 0
    
    # First request to get total count
    initial_params = {
        "assay_chembl_id": assay_id,
        "standard_type": "IC50",
        "limit": limit,
        "offset": offset
    }
    initial_response = requests.get(ACTIVITY_URL, params=initial_params)
    initial_response.raise_for_status()
    initial_page = initial_response.json()
    total = get_total(initial_page)
    # print(f"total assays available: {total}")
    
    activities.extend(initial_page.get("activities", []))

    total = initial_page.get("page_meta", {}).get("total_count", 0)
    if total == 0:
        return None
    # with tqdm(desc=f"Fetching IC50 for {assay_id}") as pbar:
        # pbar.total = total
        # tab code below for pbar
    while True:
        params = {
            "assay_chembl_id": assay_id,
            "standard_type": "IC50",
            "limit": limit,
            "offset": offset
        }
        response = requests.get(ACTIVITY_URL, params=params)
        # print("Querying:", response.url)
        response.raise_for_status()
        page = response.json()
        data = page.get("activities", [])
        
        activities.extend(data)
        
        # pbar.update(len(data))
        if len(data) < limit:
            break
        offset += limit

    return activities

#%%

# initialize dictionary of cell line - cell chembl id
cell_line_name_2_chembl, cell_line_chembl_2_name = download_cell_line_map()

#%% Find any cell with MCF7 in the label
cell_lines = ["MCF7","HepG2","HT-29"] # correct nomenclature for chembl
for cell_name, cell_id in list(cell_line_name_2_chembl.items()):
    if "MCF7" in cell_name:
        print(cell_name, "->", cell_id)
for cell_name, cell_id in list(cell_line_name_2_chembl.items()):
    if "HepG2" in cell_name:
        print(cell_name, "->", cell_id)
for cell_name, cell_id in list(cell_line_name_2_chembl.items()):
    if "HT-29" in cell_name:
        print(cell_name, "->", cell_id)

#%% fetch assays for cell lines
assays_of_cell = {}
for cell_name in cell_lines:
    cell_chembl_id = cell_line_name_2_chembl[cell_name]   
    assays_of_cell[cell_name] = fetch_assays_for_cell(cell_chembl_id, limit=1000)
#%%
cell_name = "HT-29"
cell_chembl_id = cell_line_name_2_chembl[cell_name]   
assays_of_cell[cell_name] = fetch_assays_for_cell(cell_chembl_id, limit=1000)
#%% save assays
for cell_line, assays in assays_of_cell.items():
    with open(CHEMBL_INPUT_DATA_DIR+cell_line+'_assays.txt', 'w') as f:
        string = '\n'.join(assays)
        f.write(string)
#%% load assays
# TODO implement

#%%OPTIONAL print total number of activities for assays:
# for cell_name in cell_lines[:1]:

#     # discard duplicate activity ids
#     activity_ids = []
#     all_acts = 0
#     print(f"\nProcessing cell line: {cell_name}, id {cell_chembl_id}")
#     for assay_id in tqdm(assays_of_cell[cell_name], desc="Fetching IC50 activities"):
#         params = {"assay_chembl_id": assay_id, "standard_type": "IC50", "limit": 1000}
#         response = requests.get("https://www.ebi.ac.uk/chembl/api/data/activity.json", params=params)
#         page = response.json()
#         total = get_total(page)
#         print(assay_id, total)
#         all_acts+=total
#     print(all_acts, cell_name)
#%%
#%% Get IC50 activity data for all cell lines
ic50_of_cell=defaultdict(list)
for cell_name in ["HT-29"]:#cell_lines[1:2]:

    # discard duplicate activity ids
    activity_ids = set()
    document_chembl_ids = set()
    
    print(f"\nProcessing cell line: {cell_name}, id {cell_chembl_id}")
    for n, assay_id in enumerate(tqdm(assays_of_cell[cell_name], desc="Fetching IC50 activities")):
        
        # get all chembl activity data for assay_id
        ic50_data_list = fetch_ic50_data_for_assay(assay_id)
        # print(f"Retrieved {len(ic50_data)} IC50 records")
        
        if ic50_data_list: # is None if no activity is recorded for given assay_id
            for ic50_data in ic50_data_list:
                activity_id=ic50_data['activity_id']
                document_chembl_id=ic50_data['document_chembl_id']
                if (not activity_id in activity_ids) and (not document_chembl_id in document_chembl_ids):
                    activity_ids.add(activity_id)
                    ic50_of_cell[cell_name].append(ic50_data)
            
#%% Convert to DataFrame and save
df_of={}
for cell_line in cell_lines[1:]:
    print(cell_line)
    df_of[cell_line] = pd.DataFrame(ic50_of_cell[cell_line])
    df_of[cell_line].to_csv(CHEMBL_INPUT_DATA_DIR+cell_line+"_activity_second_part.tsv", sep='\t', index=False)
    print(cell_line, f"Data saved to {CHEMBL_INPUT_DATA_DIR}{cell_line}_activity.tsv")


