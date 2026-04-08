import importlib.resources as pkg_resources
import os
from datetime import datetime

import pandas as pd
from fastmcp import FastMCP
from fastmcp.exceptions import ToolError

from celltypeschema import CellMarkerDB, SingleCellBaseParam

db_mcp = FastMCP("BioMarkerMCP-DB-Server")


@db_mcp.tool(tags={"celltype"})
def query_cancer_celltype_marker(request: CellMarkerDB):
    """query the celltype marker in cancer from cellmarker database"""
    try:
            db_file = "C:/Users/tangbx/AppData/Local/saya/dingent/plugins/cellmarker/Cell_marker_All.csv"
       
            if not os.path.exists(db_file):
                raise FileNotFoundError(f"Database file not found at {db_file}")
            db_df = pd.read_csv(db_file)
    except Exception as e:
        raise ToolError(f"Failed to load database: {str(e)}")

    # Apply filters based on provided fields
    if request.species:
        db_df = db_df.loc[db_df["species"].isin(request.species)]

    if request.tissue_class:
        db_df = db_df.loc[db_df["tissue_class"].isin(request.tissue_class)]

    if request.tissue_type:
        db_df = db_df.loc[db_df["tissue_type"].isin(request.tissue_type)]

    if request.cancer_type:
        cancer_types = db_df["cancer_type"].unique()
        db_df = db_df.loc[db_df["cancer_type"].isin(request.cancer_type)]
        if len(db_df) == 0:
            return f"404 NOT FOUND ERROR: No records found for cancer_type: {request.cancer_type}, available cancer_types: {cancer_types}"

    if request.cell_type:
        db_df = db_df.loc[db_df["cell_type"].isin(request.cell_type)]

    if request.cell_name:
        cellnames = db_df.loc[:, "cell_name"].unique()
        db_df = db_df.loc[db_df["cell_name"].isin(request.cell_name)]
        if len(db_df) == 0:
            return f"404 NOT FOUND ERROR: No records found for cell_name: {request.cell_name}, available cell_names: {cellnames}"

    if request.Symbol:
        db_df = db_df.loc[db_df["Symbol"].isin(request.Symbol)]

    if request.Genetype:
        db_df = db_df.loc[db_df["Genetype"].isin(request.Genetype)]

    if request.GeneID:
        db_df = db_df.loc[db_df["GeneID"].isin(request.GeneID)]

    # Get filtered results
    result = db_df.loc[:, request.show_columns].head(request.show_num)

    # Create output directory if it doesn't exist
    output_dir = "query_results"
    os.makedirs(output_dir, exist_ok=True)

    # Generate filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"cellmarker_query_{timestamp}.csv"
    output_path = os.path.join(output_dir, filename)

    # Write results to CSV file
    db_df.loc[:, request.show_columns].to_csv(output_path, index=False)

    return {
        "query_parameters": request.model_dump(),
        "query_result": {
            f"head_{request.show_num}_records": result.to_dict(),
            "full_records_output_file": output_path,
        },
    }


@db_mcp.tool(tags={"celltype"})
def query_normal_celltype_marker(request: SingleCellBaseParam):
    """query the singlecellbase from singlecellbase database"""
    try:
            db_file = "C:/Users/tangbx/AppData/Local/saya/dingent/plugins/cellmarker/singleCellBase_20230904_ALL.txt"
       
            if not os.path.exists(db_file):
                raise FileNotFoundError(f"Database file not found at {db_file}")
            db_df = pd.read_csv(db_file, sep="\t")
    except Exception as e:
        raise ToolError(f"Failed to load database: {str(e)}")

    if request.kingdom:
        if request.kingdom not in db_df["kingdom"].unique():
            return f"404 NOT FOUND ERROR: No records found for kingdom: {request.kingdom}, available kingdoms: {db_df['kingdom'].unique()}"
        db_df = db_df.loc[db_df["kingdom"] == request.kingdom]

    if request.phylum:
        if request.phylum not in db_df["phylum"].unique():
            return f"404 NOT FOUND ERROR: No records found for phylum: {request.phylum}, available phylums: {db_df['phylum'].unique()}"
        db_df = db_df.loc[db_df["phylum"] == request.phylum]

    if request.classes:
        if request.classes not in db_df["classes"].unique():
            return f"404 NOT FOUND ERROR: No records found for classes: {request.classes}, available classes: {db_df['classes'].unique()}"
        db_df = db_df.loc[db_df["classes"] == request.classes]

    if request.species:
        if request.species not in db_df["species"].unique():
            return f"404 NOT FOUND ERROR: No records found for species: {request.species}, available species: {db_df['species'].unique()}"
        db_df = db_df.loc[db_df["species"] == request.species]

    if request.sample_type:
        if request.sample_type not in db_df["sample_type"].unique():
            return f"404 NOT FOUND ERROR: No records found for sample_type: {request.sample_type}, available sample_types: {db_df['sample_type'].unique()}"
        db_df = db_df.loc[db_df["sample_type"] == request.sample_type]

    if request.tissue_type:
        if request.tissue_type not in db_df["tissue_type"].unique():
            return f"404 NOT FOUND ERROR: No records found for tissue_type: {request.tissue_type}, available tissue_types: {db_df['tissue_type'].unique()}"
        db_df = db_df.loc[db_df["tissue_type"] == request.tissue_type]

    if request.cell_type:
        if request.cell_type not in db_df["cell_type"].unique():
            return f"404 NOT FOUND ERROR: No records found for cell_type: {request.cell_type}, available cell_types: {db_df['cell_type'].unique()}"
        db_df = db_df.loc[db_df["cell_type"] == request.cell_type]

    # Get filtered results
    result = db_df.loc[:, request.show_columns].head(request.show_num)

    # Create output directory if it doesn't exist
    output_dir = "query_results"
    os.makedirs(output_dir, exist_ok=True)

    # Generate filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"singlecellbase_query_{timestamp}.csv"
    output_path = os.path.join(output_dir, filename)

    # Write results to CSV file
    db_df.loc[:, request.show_columns].to_csv(output_path, index=False)

    return {
        "query_parameters": request.model_dump(),
        "query_result": {
            f"head_{request.show_num}_records": result.to_dict(),
            "full_records_output_file": output_path,
        },
    }

if __name__ == "__main__":
    db_mcp.run()