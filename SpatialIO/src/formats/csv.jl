# SpatialIO/src/formats/csv.jl
# CSV helpers (tissue positions, CosMx flat files, etc.)

"""
    read_tissue_positions_csv(path::String) -> DataFrame

Parse a Visium `tissue_positions.csv` (or `.csv.gz`) into a DataFrame with
columns: barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres.
"""
function read_tissue_positions_csv(path::String)
    # TODO Phase 1: implement (handle both header and no-header variants from
    # different SpaceRanger versions)
    error("read_tissue_positions_csv: not yet implemented (Phase 1 deliverable)")
end
