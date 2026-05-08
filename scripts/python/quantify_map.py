# Created By: Faisal Isam
# Last Updated: May 4, 2026
# Input: taxa directory, mine directory, ecoreg file
# Required File: "ca.json" 

import os
import glob
import sys
import itertools                        

from math import pi                     # For creating pi charts
import pandas as pd                     # For working with data
import geopandas as gpd                 # For geospatial data manipulation
from shapely.geometry import Point      # To create geometry objects from lat/lon

# Tools for creating the dashboard
import panel as pn
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from bokeh.palettes import Category20
from bokeh.transform import cumsum
from bokeh.models import Legend, LegendItem

# Function to load and combine csv files
def load_csvs(directory, label):

    # Checks if the directory exists 
    if not os.path.isdir(directory):
        print(f"Error: Directory {directory} does not exist.")
        sys.exit(1)

    # Retrieves all files in the given directory
    files = glob.glob(os.path.join(directory, '*.csv'))
    if not files:
        print(f"No CSV files found in {directory}")
        sys.exit(1)

    dfs = []    # Temporary list to store each CSV as a DataFrame

    # Loop to read each file in the directory
    for i, f in enumerate(files, 1):
        print(f"[{i}/{len(files)}] Reading {label} file: {os.path.basename(f)}")
        try:
            df = pd.read_csv(f, encoding='latin1', low_memory=False)    # Reads each file
            df.columns = df.columns.str.strip().str.lower()             # Converts everything to lowercase and removes any leading or trailing characters
            dfs.append(df)                                              # Adds the processed file into the dataframe
        except Exception as e:
            print(f"Warning: Could not read {f}, skipping. Error: {e}")

    # Combines all the files into a single data frame
    combined = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(combined)} rows of {label} data from {len(dfs)} files")

    return combined

# Function to process the taxa data frame
def preprocess_taxa(df):

    # Only run cleaning if required columns exist
    if {"lat", "lon", "family_name"}.issubset(df.columns):
        # Remove rows with missing or invalid coordinates or family names
        df = df.dropna(subset=["lat", "lon", "family_name"])
        df = df[df["lat"] != "N/A"]
        df = df[df["lon"] != "N/A"]

        # Converts coordinates to numeric values, coercing errors to NaN
        df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
        df["lon"] = pd.to_numeric(df["lon"], errors="coerce")

        # Drop rows that still have missing coordinate values
        df = df.dropna(subset=["lat", "lon"])

        # Remove duplicate sampling points (same family and location)
        df = df.drop_duplicates(subset=["family_name", "lat", "lon"])

    return df

# Function to process the mine data frame
def preprocess_mine(df):

    # Initialize coordinate storage and counters
    mine_coords, pper_coords, cs_coords = [], [], []
    total_cs_sites, total_mine_sites, total_pper_sites = 0, 0, 0

    # EcoReg DNA Research Sites 
    if {"station_latitude", "station_longitude", "metal", "mine"}.issubset(df.columns):

        # Remove duplicates and invalid coordinates
        cs_df = df.dropna(subset=["station_latitude", "station_longitude", "metal", "mine"])
        cs_df = cs_df[cs_df["station_latitude"] != "N/A"]
        cs_df = cs_df[cs_df["station_longitude"] != "N/A"]
        cs_df = cs_df.drop_duplicates(subset=["mine", "station_latitude", "station_longitude", "metal"])

        # Iterate through each valid record
        for _, row in cs_df.iterrows():
            lat, lon = row.get("station_latitude"), row.get("station_longitude")
            if pd.notna(lat) and pd.notna(lon):
                metal = row.get("metal")
                site_name = row.get("mine")
                if pd.notna(metal):
                    try:
                        # EcoReg site count is fixed by project design
                        total_cs_sites = 6
                        cs_coords.append({"name": str(site_name), "lat": float(lat), "lon": float(lon), "metal": str(metal)})
                    except:
                        pass

    # MDMER Mines
    if {"latitude", "longitude", "metal", "facility name/nom de l'installation"}.issubset(df.columns):

        # Remove duplicates and invalid coordinates
        mine_df = df.dropna(subset=["latitude", "longitude", "metal", "facility name/nom de l'installation"])
        mine_df = mine_df[mine_df["latitude"] != "N/A"]
        mine_df = mine_df[mine_df["longitude"] != "N/A"]
        mine_df = mine_df.drop_duplicates(subset=["facility name/nom de l'installation", "latitude", "longitude", "metal"])

        # Iterate through each valid record
        for _, row in mine_df.iterrows():
            lat, lon = row.get("latitude"), row.get("longitude")
            if pd.notna(lat) and pd.notna(lon):
                metal = row.get("metal")
                site_name = row.get("facility name/nom de l'installation")
                if pd.notna(metal):
                    try:
                        total_mine_sites += 1
                        mine_coords.append({"name": str(site_name), "lat": float(lat), "lon": float(lon), "metal": str(metal)})
                    except:
                        pass

    # PPER Mills (coordinates stored in one text field)
    if {"coordinates", "product", "company"}.issubset(df.columns):

        # Remove duplicates and invalid coordinates
        pper_df = df.dropna(subset=["coordinates", "product", "company"])
        pper_df = pper_df[pper_df["coordinates"] != "N/A"]
        pper_df = pper_df.drop_duplicates(subset=["company", "coordinates", "product"])

        # Iterate through each valid record
        for _, row in pper_df.iterrows():
            coords = row.get("coordinates")
            product = row.get("product")
            site_name = row.get("company")
            if pd.notna(coords) and pd.notna(product):
                try:
                    parts = coords.split(",")
                    if len(parts) == 2:
                        lat_f, lon_f = float(parts[0].strip()), float(parts[1].strip())
                        total_pper_sites += 1
                        pper_coords.append({"name": str(site_name), "lat": lat_f, "lon": lon_f, "product": str(product)})
                except:
                    pass

    # Converts the lists into DataFrames
    return (
        pd.DataFrame(mine_coords),
        pd.DataFrame(pper_coords),
        pd.DataFrame(cs_coords),
        total_cs_sites,
        total_mine_sites,
        total_pper_sites,
    )

# Function to measure nearby taxa density around mines
def quantify_points(taxa_df, mine_df, site_type):
    # Convert taxa and mine dataframes to GeoDataFrames with point geometry
    # Assumes lon/lat are already cleaned + numeric upstream
    taxa_gdf = gpd.GeoDataFrame(
        taxa_df,
        geometry=[Point(xy) for xy in zip(taxa_df['lon'], taxa_df['lat'])],
        crs="EPSG:4326"     # WGS84 (lat/lon)
    )

    # Convert mine DataFrame into a GeoDataFrame with point geometry
    mine_gdf = gpd.GeoDataFrame(
        mine_df,
        geometry=[Point(xy) for xy in zip(mine_df['lon'], mine_df['lat'])],
        crs="EPSG:4326"
    )

    # Projects both datasets into a CRS that uses meters
    # Required so the 5000 buffer distance is interpreted as meters (not degrees)
    taxa_gdf = taxa_gdf.to_crs(epsg=3347)
    mine_gdf = mine_gdf.to_crs(epsg=3347)

    results = []    # List to store output for each mine

    # Loop through each mine and count nearby taxa within a 5 km radius
    for _, mine_row in mine_gdf.iterrows():
        # Create a circular buffer of 5 km around the mine/mill point
        mine_point = mine_row.geometry
        buffer = mine_point.buffer(5000)  # 5 km radius

        # Filter taxa points that fall within this buffer
        # `.within()` excludes boundary points by design (acceptable here)
        nearby_taxa = taxa_gdf[taxa_gdf.geometry.within(buffer)]
        
        # Total counts (including duplicates)
        taxa_count_total = len(nearby_taxa)
        species_total = len(nearby_taxa["species_name"].dropna())
        bins_total = len(nearby_taxa["bin_uri"].dropna())

        # Unique counts (removing duplicates)
        taxa_list_unique = list(nearby_taxa['family_name'].unique())
        taxa_count_unique = len(taxa_list_unique)

        # Remove generic placeholder species labels (e.g. "sp.")
        # Done AFTER totals so raw counts are preserved
        nearby_taxa = nearby_taxa[~nearby_taxa["species_name"].str.contains("sp", na=False)]

        # Unique species and BINs after cleaning
        species_unique = nearby_taxa["species_name"].dropna().unique().tolist()
        num_species_unique = len(species_unique)

        bins_unique = nearby_taxa["bin_uri"].dropna().unique().tolist()
        num_bins_unique = len(bins_unique)

        # Append the results for this mine into the lists
        results.append({
            'site_name': mine_row['name'],                              # Name of the site 
            'mine_lat': mine_row['lat'],                                # Latitude of the mine
            'mine_lon': mine_row['lon'],                                # Longitude of the mine
            'num_taxa_total': taxa_count_total,                         # Total taxa count (including duplicates)
            'taxa_list_total': list(nearby_taxa['family_name']),        # List of all taxa (including duplicates)
            'num_taxa_unique': taxa_count_unique,                       # Unique taxa count
            'taxa_list_unique': taxa_list_unique,                       # List of unique taxa families
            'num_records_with_species_name': species_total,             # Total species count (including duplicates)
            'species_total': list(nearby_taxa["species_name"].dropna()),# List of all species (including duplicates)
            'num_species_unique': num_species_unique,                   # Unique species count
            'species_unique': species_unique,                           # List of unique species
            'num_bins_total': bins_total,                               # Total BINs count (including duplicates)
            'bin_ids_total': list(nearby_taxa["bin_uri"].dropna()),     # List of all BINs (including duplicates)
            'num_bins_unique': num_bins_unique,                         # Unique BINs count
            'bin_ids_unique': bins_unique,                              # List of unique BINs
            'metal/product': mine_row.get('metal') if pd.notna(mine_row.get('metal')) else mine_row.get('product'),
            'site_type': site_type
        })

    # Convert results into a DataFrame for visualization
    return pd.DataFrame(results)


def summarize_ecoreg(df):

    df["family_name"] = None

    if "family_name" in df.columns:
        df["family_name"] = df["family_name"]

    if "family" in df.columns:
        df["family_name"] = df["family_name"].fillna(df["family"])


    if "species_name" not in df.columns:
        if "species" in df.columns:
            df["species_name"] = df["species"]
        else:
            df["species_name"] = None

    if "bin_uri" not in df.columns:
        if "bin" in df.columns:
            df["bin_uri"] = df["bin"]
        else:
            df["bin_uri"] = None

    df["latitude"] = pd.to_numeric(df["latitude"], errors="coerce")
    df["longitude"] = pd.to_numeric(df["longitude"], errors="coerce")

    # Normalize site names for grouping
    if "site_name" not in df.columns:
        print("Warning: site_name missing, assigning 'Unknown'")
        df["site_name"] = "Unknown"

    df["site_name"] = df["site_name"].astype(str).str.strip()
    df = df[~df["site_name"].isin(["nan", "None", ""])]
    df = df.dropna(subset=["site_name"])

    def ensure_list(x):
        if isinstance(x, list):
            return x
        if pd.isna(x) or str(x).strip() == "":
            return []
        return [str(x).strip()]

    # Wrap scalar values into lists so aggregation via sum() works cleanly
    df["family_name"] = df["family_name"].apply(ensure_list)
    df["species_name"] = df["species_name"].apply(ensure_list)
    df["bin_uri"] = df["bin_uri"].apply(ensure_list)

    def flatten_lists(series):
        result = []
        for item in series:
            if isinstance(item, list):
                result.extend(item)
            elif pd.notna(item) and str(item).strip() != "":
                result.append(str(item).strip())
        return result

    grouped = df.groupby("site_name", as_index=False).agg({
        "family_name": flatten_lists,
        "species_name": flatten_lists,
        "bin_uri": flatten_lists,
        "latitude": "first",
        "longitude": "first",
        "metal": "first" if "metal" in df.columns else "first"
    })

    results = []
    for _, row in grouped.iterrows():
        families = row["family_name"]
        species = row["species_name"]
        bins = row["bin_uri"]

        # Total counts (including duplicates)
        taxa_count_total = len(families)
        species_total = len(species)
        bins_total = len(bins)

        # Unique counts (removing duplicates)
        taxa_list_unique = list(set(families))
        taxa_count_unique = len(taxa_list_unique)

        species_unique = list(set(species))
        num_species_unique = len(species_unique)

        bins_unique = list(set(bins))
        num_bins_unique = len(bins_unique)

        results.append({
            "site_name": row["site_name"],                  # Name of the site 
            "mine_lat": row["latitude"],                    # Latitude of the mine
            "mine_lon": row["longitude"],                   # Longitude of the mine
            "num_taxa_total": taxa_count_total,             # Total taxa count
            "taxa_list_total": families,                    # List of all taxa (including duplicates)
            "num_taxa_unique": taxa_count_unique,           # Unique taxa count
            "taxa_list_unique": taxa_list_unique,           # List of unique taxa
            "num_records_with_species_name": species_total, # Total species count
            "species_total": species,                       # List of all species (including duplicates)
            "num_species_unique": num_species_unique,       # Unique species count
            "species_unique": species_unique,               # List of unique species
            "num_bins_total": bins_total,                   # Total BINs count
            "bin_ids_total": bins,                          # List of all BINs (including duplicates)
            "num_bins_unique": num_bins_unique,             # Unique BINs count
            "bin_ids_unique": bins_unique,                  # List of unique BINs
            "metal/product": row["metal"],
            "site_type": "Novel Data"
        })

    return pd.DataFrame(results)


# Function to generate a consistent colour map
def generate_colour_map(all_names):
    n = len(all_names)                                                                  # Determine the number of unique category names
    palette = Category20[20] if n <= 20 else (Category20[20] * ((n // 20) + 1))         # Use Bokeh's Category20 palette, repeat it if more than 20 items are needed

    return {name: palette[i] for i, name in enumerate(sorted(all_names))}               # Return a dictionary mapping each sorted name to a colour

# Function to create a pie chart
def make_pie_chart(data, column, title, colour_map=None):
    # Flatten nested lists into a Series
    counts = (
        pd.Series([item for sublist in data[column] for item in sublist if isinstance(sublist, list)])
        .value_counts()
        .reset_index()
    )
    counts.columns = ["name", "count"]

    # Convert counts to angles (radians)
    counts["angle"] = counts["count"] / counts["count"].sum() * 2 * pi

    # Colors (use provided map, fallback to Category20)
    if colour_map:
        counts["colour"] = counts["name"].map(colour_map).fillna("grey")
    else:
        palette = Category20[20] if len(counts) <= 20 else (Category20[20] * ((len(counts) // 20) + 1))
        counts["colour"] = palette[:len(counts)]

    n = len(counts)

    # Dynamically adjust legend size/layout based on number of categories
    if n <= 30:
        legend_cols = 1
        fig_width = 550
    elif n <= 75:
        legend_cols = 2
        fig_width = 750
    elif n <= 130:
        legend_cols = 3
        fig_width = 950
    else:
        legend_cols = 8
        fig_width = 1450

    # Initialize Bokeh figure
    p = figure(
        width=fig_width,
        height=500,
        title=title,
        toolbar_location=None,
        x_range=(-0.6, 0.6),        # Fixed range for centered pie
    )

    # Draw pie wedges
    r = p.wedge(
        x=0, y=1, radius=0.4,
        start_angle=cumsum("angle", include_zero=True),
        end_angle=cumsum("angle"),
        line_color="white",
        fill_color="colour",
        source=counts
    )

    # Tooltip (italicize species only)
    hover = HoverTool()
    if column == "species_total":
        hover.tooltips = [("Name", "<i>@name</i>"), ("Count", "@count")]
    else:
        hover.tooltips = [("Name", "@name"), ("Count", "@count")]
    p.add_tools(hover)

    # Remove axes/grid for clean pie appearance
    p.axis.visible = False
    p.grid.visible = False

    # Build explicit legend (avoids warnings)
    legend_items = [
        LegendItem(label=row["name"], renderers=[r], index=i)
        for i, row in counts.iterrows()
    ]

    legend = Legend(
        items=legend_items,
        ncols=legend_cols,
        orientation="vertical",
        glyph_width=12,
        glyph_height=12,
        spacing=2,
        label_standoff=6,

        border_line_color="grey",
        border_line_width=1,
        border_line_alpha=0.4,
        padding=10,

        label_text_font_style="italic" if column == "species_total" else "normal",
    )

    p.add_layout(legend, "right")
    legend.location = "top"

    return p

# Function to generate an HTML summary card showing key metrics
def make_summary(results_df, total_sites):

    # Flatten nested lists of unique taxa across all sites
    flattened_taxa_list = list(itertools.chain(*results_df["taxa_list_unique"]))
    
    # Convert to pandas Series before using pd.unique()
    unique_taxa = pd.Series(flattened_taxa_list).nunique()
    
    total_sites_with_nearby_taxa = len(results_df)
    total_taxa = results_df["num_taxa_total"].sum()
    avg_taxa = results_df["num_taxa_total"].mean()
    max_taxa = results_df["num_taxa_total"].max()
    
    total_bin = results_df["num_bins_total"].sum()
    avg_bin = results_df["num_bins_total"].mean()
    max_bin = results_df["num_bins_total"].max()
    
    # Flatten the list of lists in the 'bin_ids_unique' column
    flattened_bin_list = list(itertools.chain(*results_df["bin_ids_unique"]))
    
    # Compute unique BINs after flattening the list
    unique_bin = pd.Series(flattened_bin_list).nunique()

    # Build an HTML panel using inline CSS for styling
    html = f"""
    <div style="background-color:#f8f9fa;
                border-radius:20px;
                padding:20px;
                margin-top: 15px;
                box-shadow:0 2px 6px rgba(0,0,0,0.1);
                width:300px;">
        <h3 style="text-align:center; color:#2c3e50;">Summary</h3>
        <hr style="border:1px solid #e0e0e0;">
            <p style="font-size:16px;"><b>Total Sites:</b>{total_sites}</p>
        <p style="font-size:16px;"><b>Total Sites With Nearby Taxa:</b> {total_sites_with_nearby_taxa}</p>
        <p style="font-size:16px;"><b>Total Nearby Taxa:</b> {total_taxa:,}</p>
        <p style="font-size:16px;"><b>Avg Taxa per Site:</b> {avg_taxa:.2f}</p>
        <p style="font-size:16px;"><b>Max Nearby Taxa:</b> {max_taxa}</p>
        <p style="font-size:16px;"><b>Total Nearby Unique Taxa:</b> {unique_taxa}</p>
        <p style="font-size:16px;"><b>Total Nearby BINs:</b> {total_bin:,}</p>
        <p style="font-size:16px;"><b>Avg BIN per Site:</b> {avg_bin:.2f}</p>
        <p style="font-size:16px;"><b>Max Nearby BIN:</b> {max_bin}</p>
        <p style="font-size:16px;"><b>Total Nearby Unique BINs:</b> {unique_bin}</p>
    </div>
    """

    return pn.pane.HTML(html, width=300, height=270)

# Function to create a dashboard to visualize the data of mines and the surrounding taxa
def plot_points(results_df, title, total_sites, family_colour_map, species_colour_map, bin_colour_map, image_path):

    pn.extension()

    # Filter out sites that have no nearby taxa and remove duplicate coordinate entries 
    results_df = results_df[results_df['num_taxa_total'] > 0].copy()
    results_df = results_df.drop_duplicates(subset=["site_name", "mine_lat", "mine_lon"])

    # Convert list columns into readable comma-seperated strings for hover tooltipd
    results_df["species_str"] = results_df["species_unique"].apply(lambda x: ", ".join(x) if isinstance(x, list) else "")
    results_df["families_str"] = results_df["taxa_list_unique"].apply(lambda x: ", ".join(x) if isinstance(x, list) else "")

    # Sort sites by number of taxa (descending order) and assign an index label
    results_df = results_df.sort_values("num_taxa_total", ascending=False).reset_index(drop=True)
    results_df["x"] = [str(i + 1) for i in range(len(results_df))]

    # Assign colour codes by site type for visual distinction
    results_df["colour"] = results_df["site_type"].map({
        "MDMER Mines": "blue",
        "PPER Mills": "green",
        "EcoReg DNA Research Sites": "red"
    }).fillna("grey")

    # Prepare the Bokeh ColumnDataSource for interactive plotting
    x_vals = list(results_df["x"])
    source = ColumnDataSource(results_df)

    # Create bar chart showing number of taxa per site
    p1 = figure(
        x_range=x_vals,
        title=title,
        x_axis_label="Site Index",
        y_axis_label="Number of Taxa",
        width=700,
        height=500,
        tools="pan,wheel_zoom,reset,save,hover"
    )

    # Draw vertical bars representing taxa counts for each site 
    p1.vbar(x="x", top="num_taxa_total", width=0.8, source=source, color="colour", alpha=0.7)

    # Configure hover tooltips to display information
    hover = HoverTool()
    hover.attachment = "right"
    hover.tooltips = [
        ("Site", "@site_name"),
        ("# Taxa", "@num_taxa_total"),
        ("# Species", "@num_records_with_species_name"),
        ("# BINs", "@num_bins_total"),
        ("Families", "@families_str"),
        ("Species", "@species_str"),
        ("BINs", "@bin_ids_unique")
    ]
    p1.add_tools(hover)

    # Create a summary table of all sites with nearby taxa
    top_sites = results_df[["site_name", "metal/product", "num_taxa_total", "num_taxa_unique", "num_records_with_species_name", "num_species_unique", "num_bins_total", "num_bins_unique", "site_type", "mine_lat", "mine_lon"]]
    top_sites.index = range(1, len(top_sites) + 1)
    table = pn.pane.DataFrame(top_sites, height=400, width=700, name="Sites")

    # Generate pie charts summarizing family, species, and BIN
    family_pie_chart = make_pie_chart(results_df, "taxa_list_total", "Nearby Families", colour_map=family_colour_map)
    species_pie_chart = make_pie_chart(results_df, "species_total", "Nearby Species", colour_map=species_colour_map)
    bin_pie_chart = make_pie_chart(results_df, "bin_ids_total", "Nearby BINs", colour_map=bin_colour_map)

    # Put together the full dashboard layout 
    layout = pn.Column(
        pn.Row(
            p1,
            pn.Spacer(width=30),
            pn.Column("### All Sites with Nearby Taxa", table),
            pn.Spacer(width=30),
            pn.Column(make_summary(results_df, total_sites))
        ),
        pn.Spacer(height=40),
        pn.Row(
            pn.Column(family_pie_chart),
            pn.Spacer(width=5),
            pn.Column(species_pie_chart),
            pn.Spacer(width=5),
            pn.Column(bin_pie_chart),
        ),
        pn.Spacer(height=120),
        pn.Row(
            pn.Spacer(width=600),
            pn.Column(pn.pane.JPG(image_path, width=600, height=800))
        ),
        pn.Row(
            pn.Spacer(width=600),
            pn.Column(pn.pane.JPG('../../assets/Acknowledgments_Team_Funders.jpg', width=600))
        )
    )

    return layout

# Main Function
def main(taxa_dir, mine_dir, ecoreg_dir): 

    # Loads and preprocesses the taxa datasets
    taxa_data = preprocess_taxa(load_csvs(taxa_dir, "taxa"))

    # Loads all mine-related data
    mine_data = load_csvs(mine_dir, "mine")

    # Separates the mine data into Mines, PPER, and CS datasets
    mines, pper, cs, total_cs_sites, total_mine_sites, total_pper_sites = preprocess_mine(mine_data)

    # Runs the 5 km buffer analysis for each site type
    mine_results = quantify_points(taxa_data, mines, "MDMER Mines")
    pper_results = quantify_points(taxa_data, pper, "PPER Mills")
    cs_results = quantify_points(taxa_data, cs, "EcoReg DNA Research Sites")

    # Load and summarize Novel / EcoReg TSV data
    ecoreg_data = load_csvs(ecoreg_dir, "ecoreg")
    ecoreg_results = summarize_ecoreg(ecoreg_data)

    # Ensure one record per EcoReg site
    cs_results = cs_results.sort_values('num_taxa_total', ascending=False)
    cs_results = cs_results.drop_duplicates(subset=['site_name'], keep='first')

    # Merge all sites into one
    all_results = pd.concat([mine_results, pper_results, cs_results, ecoreg_results], ignore_index=True)

    # Extract all unique taxonomic identifiers across all sites for consistent colour assignment
    all_families = sorted(set([item for sublist in all_results['taxa_list_unique'] for item in sublist if isinstance(sublist, list)]))
    all_species = sorted(set([item for sublist in all_results['species_unique'] for item in sublist if isinstance(sublist, list)]))
    all_bins = sorted(set([item for sublist in all_results['bin_ids_unique'] for item in sublist if isinstance(sublist, list)]))

    # Generate consistent colour maps for all
    family_colour_map = generate_colour_map(all_families)
    species_colour_map = generate_colour_map(all_species)
    bin_colour_map = generate_colour_map(all_bins)

    # Combine into one dashboard
    tabs = pn.Tabs(
        ("MDMER Mines", plot_points(mine_results, "Nearby Taxa at MDMER Mines", total_mine_sites, family_colour_map, species_colour_map, bin_colour_map, image_path='../../assets/MDMER Mines in Canada.jpg')),
        ("PPER Mills", plot_points(pper_results, "Nearby Taxa at PPER Mills", total_pper_sites, family_colour_map, species_colour_map, bin_colour_map, image_path='../../assets/PPER Mills in Canada.jpg')),
        ("EcoReg DNA Research Sites", plot_points(cs_results, "Nearby Taxa at EcoReg DNA Research Sites", total_cs_sites, family_colour_map, species_colour_map, bin_colour_map, image_path='../../assets/EcoReg DNA Research Sites in Canada.jpg')),
        ("Novel Data to BOLD", plot_points(ecoreg_results, "Nearby Taxa at Novel Data to BOLD", 8, family_colour_map, species_colour_map, bin_colour_map, image_path='../../assets/Novel Data to BOLD in Canada.jpg'))
    )

    # Save and display dashboard
    tabs.save("taxa_dashboard_may_5_2026.html", embed=True)
    tabs.show()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 quantify_map.py <tax_directory> <mine_directory> <ecoreg_directory>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])

    