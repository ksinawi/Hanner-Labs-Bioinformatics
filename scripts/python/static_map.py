# Created By: Faisal Isam
# Last Updated: October 9, 2025
# Input: taxa directory and mine directory
# Required File: "ca.json" 

import os
import glob
import sys

import pandas as pd                     # For working with data
import geopandas as gpd                 # For geospatial data manipulation
from shapely.geometry import Point      # To create geometry objects from lat/lon

import cartopy.crs as ccrs              # Coordinate Reference System
import cartopy.feature as cfeature      # For adding land, rivers, etc
import matplotlib.pyplot as plt         # For creating the plot 
import matplotlib.patches as mpatches   # For creating the scale 

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

    dfs = []

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

    # Remove rows with missing or invalid coordinates or family names
    df = df.dropna(subset=["lat", "lon", "family_name"])
    df = df[df["lat"] != "N/A"]
    df = df[df["lon"] != "N/A"]

    # Convert coordinates to numeric values, coercing errors to NaN
    df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
    df["lon"] = pd.to_numeric(df["lon"], errors="coerce")

    # Drop rows that still have missing coordinate values
    df = df.dropna(subset=["lat", "lon"])

    # Remove duplicate sampling points (same family and location)
    df = df.drop_duplicates(subset=["family_name", "lat", "lon"])

    return df

# Function to process the mine data frames
def preprocess_mine(df):

    mine_coords = []
    pper_coords = []
    cs_coords = []

    # Iterate over each row and check for coordinate formats used by each data type
    for _, row in df.iterrows():

        # EcoReg DNA Research Sites 
        lat, lon = row.get("station_latitude"), row.get("station_longitude")    # Get values        
        if pd.notna(lat) and pd.notna(lon):                                     # Check for valid values
            metal = row.get("metal")                                            # Get corresponding metal
            if pd.notna(metal):                                                 # Check for valid metal
                # Convert to proper variable type and add to list 
                try:
                    cs_coords.append({"lat": float(lat), "lon": float(lon), "metal": str(metal)}) 
                except:
                    pass

        # MDMER Mines
        # Same steps as DNA Research Sites
        lat, lon = row.get("latitude"), row.get("longitude")
        if pd.notna(lat) and pd.notna(lon):
            metal = row.get("metal")
            if pd.notna(metal):
                try:
                    mine_coords.append({"lat": float(lat), "lon": float(lon), "metal": str(metal)})
                except:
                    pass

        # PPER Mils (coordinates stored in one text field)
        coords = row.get("coordinates")                                      # Get coordinates field
        if pd.notna(coords) and coords != "N/A":                             # Check for valid coordinates
            product = row.get("product")                                     # Get corresponding PPER
            if pd.notna(product):                                            # Check for valid PPER
                # Splits field by the comma to get the lat and lon 
                try:
                    parts = coords.split(",")
                    if len(parts) == 2:
                        lat_f, lon_f = float(parts[0].strip()), float(parts[1].strip())
                        pper_coords.append({"lat": lat_f, "lon": lon_f, "product": str(product)})
                except:
                    pass

    # Converts the lists into DataFrames
    return pd.DataFrame(mine_coords), pd.DataFrame(pper_coords), pd.DataFrame(cs_coords)

# Function to filter points to only those within the Canadian boundaries
def filter_points(df, lat_col="lat", lon_col="lon", country_geo="ca.json"):

    # Load Canada GeoJSON and ensure WGS84 (lat/lon) coordinate system
    canada_json = gpd.read_file(country_geo).to_crs(epsg=4326)
    df = df.dropna(subset=[lat_col, lon_col])

    # Create a GeoDateFrame from the latitude and longitude
    gdf = gpd.GeoDataFrame(
        df,
        geometry=[Point(xy) for xy in zip(df[lon_col], df[lat_col])],
        crs="EPSG:4326"
    )

    # Keep only points that fall within the Canadian boundary
    gdf = gdf[gdf.within(canada_json.union_all())]

    return pd.DataFrame(gdf.drop(columns="geometry"))

# Function to draw a north arrow on the map
def add_north_arrow(ax, x=0.05, y=0.8, size=20):

    ax.text(x, y+0.07, 'N', transform=ax.transAxes,
        ha='center', va='center', fontsize=size, fontweight='bold')

    ax.annotate('', xy=(x, y+0.05), xytext=(x, y-0.02),
    arrowprops=dict(facecolor='black', width=5, headwidth=15),
        xycoords=ax.transAxes)

# Function to add a scale bar
def add_scale_bar(ax, length_km=500, location=(0.05, 0.05)):

    x, y = location
    bar_height = 0.01

    # Draw alternating black/white rectanlges to represent the scale bar
    for i in range(5):
        ax.add_patch(mpatches.Rectangle(
            (x + i*0.05, y), 0.05, bar_height,
            transform=ax.transAxes,
            facecolor='black' if i % 2 == 0 else 'white',
            edgecolor='black'
        ))

    # Add numeric labels for 0 and the total distance
    ax.text(x, y - 0.025, "0", transform=ax.transAxes,
        ha="center", va="top", fontsize=9)
        
    ax.text(x + 0.25, y - 0.025, f"{length_km} km",
        transform=ax.transAxes,
        ha="center", va="top", fontsize=9)

# Function to plot all datasets onto the static map
def plot_map(taxa_df, mines_df, cs_df, pper_df, geojson_file="ca.json"):

    # Load Canadian boundary geometry
    canada = gpd.read_file(geojson_file)

    # Generate province centroids for labeling
    canada_proj = canada.to_crs(epsg=3347)
    canada_proj['centroid'] = canada_proj.geometry.centroid
    canada_proj['lon'] = canada_proj.centroid.x
    canada_proj['lat'] = canada_proj.centroid.y
    canada_proj = canada_proj.to_crs(epsg=4326)

    # Create the figure and map axis
    fig = plt.figure(figsize=(20, 18))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Zoom to Canada bounding box
    minx, miny, maxx, maxy = canada.total_bounds
    ax.set_extent([minx-2, maxx+2, miny-2, maxy+2], crs=ccrs.PlateCarree())

    # Draw base map features
    ax.add_feature(cfeature.LAND, facecolor='white')
    ax.add_feature(cfeature.OCEAN, facecolor='lightgrey')
    ax.add_feature(cfeature.LAKES, alpha=0.4, edgecolor='lightgrey', facecolor='lightgrey')
    ax.add_feature(cfeature.BORDERS, edgecolor='lightgrey')
    ax.add_feature(cfeature.RIVERS, linewidth=0.5, edgecolor='lightgrey')
    ax.add_feature(cfeature.COASTLINE, edgecolor='lightgrey')

    # Add province borders and labels 
    canada.boundary.plot(ax=ax, edgecolor='lightgrey', linewidth=1, transform=ccrs.PlateCarree())
    for idx, row in canada.iterrows():
        centroid = row['geometry'].centroid 
        code = row['id'][-2:] # Last two characters represent the province symbol
        ax.text(centroid.x, centroid.y, code,
            fontsize=7, color='grey',
            ha='center', va='center',
            transform=ccrs.PlateCarree())

    # Plot Taxa
    if not taxa_df.empty:
        ax.scatter(taxa_df['lon'], taxa_df['lat'], s=25, c='#B0D609', marker='o', alpha=0.8,
        label='Important EEM Taxa in BOLD', transform=ccrs.PlateCarree())

    # Plot Mines
    if not mines_df.empty:
        ax.scatter(mines_df['lon'], mines_df['lat'], s=40, c='#486AD6', marker='s', alpha=0.8,
        label='MDMER Mines', transform=ccrs.PlateCarree())

    # Plot PPER Mills
    if not pper_df.empty:
        ax.scatter(pper_df['lon'], pper_df['lat'], s=40, c='#44ADD6', marker='D', alpha=0.8,
        label='PPER Mills', transform=ccrs.PlateCarree())

    # Plot EcoReg Sites
    if not cs_df.empty:
        ax.scatter(cs_df['lon'], cs_df['lat'], s=40, c='#D6581E', marker='^', alpha=0.8,
        label='EcoReg DNA Research Sites', transform=ccrs.PlateCarree())

    # Add legend, compass, scale, and title
    add_north_arrow(ax, x=0.07, y=0.88)
    add_scale_bar(ax, length_km=2000, location=(0.05, 0.05))
    ax.legend(loc='upper right', fontsize=8)
    plt.title("EEM Benthic Research Sampling", fontsize=20)
    plt.savefig("full_static_map_nov_13_2025.pdf")
    plt.show()

# Main Function
def main(taxa_dir, mine_dir):

    # Load taxa and mine data sets from provided directories
    taxa_data = preprocess_taxa(load_csvs(taxa_dir, "taxa"))
    mine_data = load_csvs(mine_dir, "mine")

    # Extract differernt subsets (mines, mills, cs)
    mines, pper, cs = preprocess_mine(mine_data)
     
    # Keep only taxa points located within Canada
    taxa_data = filter_points(taxa_data, "lat", "lon", "ca.json")

    # Create and display the final static map
    plot_map(taxa_data, mines, cs, pper, geojson_file="ca.json")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python static_map.py <tax_directory> <mine_directory>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])