import os
import sys
import glob
import pandas as pd

# --- Site mapping from Field ID prefix ---
SITE_MAP = {
    "NR-NEA": "New Britannia",
    "NE-SWB": "New Britannia",
    "MIRL": "Flin Flon Lake",
    "SLMB": "Flin Flon Lake",
    "ROSL": "Flin Flon Lake",
    "SLNW": "Flin Flon Lake",
    "WR-FC": "Wesdome",
    "WR-ML": "Wesdome",
    "WE-MIC": "Wesdome",
    "AR-CB": "Anderson",
    "AR-PB": "Anderson",
    "AE-AB": "Anderson",
    "KR-AC": "Kirkland Lake",
    "KE-AC": "Kirkland Lake",
}

SITE_INFO = {
    "New Britannia": {"lat": 49.0, "lon": -95.0, "metal": "Gold"},
    "Flin Flon Lake": {"lat": 54.8, "lon": -101.9, "metal": "Zinc, Copper, Iron, Gold"},
    "Wesdome": {"lat": 48.0, "lon": -80.0, "metal": "Gold"},
    "Anderson": {"lat": 50.0, "lon": -85.0, "metal": "Gold, Copper, Zinc"},
    "Kirkland Lake": {"lat": 48.1, "lon": -79.9, "metal": "Gold"},
}

def load_csvs(directory):
    files = glob.glob(os.path.join(directory, '*.csv'))
    if not files:
        print(f"No CSV files found in {directory}")
        sys.exit(1)

    dfs = []
    for f in files:
        df = pd.read_csv(f, encoding='latin1', low_memory=False)
        df.columns = df.columns.str.strip().str.lower()
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)

def process_data(df):
    if "site code" not in df.columns:
        print("Error: 'field id' column not found.")
        sys.exit(1)

    df["site_prefix"] = df["site code"].str.extract(r"^([A-Z]+-[A-Z]+|[A-Z]+)")

    # Map to site name
    df["site_name"] = df["site_prefix"].map(SITE_MAP)

    # Assign coordinates + metal
    df["latitude"] = df["site_name"].map(lambda x: SITE_INFO.get(x, {}).get("lat"))
    df["longitude"] = df["site_name"].map(lambda x: SITE_INFO.get(x, {}).get("lon"))
    df["metal"] = df["site_name"].map(lambda x: SITE_INFO.get(x, {}).get("metal"))

    df.drop(columns=["site_prefix"], inplace=True)

    return df

def main(main_dir):
    df = load_csvs(main_dir)
    df = process_data(df)
    df.to_csv("NEW_ESRBI_BOLD.csv", index=False)
    print("Done. Output saved as NEW_ESRBI_BOLD.csv")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python combine_csvs.py <directory>")
        sys.exit(1)

    main(sys.argv[1])