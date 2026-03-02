import os
import json
import pandas as pd
from core.simulation import SimulationManager

def main():
    print("=== ALP Propagation Simulation Pipeline ===")
    
    config_path = "config.json"
    if not os.path.exists(config_path):
        print(f"[ERROR] Configuration file not found at {config_path}")
        return

    with open(config_path, 'r') as f:
        config = json.load(f)

    # Ensure output directories exist
    os.makedirs("results", exist_ok=True)
    os.makedirs("data", exist_ok=True)

    catalog_path = config["file_paths"]["source_catalog_output"]
    if not os.path.exists(catalog_path):
        print(f"[ERROR] Source catalog not found. Please run utils/generate_source_catalog.py first.")
        return

    sources_df = pd.read_csv(catalog_path)
    print(f"[INFO] Loaded {len(sources_df)} sources from catalog.")

    # Initialize and run simulation
    manager = SimulationManager(config)
    results_df = manager.run_full_scan(sources_df)

    output_path = config["file_paths"]["final_events_output"]
    results_df.to_csv(output_path, index=False)
    print(f"[INFO] Simulation complete. Results saved to: {output_path}")

if __name__ == "__main__":
    main()