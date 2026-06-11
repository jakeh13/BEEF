import os

def make_output_folders(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/layout", exist_ok=True)
    os.makedirs(f"{output_dir}/layout/cross_sections", exist_ok=True)
    os.makedirs(f"{output_dir}/layout/monitors", exist_ok=True)
    os.makedirs(f"{output_dir}/data", exist_ok=True)
    os.makedirs(f"{output_dir}/logs", exist_ok=True)
    os.makedirs(f"{output_dir}/fields", exist_ok=True)
    