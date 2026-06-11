import os

def resolve_output_folder(base: str) -> str:
    """Return base if it does not exist, otherwise base_00, base_01, …"""
    if not os.path.exists(base):
        return base
    for i in range(100):
        candidate = f"{base}_{i:02d}"
        if not os.path.exists(candidate):
            return candidate
    raise RuntimeError(f"Could not find a free output folder after 100 attempts for base '{base}'")

def make_output_folders(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/layout", exist_ok=True)
    os.makedirs(f"{output_dir}/layout/cross_sections", exist_ok=True)
    os.makedirs(f"{output_dir}/layout/monitors", exist_ok=True)
    os.makedirs(f"{output_dir}/data", exist_ok=True)
    os.makedirs(f"{output_dir}/logs", exist_ok=True)
    os.makedirs(f"{output_dir}/fields", exist_ok=True)
    