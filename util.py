TARGET_TYPE_MAPPING = {
    0: "Flat",
    1: "Cylindrical",
    2: "Elliptical",
}

def get_target_type_name(target_type: int) -> str:
    return TARGET_TYPE_MAPPING.get(target_type, "Unknown")