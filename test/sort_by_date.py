# import re

# # List of names
# names = [
#     "SENTINEL2A_20210312-194424-871_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20210605-194425-892_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20211107-194428-015_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20210113-074924-099_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20210707-074924-705_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20211207-194425-938_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20210209-073926-725_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20210807-074426-679_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20211228-193926-013_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20210405-073924-400_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20210907-073930-033_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20210507-193925-576_L3A_T38SQG_C_V1-0",
#     "SENTINEL2X_20211008-194429-848_L3A_T38SQG_C_V1-0"
# ]

# # Define a key function to extract the date portion
# def extract_date(name):
#     match = re.search(r'_(\d{8})-', name)
#     if match:
#         return match.group(1)
#     return ""

# # # Sort the names based on the extracted date
# # sorted_names = sorted(names, key=extract_date)

# # Using a lambda function within the sorted function
# sorted_names = sorted(names, key=lambda name: re.search(r'_(\d{8})-', name).group(1) if re.search(r'_(\d{8})-', name) else "")

# # Print the sorted names
# for name in sorted_names:
#     print(name)


import re

original_string = "SENTINEL2A_20210312-194424-871_L3A_T38SQG_C_V1-0_FRC_B3.tif"
new_string = re.sub(r'B\d+', 'NRGB', original_string)

print("Original String:", original_string)
print("Modified String:", new_string)