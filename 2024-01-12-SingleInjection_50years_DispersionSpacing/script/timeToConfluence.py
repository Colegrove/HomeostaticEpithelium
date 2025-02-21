
## Hunter Colegrove
## 17 Jan 2024
## Plot time to confluence using HomeostaticEpithelium vis output files.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np


directory_path=f"/net/feder/vol1/project/HomeostaticEpithelium/dat/2024-01-12-SingleInjection_50years_DispersionSpacing/results/"

def confluence_or_loss(prefix, tConfluence):
    file_path = f"{directory_path}{prefix}.{tConfluence}.txt"
    #print(file_path)
    with open(file_path, 'rb') as file:
        try: 
            # Move to the end of the file
            file.seek(-2, os.SEEK_END)
            while file.read(1) != b'\n':
                file.seek(-2, os.SEEK_CUR)
        except OSError:
            file.seek(0)

        # Read the last line
        last_line = file.readline().decode()
        if(last_line != 'Confluence' and last_line != 'Loss'):
            last_line = 'Ongoing'

    return last_line
    

def find_max_suffix(directory):
    files = os.listdir(directory)
    file_suffixes = {}

    for file in files:
        # Split the filename into parts
        parts = file.split('.')
        # Exclude non-VisFiles
        if(parts[0].split("_")[0] != "VisFile"):
            continue
        
        # Extract the prefix and suffix
        prefix = '.'.join(parts[:-2])
        suffix = int(parts[-2])

        # Update the maximum suffix for each prefix and add params to dictionary
        params = prefix.split("_")
        block = float(params[2])
        sigma = float(params[6])
        dose = int(params[8])
        replicate = int(params[11])
        if prefix in file_suffixes:
            file_suffixes[prefix] = [block, sigma, dose, max(int(file_suffixes[prefix][3]), suffix), replicate]
        else:
            file_suffixes[prefix] = [block, sigma, dose, suffix, replicate]
        
    return file_suffixes


########
### Main
########

tConfluence_dict = find_max_suffix(directory_path)

## Add if simulation reached confluence or loss or neither
for prefix, params in tConfluence_dict.items():
    conf_loss = confluence_or_loss(prefix, params[3])
    params.append(conf_loss)
    tConfluence_dict[prefix] = params

df = pd.DataFrame.from_dict(tConfluence_dict, orient='index', columns=['Block_prob', 'Sigma', 'Dose', 'tConfluence', 'Replicate', 'Confluence'])
df = df.sort_values(by=['Block_prob', 'Sigma', 'Dose', 'Replicate'])
df['tConfluence'] = df['tConfluence']/(360*4.5)

df.to_csv('/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-12-SingleInjection_50years_DispersionSpacing/dat/timeToConfluence.csv')

# plot_df = df[df['Confluence'] == "Confluence"]
# plot_df = df[df['Confluence'] == "Loss"]
# print(plot_df)


# ########
# ### Plot
# ########

# sns.set_theme(style="ticks")

# # Initialize the figure
# f, ax = plt.subplots(figsize=(7, 6))

# # boxplot 
# sns.boxplot(plot_df, x="Block_prob", y="tConfluence", fill=False, whis=(0,100))

# # Add in points to show each observation
# sns.stripplot(plot_df, x="Block_prob", y="tConfluence")

# # Tweak the visual presentation
# ax.xaxis.grid(True)
# ax.set(ylabel="")
# sns.despine(trim=True, left=True)
# yticks = np.arange(0, 51, 5)
# ax.set_yticks(yticks)

# plt.xlabel("Blocking Probability")
# #plt.ylabel("Years to Confluence - (80% Tissue Coverage)")
# plt.ylabel("Years to Loss")
# plt.savefig(f"/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-01-05-SingleInjection_50years_blockingProb/results/tConfluence.png")