
## Hunter Colegrove
## 19 Sep 2023
## Plot colony size distributions using HomeostaticEpithelium vis output files.
## Use to plot results from a multiple simulation experiment at multiple timepoints.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


XMAX_BLOCK = 1000 # set x-axis maximum value (this will depend on the data)
#XMAX_GROWTH = 150
XTICK_BLOCK = 100 # set x-axis ticks interval (this will depend on the data)
#XTICK_GROWTH = 15
#COLOR = "#AA4499" # (teal) set the color scheme for the plot. opacity decremented based on time points
#COLOR = "#332288" # (indigo)

# Read in vis file and convert to dataframe
def read_vis_data(vis_path):
    try:
        df = pd.read_csv(vis_path, delimiter="\t", header=None)
        df = df.rename(columns={0: 'x', 1: 'z', 2: 'y', 3: 'h', 4: 's', 5: 'v', 6: 'alpha', 7: 'mut', 8: 'injSite'})
        return df
    except:
        print(f"File: {vis_path} not found.")

# Look only in the basal layer
def basal_only(df):
    df = df[df.y < 1]
    return df

# Remove x,z,y coordinates from dataframe
def remove_coords(df):
    df = df.drop(columns=['x', 'z', 'y'])
    return df

# Remove Reference cells and dead cells from the colony counts
def remove_ref_dead(df):
    df = df[(df['mut'] != 0.0) & (df['mut'] != -1)]
    df = df[(df['injSite'] >= 0)]
    # df = df[(df['h'] != 1.0) & (df['s'] != 1.0) & (df['v'] != 1.0)]
    # df = df[(df['alpha'] != 0.0)]
    return df

# Compute the amount of colonies greater or equal to a given x-axis value
def compute_colonies(x_axis, colony_sizes):
    y_values = [sum(value >=x for value in colony_sizes) for x in x_axis]
    return y_values

# Compute the average colony size across all replicates
def compute_average_colonies(all_y_values, replicates):
    ave_colony_counts = [sum(y_values) / replicates for y_values in zip(*all_y_values)]
    return ave_colony_counts

# Describes the opacity pallette for the plot color scheme
def opacity_palette_timeseries(num_timesteps):
    opacity_levels = [(i+1) / num_timesteps for i in range(num_timesteps)]
    #print(opacity_levels)
    palette = [alpha for alpha in opacity_levels]
    return palette

# Describes the opacity palette for the plot color scheme
def opacity_palette_parameter(num_params):
    opacity_levels = [(i+1) / num_params for i in range(num_params)]
    #print(opacity_levels)
    palette = [alpha for alpha in opacity_levels]
    return palette
 
# # Plot colony size distribution - all parameters at final timepoint
# def plot_parameter_block(block_values, colony_sizes, key):
#     x_axis_seq = range(1,XMAX_BLOCK)
#     x_axis_ticks = [i if i == 1 else i-1 for i in range(1, XMAX_BLOCK+XTICK_BLOCK, XTICK_BLOCK)]
#     y_axis_seq = range(0,11,1) # based on the number of cells correted
#     num_params = len(block_values)
#     #block_values = block_values[::-1]
#     #palette = opacity_palette_parameter(num_params)
#     palette = ['#332288', '#DDCC77', '#44AA99', '#AA4499', '#88CCEE', '#CC6677', "#999933"]

#     sns.set_style("whitegrid")
#     plt.figure(figsize=(10, 6))
#     sns.set(font_scale = 1.6)
#     sns.set_style(rc = {'axes.facecolor': 'DDDDDD'})
#     # for i, parameter, y_values in zip(palette, block_values, colony_sizes):
#     #     sns.lineplot(x=x_axis_seq, y=y_values, marker='o', label=f"Blocking Prob. = {parameter}", color=COLOR, alpha=i)
#     print(f"\n{len(palette)} {len(block_values)} {len(colony_sizes)}\n")
#     for i, parameter, y_values in zip(palette, block_values, colony_sizes):
#         print(f"x_axis_seq: {len(x_axis_seq)}, y_values: {len(y_values)}\n")
#         sns.lineplot(x=x_axis_seq, y=y_values, marker = 'o', label=f"Blocking Prob. = {parameter}", color=i)
#     print("check")
#     # Reverse the order of the legend
#     handles, labels = plt.gca().get_legend().legend_handles, plt.gca().get_legend().texts
#     handles = handles[::-1]
#     labels = [text.get_text() for text in labels][::-1]  # Extract text and reverse the order
#     plt.legend(handles, labels)

#     plt.xticks(x_axis_ticks)
#     plt.yticks(y_axis_seq)
#     plt.xlabel("clone size (X)")
#     plt.ylabel("number of clones >= X")
#     #plt.title(f"Colony size distributions of 10 corrected cells in one injection site.\nTime {key/4.5} Days")
#     plt.title(f"Colony size distributions of 10 corrected cells in one injection site.\nTime 1 Year")
#     #plt.legend()
#     plt.savefig(f"/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-04-10-SingleInjection_Confluence_blockProb/resultsColonySizeDist_block_{key}.png")


######
## Main
######

years = 50
timesteps = [2] # make a list of the model vis output timepoints

for i in range(1,years*2):
    time = int(i * 182*4.5)
    timesteps.append(time)

all_y_values = [] # values to be stored during the loop
block_values = [0, 0.001, 0.01, 0.1, 0.2, 0.5, 1] # values from experiment taken from snakefile
sigma = 2
dose = 10
growth_val = 1.0
replicates = [x for x in range(0,100)] # num replicates from experiment taken from snakefile
plot_list = [] # list containing average across replicates to be used for plotting 
plot_dict = {} # contains the plot_list as value and the selective advantage parameter as the key
cloneSize_dict = {}


for block_val in block_values:
    for timestep in timesteps:
        colony_sizes = []
        for replicate in replicates:
            vis_path=f"/net/feder/vol1/project/HomeostaticEpithelium/dat/2024-04-10-SingleInjection_Confluence_blockProb/results/VisFile_block_{block_val}_growth_1.0_sigma_{sigma}_dose_{dose}_correctionTime_1_{replicate}.{timestep}.txt"
            ##### Data handling functions
            df = read_vis_data(vis_path)
            if df is None:
                continue
            basal_df = basal_only(df)
            df_hsv = remove_coords(basal_df)
            df_colony = remove_ref_dead(df_hsv)
            colony_size = df_colony.value_counts().to_list() # generate colony counts as a list
            sum = df_colony.value_counts().sum()
            colony_sizes.append(sum)
        print(f'length of colony sizes: {len(colony_sizes)}')
        key = f'{block_val}_{timestep}'   
        if(len(colony_sizes) > 0): 
            cloneSize_dict[key] = [colony_sizes, block_val, timestep, sigma, dose]

        #y_values_ave = compute_average_colonies(all_y_values, len(replicates)) # average colony counts across replicates
        #plot_list.append(y_values_ave)
        #y_values = compute_colonies(range(1,XMAX_BLOCK), colony_sizes) # sum colony counts for y values
        #all_y_values.append(y_values) # add colony count to list
        #all_y_values = []
        #print(f'{growth_val}.{timestep} loop completed')
        #print(f'{block_val}.{timestep} loop completed')
    #plot_dict[growth_val] = plot_list
    #plot_dict[block_val] = plot_list
    #print(plot_dict)
    #plot_list = []
print("cloneSize dict: ") 
print(cloneSize_dict)                 

df = pd.DataFrame.from_dict(cloneSize_dict, orient='index', columns=['cloneSizes', 'blockProb', 'year', 'sigma', 'dose'])
df['year'] = df['year']/(364*4.5)

df.to_csv('/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/2024-04-10-SingleInjection_Confluence_blockProb/dat/cloneSizes.csv')



##### Plotting colony size distributions of each parameter across collected timepoints
#for key in plot_dict:
#    plot_timeseries(timesteps, plot_dict[key], key)

##### Plotting colony size distributions of all parameters at a single timepoint
# plot_dict_params = {} # new dictionary that holds as a key the final timepoint and as a value a list of lists of all parameter average simulations
# plot_list_params = []
# final_t_position = len(timesteps) - 1
# final_t = timesteps[final_t_position]
# for key in plot_dict:
#     plot_list_params.append(plot_dict[key][final_t_position])
# plot_dict_params[final_t] = plot_list_params

# print(plot_dict_params)
#print(plot_dict_params[final_t])
#plot_parameter_block(block_values, plot_dict_params[final_t], final_t)
#plot_parameter_growth(growth_values, plot_dict_params[final_t], final_t)


