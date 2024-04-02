# Summarises the output of MobSuite for each plasmid.
# Outputs a tsv file with plasmid name, mobility, relaxase type, orit type, GC content, size, and inclusion group.
# The output file is plasmid_metadata.tsv in the plasmid_grouping directory.
#
# The input files are located in the mobsuite directory.
# The input files should be named <plasmid>.tsv where <plasmid> is the name of the plasmid.

import pandas as pd
import os

# Read in all of the MobSuite output files and extract metadata for each plasmid

mobsuite_dir = 'mobsuite'  # location of MobSuite output files
mobsuite_in_list = []  # list of MobSuite output files
for plasmid_tsv in os.listdir(mobsuite_dir):
    mobsuite_in_list.append(f'{mobsuite_dir}/{plasmid_tsv}')
    print(f'Added {mobsuite_dir}/{plasmid_tsv} to the list')

plasmid_list = []  # list of plasmid names extracted from the MobSuite file names
for _plasmid in mobsuite_in_list:
    plasmid_list.append(_plasmid.replace('mobsuite/', '').replace('.tsv', ''))
    print(f'Added {_plasmid.replace("mobsuite/", "").replace(".tsv", "")} to the list')

dir_out = 'plasmid_grouping'  # output directory for plasmid metadata files

plasmid_inc_out = f'{dir_out}/plasmid_inc.tsv'  # output file with plasmid inclusion group information
metadata_out = f'{dir_out}/metadata.tsv'  # output file with plasmid metadata

# initialize lists to store extracted metadata for each plasmid
sample_id = []
mobility = []
relaxase = []
oriT = []
gc = []
size = []
inc_group = []

# iterate through all of the MobSuite output files
for mobsuite in mobsuite_in_list:
    print(f'Reading {mobsuite}')
    mob_df = pd.read_csv(mobsuite, sep='\t')  # read in the MobSuite output file
    sample_id.append(mob_df['sample_id'][0])  # add the sample ID to the list
    print(f'Added {mob_df["sample_id"][0]} to the sample ID list')
    mobility.append(mob_df['predicted_mobility'][0])  # add mobility to the list
    print(f'Added {mob_df["predicted_mobility"][0]} to the mobility list')
    relaxase.append(mob_df['relaxase_type(s)'][0])  # add relaxase type to the list
    print(f'Added {mob_df["relaxase_type(s)"][0]} to the relaxase type list')
    oriT.append(mob_df['orit_type(s)'][0])  # add orit type to the list
    print(f'Added {mob_df["orit_type(s)"][0]} to the orit type list')
    gc.append(mob_df['gc'][0])  # add GC content to the list
    print(f'Added {mob_df["gc"][0]} to the GC content list')
    size.append(mob_df['size'][0])  # add plasmid size to the list
    print(f'Added {mob_df["size"][0]} to the plasmid size list')
    inc_group.append(mob_df['rep_type(s)'][0])  # add inclusion group to the list
    print(f'Added {mob_df["rep_type(s)"][0]} to the inclusion group list')

# create a dataframe with all of the extracted metadata
metadata_df = pd.DataFrame({'sample_id': sample_id, 'mobility': mobility, 'relaxase': relaxase,
                            'oriT': oriT, 'gc': gc, 'size': size, 'inc_group': inc_group})

try:
    os.mkdir(dir_out)  # create the output directory if it doesn't already exist
except:
    pass

metadata_df.to_csv(metadata_out, sep='\t', index=False)  # write the metadata to a tab-delimited file

