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
    inc_group.append(mob_df['rep_type(s)'][0])  # add incompatibility group to the list
    print(f'Added {mob_df["rep_type(s)"][0]} to the incompatibility group list')

# create a dataframe with all of the extracted metadata
metadata_df = pd.DataFrame({'sample_id': sample_id, 'mobility': mobility, 'relaxase': relaxase,
                            'oriT': oriT, 'gc': gc, 'size': size, 'inc_group': inc_group})

try:
    os.mkdir(dir_out)  # create the output directory if it doesn't already exist
except:
    pass

metadata_out = f'{dir_out}/metadata.tsv'  # output file with plasmid metadata

metadata_df.to_csv(metadata_out, sep='\t', index=False)  # write the metadata to a tab-delimited file

# create a dataframe with incompatibility group information for each plasmid
inc_list = []  # list of incompatibility groups
print(f'Creating incompatibility group dataframe for {len(inc_group)} plasmids')
for _inc in inc_group:  # loop through each incompatibility group in the list
    print(f'Processing incompatibility group {_inc}')
    for inc in _inc.split(','):  # loop through each incompatibility group in the comma-separated string
        if inc in inc_list:  # if the incompatibility group is already in the list,
            print(f'\t{inc} already in the list')  # print a message
            pass  # and do nothing
        else:  # otherwise,
            inc_list.append(inc)  # add the incompatibility group to the list
            print(f'\tAdded {inc} to the list')  # print a message

column_headers = ['Plasmid'] + inc_list  # create a list of column headers
print(f'Column headers: {column_headers}')  # print the list of column headers

inc_group_df = pd.DataFrame(columns=column_headers)  # create a dataframe with the column headers
inc_group_df['Plasmid'] = plasmid_list  # add a column for the plasmid name
print(f'Inc group dataframe:\n{inc_group_df}')  # print the dataframe

n = 0  # initialize a counter
for plasmid in inc_group_df['Plasmid']:  # loop through each plasmid in the dataframe
    print(f'Processing plasmid {plasmid}')
    tsv_in = pd.read_csv(f'mobsuite/{plasmid}.tsv', sep='\t')  # read in the MobSuite output file for that plasmid
    rep_type = str(tsv_in['rep_type(s)'][0])  # extract the incompatibility group(s) from the file
    rep_list = rep_type.split(',')  # split the incompatibility group(s) into a list
    for inc in inc_list:  # loop through each incompatibility group in the list
        print(f'\tProcessing incompatibility group {inc} for plasmid {plasmid}')
        if inc in rep_list:  # if the incompatibility group is in the list of incompatibility groups for the plasmid,
            inc_group_df[inc][n] = 1  # set the value in the dataframe to 1
            print(f'\t\tSetting value to 1')
        else:  # otherwise,
            inc_group_df[inc][n] = 0  # set the value in the dataframe to 0
            print(f'\t\tSetting value to 0')
    n += 1  # increment the counter

plasmid_inc_out = f'{dir_out}/plasmid_inc.tsv'  # output file with plasmid incompatibility group information

print(f'Writing incompatibility group dataframe to {plasmid_inc_out}')
inc_group_df.to_csv(plasmid_inc_out,sep='\t', index=False)  # write the dataframe to the output file