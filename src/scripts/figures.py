import matplotlib.pyplot as plt
import pandas as pd
from constant import *
import seaborn as sns
import numpy as np

sns.set()
sns.color_palette("bright")


def main():
    marker_symbols = ['o', 'v', 's', 'p', 'X', '*', 'D', '<', '>']
    markers = dict(zip(ANCHOR_TYPES, marker_symbols))
    chain_results = {}
    anchor_results = {}
    for _, _, files in os.walk(RESULT_FOLDER):
        for fi in files:
            print(f'{RESULT_FOLDER}{fi}')
            file_name_parts = fi.split('-')
            try:
                k = file_name_parts[-2]
                genome = file_name_parts[-1]
            except:
                continue
            key = f'k: {k}'
            if 'anchor' in fi:
                # is an anchor summary
                anchor_df = pd.read_csv(f'{RESULT_FOLDER}{fi}')
                anchor_results[key] = anchor_df
            if 'chain' in fi:
                chain_df = pd.read_csv(f'{RESULT_FOLDER}{fi}')
                chain_results[key] = chain_df
    x = []

    for name in chain_results.keys():
        x.append(name)
    
    x_values = np.arange(len(x))
  
    for mode in ['total', 'tidy']:
        for result_type in chain_df.columns[3:]:
            jitter = -0.2
            for anchor_type in ANCHOR_TYPES:
                plt.scatter([val_x + jitter  for val_x in x_values], 
                            [df.groupby(['mode', 'type']).mean()[
                            result_type][mode][anchor_type] for name, df in chain_results.items()],
                            label=anchor_type,
                            marker=markers[anchor_type])
                jitter += 0.4/7
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.xticks(x_values, x)
            plt.savefig(
                f'{RESULT_FOLDER}/chain-{result_type}-{mode}.svg', bbox_inches='tight')
            plt.clf()

    for mode in ['total', 'tidy']:
        for result_type in anchor_df.columns[2:]:
            jitter = -0.2
            for anchor_type in ANCHOR_TYPES:
                plt.scatter([val_x + jitter  for val_x in x_values], [df.groupby(['mode', 'type']).mean()[
                            result_type][mode][anchor_type] for name, df in anchor_results.items()],
                            label=anchor_type.replace("-", " "),
                            marker=markers[anchor_type])
                jitter += 0.4/7
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.xticks(x_values, x)
            plt.savefig(
                f'{RESULT_FOLDER}/anchor-{result_type}-{mode}.svg', bbox_inches='tight')
            plt.clf()


    reads = get_reads(READS_DIR)
    total_read_lengths = 0
    total_number_of_errors = 0
    total_error_prob = 0
    for i in reads.keys():
        read_length, _, _, number_of_errors, total_error_probability = get_read_properties(
            READ_PATH.format(i))
        total_read_lengths += read_length
        total_number_of_errors += number_of_errors
        total_error_prob += total_error_probability
    print(f'average read length: {total_read_lengths/len(reads)}')
    print(f'average number of errors: {total_number_of_errors/ len(reads)}')
    print(f'average error probability: {total_error_prob/len(reads)}')


if __name__ == '__main__':

    main()
