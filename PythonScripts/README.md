# Python scripts for transfer function comparison


## ```TFCompare_CosmicIC_CLASS.py```
This script compares and visualizes normalized data from two different sources (e.g., `cosmicic` and `class` simulation outputs). It generates log-log plots for the transfer function.

### Usage
The script requires two input files via command-line arguments. Each file is expected to be a text-based format with at least 3 columns.
```
python3 plot_comparison.py --file_cosmicic path/to/cosmicic_data.txt --file_class path/to/class_data.txt
```

### Arguments
* ```--file_cosmicic```: Path to the first data file (rendered with circular markers).
* ```--file_class```: Path to the second data file (rendered as a solid line).

## Outputs
Upon successful execution, the script saves two images in the current directory:
1.  ```fig1_col2.png```: Comparison of the normalized second column.
2.  ```fig2_col3.png```: Comparison of the normalized third column.
