"""
This script generates a LaTeX table which presents the numerical scores
comparing multizone simulations to APOGEE data.
"""

import pandas as pd
import paths
from score_multizone_outputs import main as gen_scores
from score_multizone_outputs import DTD_LIST

def main():
    # Import numerical scores of multi-zone outputs, or generate if missing
    csv_file = paths.output / 'scores.csv'
    if not csv_file.exists():
        gen_scores()
    df = pd.read_csv(csv_file, index_col=[0, 1])
    # Drop bimodality column (no numerical scores)
    df.drop('bimodality', axis=1, inplace=True)
    # Round scores
    df = df.round(3)
    for col in df.columns:
        df[col] = df[col].astype('str')
    # Convert to LaTeX format and write to output file
    # Fancy row labels
    dtd_labels = ['Two-population', '($t_p=0.05$ Gyr)', '', '',
                  'Power-law', '($\\alpha=-1.4$)', '', '',
                  'Power-law', '($\\alpha=-1.1$)', '', '',
                  'Exponential', '($\\tau=1.5$ Gyr)', '', '',
                  'Exponential', '($\\tau=3.0$ Gyr)', '', '',
                  'Plateau', '($W=0.3$ Gyr)', '', '',
                  'Plateau', '($W=1.0$ Gyr)', '', '',
                  'Triple-system', '($t_{\\rm rise}=0.5$ Gyr)', '', '']
    sfh_labels = ['Inside-out', 'Late-burst', 'Early-burst', 'Two-infall']
    df.reset_index(drop=False, inplace=True)
    df['DTD'] = dtd_labels
    df['SFH'] = sfh_labels * len(DTD_LIST)
    latex_table = df.style.hide(axis=0).to_latex()
    # Remove tabular environment & add horizontal lines
    rows = latex_table.split('\n')[2:-2]
    for i in range(3, len(rows)-4, 4):
        rows[i] = rows[i].replace('\\\\', '\\\\ \n\\hline')
    latex_table = '\n'.join(rows)
    # Import table header and footer
    with open(paths.scripts / 'scores_table_header.txt', 'r') as f:
        header_footer = f.read()
        header, footer = header_footer.split('===')
    # Replace tabular environment with deluxetable
    latex_table = header + latex_table + footer
    with open(paths.output / 'scores_table.tex', 'w') as f:
        f.write(latex_table)


if __name__ == '__main__':
    main()
