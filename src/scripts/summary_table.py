"""
This script generates a LaTeX table which summarizes the results of the 
multizone simulations across multiple parameter spaces.
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
    # Convert to LaTeX format and write to output file
    # Mask numerical scores with no / meh / yes marks
    for col in df.columns:
        if col != 'bimodality':
            df[col] = df[col].mask(df[col] < df[col].quantile(0.33), 
                               other='\yes').mask(
                                   (df[col] >= df[col].quantile(0.33)) & 
                                   (df[col] < df[col].quantile(0.67)), 
                               other='\meh').mask(
                                   df[col] >= df[col].quantile(0.67),
                               other='\\no')
    # Mask bimodality booleans with no / yes marks
    df['bimodality'] = df['bimodality'].mask(df['bimodality'].astype('bool'), 
                                        other='\yes').mask(
                                             ~df['bimodality'].astype('bool'),
                                        other='\\no')
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
    # New column labels
    df.columns = ['DTD', 'SFH', 'MDF', '[O/Fe] DF', '[O/Fe] Bimodality', 
                  '[Fe/H]--[O/Fe]', 'Age--[O/Fe]']
    latex_table = df.to_latex(column_format='ll|ccccc', index=False)
    # Replace \toprule, \midrule, \bottomrule with \hline
    latex_table = latex_table.replace('\\toprule', '\\hline\\hline')
    latex_table = latex_table.replace('\\midrule', '\\hline')
    latex_table = latex_table.replace('\\bottomrule', '\\hline')
    # Add horizontal lines in between DTD sections
    rows = latex_table.split('\n')
    for i in range(3, len(rows)-4, 4):
        rows[i] = rows[i].replace('\\\\', '\\\\ \n\\hline')
    latex_table = '\n'.join(rows)
    # Write to output
    with open(paths.output / 'summary_table.tex', 'w') as f:
        f.write(latex_table)


if __name__ == '__main__':
    main()
