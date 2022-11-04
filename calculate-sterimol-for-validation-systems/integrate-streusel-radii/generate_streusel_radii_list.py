import argparse
import pandas as pd


def _get_command_line_arguments() -> argparse.Namespace():
    parser = argparse.ArgumentParser()
    parser.add_argument("radii_csv")
    return parser.parse_args()


def main() -> None:
    cli_args = _get_command_line_arguments()

    radii_df = pd.read_csv(cli_args.radii_csv)

    # extract neutral atoms
    neutral_atoms = radii_df[radii_df['charge state'] == 0]
    # remove duplicates
    neutral_atoms.duplicated(subset=['atom'],keep='last')
    # save to csv
    neutral_atoms.to_csv('streusel_neutral_atom_radii.csv', index=False)

    with open('streusel_neutral_atom_radii_dict.txt', 'w') as radii_dictionary:
        radii_dict = []

        radii_dict.append('streusel_radi: typing.Final = {')

        for index, row in neutral_atoms.iterrows():
            radii_dict.append(f"    '{row['atom']}': {row['radius']},")
        radii_dict.append('}')
        radii_dict_content = '\n'.join(radii_dict)
        radii_dictionary.write(f'{radii_dict_content}\n')



if __name__ == '__main__':
    main()
