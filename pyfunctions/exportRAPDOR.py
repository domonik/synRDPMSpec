from RAPDOR.datastructures import RAPDORData
import copy
import numpy as np

def main(file, outfile):
    rapdor_data = RAPDORData.from_file(file)
    def sort_fct(name):
        remain, frac = name.split("-")
        frac = int(frac)
        rep = int(remain[-1])
        return frac, rep, remain
    data_cols = rapdor_data._data_cols.tolist()
    data_cols.sort(key=sort_fct)

    df = rapdor_data.df[["RAPDORid"] + data_cols]
    df[data_cols] = 2 ** df[data_cols]
    df = df.fillna(0)

    df = df.rename({"RAPDORid": "name"}, axis=1)
    df.to_csv(outfile, sep=";", index=False)



if __name__ == '__main__':
    file = snakemake.input.json
    outfile = snakemake.output.csv
    main(file, outfile)

