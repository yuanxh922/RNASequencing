import pandas as pd
import numpy as np

from collections import namedtuple
ShiftHit = namedtuple('ShiftHit', 'shift hit')

def match_dfs(df_src, df_dst, ppm=10, shift=0):
    """
    find the subset contains common Mass values that exist in both df_src and df_dst
    """
    def _find_mass(df, mass, ppm=10):
        if df.empty:
            return df
        df = df[(df.Mass < mass+1) & (df.Mass > mass-1)]
        if df.shape[0] == 0:
            return df
        
        df_ppm = abs(1E6 * (df.Mass - mass)) / mass
        mask = df_ppm < ppm
        df_found = df[mask].copy()
        df_found['PPM'] = df_ppm[mask]
        return df_found
    
    df_src = df_src.copy()
    if shift != 0:
        df_src.Mass += shift
    idxs = list()
    for idx, row in df_src.iterrows():
        mass = row.Mass
        df_res = _find_mass(df_dst, mass, ppm)
        if not df_res.empty:
            df_src.loc[idx, 'Match'] = True
            idxs.extend(list(df_res.index))
    
    idxs = list(set(idxs))
    df_common = df_dst[df_dst.index.isin(idxs)]
    return df_common.copy()

def test_func(df_sample, shift):
    print('shift', shift)
    return '{}'.format(shift)

def func(df_sample, shift):
    dfm = match_dfs(df_sample, df_sample, shift=shift)
    sh = ShiftHit(shift=shift, hit=dfm.shape[0])
    return sh

def thermo_df(df, key_rows_only=True):
    df = df.rename(columns={'Monoisotopic Mass': 'Mass', 'Apex RT': 'RT', 'Sum Intensity': 'Vol',
                           'Relative Abundance': 'RA', 'Fractional Abundance': 'FA'})
    if key_rows_only:
        try:
            vols = ['Mass', 'RT', 'Vol', 'RA', 'FA']
            df = df[vols].dropna()
        except KeyError as err:
            vols = ['Mass', 'RT', 'Vol']
            df = df[vols]
        df = df.astype('float64')
    return df

def load_data(fpath, csv_format=False):
    func = pd.read_csv if csv_format else pd.read_excel
    df = func(fpath)
    df = thermo_df(df)
    return df

if __name__ == '__main__':
    import multiprocessing

    df_cmc = load_data('/Users/xyuan/Documents/SeqDataSets/Modifications/181227s07_100.xls')
    df_sample = df_cmc[(df_cmc.Mass<8000)&(df_cmc.RT<15)]

    PROCESSES = 4
    params = [(df_sample, shift) for shift in np.arange(251.0, 251.3, 0.1)]
    with multiprocessing.Pool(PROCESSES) as pool:
        results = pool.starmap(func, params)

    results = [result for result in results]
    df_shift_hits = pd.DataFrame(results)
    print(df_shift_hits)
    df_shift_hits.to_excel('./ShiftHits.xlsx')
    print("Done. Please find the result file 'ShiftHits.xlsx' under the current directory.")
