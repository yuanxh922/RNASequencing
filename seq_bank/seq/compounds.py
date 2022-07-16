import attr
import pandas as pd

@attr.s
class Compound:
    mass = attr.ib(default=0)
    rt = attr.ib(default=0)
    vol = attr.ib(default=0)
    width = attr.ib(default=1.0)
    qs = attr.ib(default=100)

    @classmethod
    def cpds2df(cls, cpds):
        col_mass = [cpd.mass for cpd in cpds]
        col_rt = [cpd.rt for cpd in cpds]
        col_vol = [cpd.vol for cpd in cpds]
        col_width = [cpd.width for cpd in cpds]
        col_qs = [cpd.qs for cpd in cpds]
        df = pd.DataFrame({
            'Mass': col_mass,
            'RT': col_rt,
            'Vol': col_vol,
            'Width': col_width,
            'Quality Score': col_qs
        })
        df['Cpd'] = range(1, len(df) + 1)
        return df

class SubCompound(Compound):
    pass