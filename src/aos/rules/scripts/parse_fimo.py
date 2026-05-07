import pandas as pd
from pathlib import Path

# input tsvs
tsvs = [Path(tsv) for tsv in snakemake.input.tsvs]
comp_group_dic = {}
for tsv in tsvs:
    comp = tsv.parts[2]
    if comp not in comp_group_dic:
        comp_group_dic[comp] = set()
    comp_group_dic[comp].add(tsv)

def parse_fimo(tsv):
    motnames = {}
    try:
        df = pd.read_table(tsv, comment='#')
    except pd.errors.EmptyDataError:
        return motnames
    motifs = list(df['motif_id'].unique())
    comp = tsv.parts[2]
    group = tsv.parts[3]
    for mot in motifs:
        _tdf = df[df['motif_id'] == mot]
        of = Path(f'footprints/plotaggregate/{comp}/bedfiles/{mot}-{group}.bed')
        with open(of, 'w') as o:
            for _, row in _tdf.iterrows():
                motid = row['motif_id']
                motidname = row['motif_alt_id']
                seqnam = row['sequence_name']
                chrom = seqnam.split(':')[0]
                start = int(seqnam.split(':')[1].split('-')[0])
                motstart = int(row['start'])
                motend = int(row['stop'])
                o.write(f"{chrom}\t{start+motstart-1}\t{start+motend}\n")
                motnames[motid] = motidname
    return motnames


for comp in comp_group_dic:
    motif_names = {}
    of = Path(f'footprints/plotaggregate/{comp}/bedfiles')
    of.mkdir(exist_ok=True, parents=True)
    # delete all files present already (checkpoint)
    for f in of.iterdir():
        if f.is_file():
            f.unlink()
    for tsv in comp_group_dic[comp]:
        motif_names.update(parse_fimo(tsv))
    # Save motif names for later use in plotting
    # remove file if it exists
    _motnamefile = Path(f'footprints/plotaggregate/{comp}/{comp}_motif_names.tsv')
    if _motnamefile.exists():
        _motnamefile.unlink()
    motnames_df = pd.DataFrame.from_dict(motif_names, orient='index', columns=['motif_alt_id'])
    motnames_df.index.name = 'motif_id'
    motnames_df.to_csv(_motnamefile, sep='\t')