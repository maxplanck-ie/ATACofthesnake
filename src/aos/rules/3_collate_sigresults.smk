checkpoint collate_sigresults:
  input:
    sigresults = SIGRESULTS
  output:
    directory("motifs")
  params:
    min_sigpeaks = config['cutoffs']['min_sigpeaks'],
    perm_cutoff = config['cutoffs']['permutation_cutoff'],
    fdr_cutoff = config['cutoffs']['fdr_cutoff'],
    lfc_cutoff = config['cutoffs']['lfc_cutoff'],
  run:
    import pandas as pd
    from pathlib import Path
    for sigresult in input.sigresults:
      match Path(sigresult).parents[0]:
        case "twogroup":
          # f"twogroup/{comp}/{comp}_heatmap.done" -> f"twogroup/{comp}/{comp}_diffacc_edgeR.tsv",
          _tsv = sigresult.replace("_heatmap.done", "_diffacc_edgeR.tsv")
          _comp = Path(sigresult).parents[1].name
          df = pd.read_table(_tsv, sep='\t')
          sig = df[df['FDR'] < params.fdr_cutoff]
          if len(sig) > params.min_sigpeaks:
            output_folder = Path(f"motifs/{_comp}")
            output_folder.mkdir(parents=True, exist_ok=True)
            # Upregs
            up = sig[sig['logFC'] > params.lfc_cutoff]
            down = sig[sig['logFC'] < -params.lfc_cutoff]
            for _df in [up, down]:
              if len(_df) > 0:
                assert _df['group_assignment'].unique() == 1
                groupname = _df['group_assignment'].unique()[0]
                with open(output_folder / f"{groupname}.bed", "w") as f:
                  for peak in _df['peak_id']:
                    f.write('\t'.join(peak.split('|'))+'\n')
        case "lrt":
          #f"lrt/{comp}/{comp}_heatmap.done" -> f"lrt/{comp}/{comp}_sig_peaks_kmeans.tsv"
          _tsv = sigresult.replace("_heatmap.done", "_sig_peaks_kmeans.tsv")
          if Path(_tsv).exists():
            df = pd.read_table(_tsv, sep='\t')
            _comp = Path(sigresult).parents[1].name
            output_folder = Path(f"motifs/{_comp}")
            output_folder.mkdir(parents=True, exist_ok=True)
            for k in df['k'].unique():
              with open(output_folder / f"{_comp}_k{k}.bed", "w") as f:
                for peak in df[df['k'] == k]['peak_id']:
                  f.write('\t'.join(peak.split('|'))+'\n')
        case "gp":
          #f"gp/{comp}/{comp}_postprocess.done" -> f"gp/{comp}/{comp}_k_table.tsv"
          #f"gp/{comp}/inttest_{comp}_{interaction}_postprocess.done" -> f"gp/{comp}/inttest_{comp}_{interaction}_k_table.tsv"
          _tsv = sigresult.replace("_postprocess.done", "_k_table.tsv")
          if Path(_tsv).exists():
            df = pd.read_table(_tsv, sep='\t', index_col=0)
            assert (df["FDR"] < params.perm_cutoff).all()
            _comp = Path(sigresult).parents[1].name
            if "inttest" in sigresult:
              _comp += "_" + Path(sigresult).name.split("_")[2]
            output_folder = Path(f"motifs/{_comp}")
            output_folder.mkdir(parents=True, exist_ok=True)
            for k in df['k'].unique():
              with open(output_folder / f"{_comp}_k{k}.bed", "w") as f:
                for peak in df[df['k'] == k].index:
                  f.write('\t'.join(peak.split('|'))+'\n')
