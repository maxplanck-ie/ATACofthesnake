import os

import pandas as pd

from aos import helper


class TestIdxToMit:
    def test_fraction(self, tmp_path):
        f = tmp_path / "idx.tsv"
        f.write_text("chr1 100\nMT 50\nchr2 20\n")
        assert helper.idx_to_mit(f) == 0.29


class TestPCAColors:
    def test_no_samplesheet(self):
        assert helper.PCA_colors(None, ["s1", "s2"]) == ""
        assert helper.PCA_colors("", ["s1", "s2"]) == ""

    def test_colors_grouped_by_first_column(self, tmp_path):
        ss = tmp_path / "ss.tsv"
        ss.write_text(
            "sample\tgroup\n"
            "s1\tA\n"
            "s2\tA\n"
            "s3\tB\n"
            "s4\tB\n"
        )
        result = helper.PCA_colors(ss, ["s1", "s2", "s3", "s4"])
        assert result == '--colors "#1f77b4" "#1f77b4" "#ff7f0e" "#ff7f0e"'


class TestGetElbow:
    def test_obvious_elbow(self):
        k_range = [1, 2, 3, 4, 5]
        inertias = [100, 10, 9, 8, 7]
        assert helper.get_elbow(inertias, k_range) == 2


class TestMergeIdx:
    def test_merge(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        os.makedirs("qc")
        with open("qc/s1_ix.tsv", "w") as f:
            f.write("chr1\t10\nchr2\t20\n*\t0\n")
        with open("qc/s2_ix.tsv", "w") as f:
            f.write("chr1\t15\nchr2\t25\n*\t0\n")

        helper.merge_idx(["qc/s1_ix.tsv", "qc/s2_ix.tsv"], "qc/ixstat.tsv")

        lines = open("qc/ixstat.tsv").read().splitlines()
        assert lines == [
            "samples\ts1\ts2",
            "chr1\t10\t15",
            "chr2\t20\t25",
        ]


class TestMergeSieve:
    def test_merge(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        os.makedirs("qc")
        with open("qc/s1_sieve.txt", "w") as f:
            f.write("# comment\ninput/s1.bam\t80\t100\n")
        with open("qc/s2_sieve.txt", "w") as f:
            f.write("# comment\ninput/s2.bam\t45\t50\n")

        helper.merge_sieve(["qc/s1_sieve.txt", "qc/s2_sieve.txt"], "qc/sieve.tsv")

        lines = open("qc/sieve.tsv").read().splitlines()
        assert lines == [
            "sample\ts1\ts2",
            "surviving\t80\t45",
            "initial\t100\t50",
            "fraction\t0.8\t0.9",
        ]


class TestPeakBoundaries:
    def test_clips_peaks_to_chrom_length(self, tmp_path):
        genomefa = tmp_path / "genome.fa"
        genomefa.write_text(">chr1\nATCGATCGATCG\n")  # 12 bases

        peaks = tmp_path / "peaks.bed"
        peaks.write_text("chr1\t5\t20\nchr1\t0\t8\n")

        of = tmp_path / "out.bed"
        helper.peak_boundaries(peaks, genomefa, None, of)

        df = pd.read_csv(of, sep="\t", header=None)
        assert df.values.tolist() == [["chr1", 5, 12], ["chr1", 0, 8]]

    def test_copies_existing_peakset(self, tmp_path):
        peakset = tmp_path / "peakset.bed"
        peakset.write_text("chr1\t0\t10\n")
        of = tmp_path / "out.bed"

        helper.peak_boundaries(None, None, peakset, of)

        assert of.read_text() == peakset.read_text()


class TestPcaToMqc:
    def test_transposes_and_labels_variance(self, tmp_path):
        raw = tmp_path / "pca_data.tsv"
        raw.write_text(
            "Component\tsampleA.scalefac.bw\tsampleB.scalefac.bw\tEigenvalue\n"
            "1\t10.0\t-10.0\t300.0\n"
            "2\t2.0\t-2.0\t50.0\n"
        )
        of = tmp_path / "PCA_mqc.tsv"

        helper.pca_to_mqc(raw, of)

        header_lines = [
            line for line in of.read_text().splitlines() if line.startswith("#")
        ]
        assert any("85.7% variance" in line for line in header_lines)
        assert any("14.3% variance" in line for line in header_lines)

        df = pd.read_csv(of, sep="\t", comment="#", index_col=0)
        assert list(df.columns) == ["PC1", "PC2"]
        assert list(df.index) == ["sampleA", "sampleB"]
        assert df.loc["sampleA", "PC1"] == 10.0
        assert df.loc["sampleA", "PC2"] == 2.0
        assert df.loc["sampleB", "PC1"] == -10.0
        assert df.loc["sampleB", "PC2"] == -2.0


class TestPlotFrip:
    def test_produces_png(self, tmp_path):
        fs = tmp_path / "fripscores_mqc.tsv"
        fs.write_text("# id: \"frip_scores\"\nSample\tFRiP\ns1\t0.5\ns2\t0.7\n")
        of = tmp_path / "fripscores.png"

        helper.plotfrip(fs, of)

        assert of.exists()
        assert of.stat().st_size > 0


class TestPlotFragsize:
    def test_produces_png(self, tmp_path):
        fs = tmp_path / "fragsize.tsv"
        fs.write_text(
            "Sample\tSize\tOccurrences\n"
            "input/s1.bam\t100\t3\n"
            "input/s1.bam\t150\t2\n"
            "input/s1.bam\t200\t1\n"
            "input/s2.bam\t110\t2\n"
            "input/s2.bam\t160\t3\n"
            "input/s2.bam\t210\t1\n"
        )
        of = tmp_path / "fragmentsizes.png"

        helper.plotfragsize(fs, of)

        assert of.exists()
        assert of.stat().st_size > 0


class TestPlotIxs:
    def test_produces_png(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        os.makedirs("figures")
        ixs = tmp_path / "ixstat.tsv"
        ixs.write_text("samples\ts1\ts2\nchr1\t80\t90\nchrM\t20\t10\n")

        helper.plotixs(ixs, "chrM")

        assert (tmp_path / "figures" / "mitofraction.png").exists()


class TestPlotSieve:
    def test_produces_png(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        os.makedirs("figures")
        sieve = tmp_path / "sieve.tsv"
        sieve.write_text(
            "sample\ts1\ts2\n"
            "surviving\t80\t45\n"
            "initial\t100\t50\n"
            "fraction\t0.8\t0.9\n"
        )

        helper.plotsieve(sieve)

        assert (tmp_path / "figures" / "alignmentsieve.png").exists()
