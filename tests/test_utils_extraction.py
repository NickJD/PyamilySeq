import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = Path(__file__).with_name('utils.py')
spec = importlib.util.spec_from_file_location('utils_under_test', MODULE_PATH)
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)


class ExtractionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write_combined(self, directory, genome, sequence='ATGAAATTT', gene_id='cds1'):
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f'{genome}.gff3'
        path.write_text(
            '##gff-version 3\n'
            f'chr1\tPyrodigal\tCDS\t1\t{len(sequence)}\t.\t+\t0\tID={gene_id}\n'
            '##FASTA\n'
            '>chr1\n'
            f'{sequence}\n',
            encoding='utf-8',
        )
        return path

    def test_combined_gff_uses_sorted_genome_order_and_translation(self):
        inp = self.root / 'combined'
        # Deliberately create files in reverse lexical order.
        self.write_combined(inp, 'Z_genome')
        self.write_combined(inp, 'A_genome')
        out = self.root / 'out' / 'genes_dna.fasta'
        out.parent.mkdir()

        utils.read_combined_files(str(inp), '.gff3', 'CDS', str(out), True, False)

        headers = [line for line in out.read_text().splitlines() if line.startswith('>')]
        self.assertEqual(headers, ['>A_genome|cds1', '>Z_genome|cds1'])
        aa = out.with_name('genes_aa.fasta').read_text()
        self.assertIn('>A_genome|cds1\nMKF', aa)
        self.assertIn('>Z_genome|cds1\nMKF', aa)

    def test_combined_gff_accepts_reverse_strand(self):
        inp = self.root / 'combined'
        inp.mkdir()
        target = 'ATGAAATTT'
        genomic = utils.reverse_complement(target)
        (inp / 'G1.gff3').write_text(
            '##gff-version 3\n'
            'chr1\tPyrodigal\tCDS\t1\t9\t.\t-\t0\tID=minus1\n'
            '##FASTA\n>chr1\n' + genomic + '\n',
            encoding='utf-8',
        )
        out = self.root / 'genes_dna.fasta'
        utils.read_combined_files(str(inp), '.gff3', 'CDS', str(out), False, False)
        self.assertIn('ATGAAATTT', out.read_text())

    def test_invalid_dna_preserves_preexisting_outputs_and_cleans_temps(self):
        inp = self.root / 'combined'
        self.write_combined(inp, 'G1', sequence='ATGAAAXTT')
        out = self.root / 'genes_dna.fasta'
        aa = self.root / 'genes_aa.fasta'
        out.write_text('OLD DNA\n')
        aa.write_text('OLD AA\n')

        with self.assertRaisesRegex(ValueError, 'STRICT NUCLEOTIDE VALIDATION FAILED'):
            utils.read_combined_files(str(inp), '.gff3', 'CDS', str(out), True, False)

        self.assertEqual(out.read_text(), 'OLD DNA\n')
        self.assertEqual(aa.read_text(), 'OLD AA\n')
        self.assertEqual(list(self.root.glob('.*.tmp')), [])

    def test_paired_gff_fasta_sorted_and_validated(self):
        inp = self.root / 'paired'
        inp.mkdir()
        for genome in ('Z_genome', 'A_genome'):
            (inp / f'{genome}.gff3').write_text(
                '##gff-version 3\n'
                'chr1\tPyrodigal\tCDS\t1\t9\t.\t+\t0\tID=cds1\n',
                encoding='utf-8',
            )
            (inp / f'{genome}.fna').write_text('>chr1\nATGAAATTT\n', encoding='utf-8')
        out = self.root / 'genes_dna.fasta'
        utils.read_separate_files(str(inp), '.gff3', '.fna', 'CDS', str(out), True, False)
        headers = [line for line in out.read_text().splitlines() if line.startswith('>')]
        self.assertEqual(headers, ['>A_genome|cds1', '>Z_genome|cds1'])
        self.assertIn('>A_genome|cds1\nMKF', (self.root / 'genes_aa.fasta').read_text())

    def test_paired_invalid_dna_preserves_existing_outputs(self):
        inp = self.root / 'paired'
        inp.mkdir()
        (inp / 'G1.gff3').write_text(
            '##gff-version 3\nchr1\tX\tCDS\t1\t9\t.\t+\t0\tID=g1\n',
            encoding='utf-8',
        )
        (inp / 'G1.fna').write_text('>chr1\nATGAAAXTT\n', encoding='utf-8')
        out = self.root / 'genes_dna.fasta'
        aa = self.root / 'genes_aa.fasta'
        out.write_text('OLD DNA\n')
        aa.write_text('OLD AA\n')

        with self.assertRaisesRegex(ValueError, 'STRICT NUCLEOTIDE VALIDATION FAILED'):
            utils.read_separate_files(str(inp), '.gff3', '.fna', 'CDS', str(out), True, False)

        self.assertEqual(out.read_text(), 'OLD DNA\n')
        self.assertEqual(aa.read_text(), 'OLD AA\n')

    def test_paired_missing_contig_fails_instead_of_silently_dropping_gene(self):
        inp = self.root / 'paired'
        inp.mkdir()
        (inp / 'G1.gff3').write_text(
            '##gff-version 3\nmissing\tX\tCDS\t1\t3\t.\t+\t0\tID=g1\n',
            encoding='utf-8',
        )
        (inp / 'G1.fna').write_text('>chr1\nATG\n', encoding='utf-8')
        out = self.root / 'genes_dna.fasta'
        with self.assertRaisesRegex(ValueError, 'absent from the associated FASTA'):
            utils.read_separate_files(str(inp), '.gff3', '.fna', 'CDS', str(out), False, False)
        self.assertFalse(out.exists())

    def test_fasta_file_input_is_sorted_and_validated_when_translated(self):
        inp = self.root / 'fasta'
        inp.mkdir()
        (inp / 'Z_genome.fna').write_text('>g1\nATGAAATTT\n', encoding='utf-8')
        (inp / 'A_genome.fna').write_text('>g1\nATGAAATTT\n', encoding='utf-8')
        out = self.root / 'genes_dna.fasta'
        utils.read_fasta_files(str(inp), '.fna', str(out), True, False)
        headers = [line for line in out.read_text().splitlines() if line.startswith('>')]
        self.assertEqual(headers, ['>A_genome|g1', '>Z_genome|g1'])
        self.assertIn('>A_genome|g1\nMKF', (self.root / 'genes_aa.fasta').read_text())

    def test_translate_dna_to_aa_is_strict_and_atomic(self):
        dna = self.root / 'input.fasta'
        aa = self.root / 'output.fasta'
        dna.write_text('>g1|cds1\nATGAAAXTT\n', encoding='utf-8')
        aa.write_text('OLD OUTPUT\n')
        with self.assertRaisesRegex(ValueError, 'STRICT NUCLEOTIDE VALIDATION FAILED'):
            utils.translate_dna_to_aa(str(dna), str(aa))
        self.assertEqual(aa.read_text(), 'OLD OUTPUT\n')

    def test_sorted_input_file_discovery(self):
        inp = self.root / 'in'
        inp.mkdir()
        for name in ('C.gff3', 'A.gff3', 'B.gff3'):
            (inp / name).write_text('', encoding='utf-8')
        self.assertEqual(
            [Path(p).name for p in utils._sorted_input_files(str(inp), '.gff3')],
            ['A.gff3', 'B.gff3', 'C.gff3'],
        )

    def test_gzip_fallback_is_sorted_and_readable(self):
        inp = self.root / 'compressed'
        inp.mkdir()
        for genome in ('Z', 'A'):
            with gzip.open(inp / f'{genome}.gff3.gz', 'wt', encoding='utf-8') as out:
                out.write(
                    '##gff-version 3\n'
                    'chr1\tX\tCDS\t1\t9\t.\t+\t0\tID=g1\n'
                    '##FASTA\n>chr1\nATGAAATTT\n'
                )
        dna = self.root / 'genes_dna.fasta'
        utils.read_combined_files(str(inp), '.gff3', 'CDS', str(dna), False, False)
        headers = [line for line in dna.read_text().splitlines() if line.startswith('>')]
        self.assertEqual(headers, ['>A|g1', '>Z|g1'])


    def test_missing_one_of_multiple_pairs_aborts_without_committing(self):
        inp = self.root / 'paired'
        inp.mkdir()
        (inp / 'A_present.gff3').write_text(
            '##gff-version 3\nchr1\tX\tCDS\t1\t3\t.\t+\t0\tID=g1\n',
            encoding='utf-8',
        )
        (inp / 'A_present.fna').write_text('>chr1\nATG\n', encoding='utf-8')
        (inp / 'B_missing.gff3').write_text(
            '##gff-version 3\nchr1\tX\tCDS\t1\t3\t.\t+\t0\tID=g1\n',
            encoding='utf-8',
        )
        out = self.root / 'genes_dna.fasta'
        out.write_text('OLD DNA\n')
        with self.assertRaisesRegex(FileNotFoundError, 'Missing paired FASTA files'):
            utils.read_separate_files(str(inp), '.gff3', '.fna', 'CDS', str(out), False, False)
        self.assertEqual(out.read_text(), 'OLD DNA\n')

    def test_invalid_utf8_in_gff_fails_without_replacing_outputs(self):
        inp = self.root / 'combined'
        inp.mkdir()
        path = inp / 'G1.gff3'
        path.write_bytes(
            b'##gff-version 3\n'
            b'chr1\tX\tCDS\t1\t9\t.\t+\t0\tID=g1\n'
            b'##FASTA\n>chr1\nATGAAA\xc5TT\n'
        )
        out = self.root / 'genes_dna.fasta'
        aa = self.root / 'genes_aa.fasta'
        out.write_text('OLD DNA\n')
        aa.write_text('OLD AA\n')
        with self.assertRaises(UnicodeDecodeError):
            utils.read_combined_files(str(inp), '.gff3', 'CDS', str(out), True, False)
        self.assertEqual(out.read_text(), 'OLD DNA\n')
        self.assertEqual(aa.read_text(), 'OLD AA\n')

    def test_coordinate_out_of_bounds_is_reported(self):
        inp = self.root / 'combined'
        inp.mkdir()
        (inp / 'G1.gff3').write_text(
            '##gff-version 3\n'
            'chr1\tX\tCDS\t1\t10\t.\t+\t0\tID=g1\n'
            '##FASTA\n>chr1\nATGAAATTT\n',
            encoding='utf-8',
        )
        out = self.root / 'genes_dna.fasta'
        with self.assertRaisesRegex(ValueError, 'only 9 bp long'):
            utils.read_combined_files(str(inp), '.gff3', 'CDS', str(out), False, False)
        self.assertFalse(out.exists())

    def test_read_genomes_returns_deterministic_sorted_order(self):
        fasta = self.root / 'genomes.fasta'
        fasta.write_text(
            '>gene3|GenomeC\nATG\n'
            '>gene1|GenomeA\nATG\n'
            '>gene2|GenomeB\nATG\n',
            encoding='utf-8',
        )
        self.assertEqual(
            utils.read_genomes_from_fasta(str(fasta)),
            ['GenomeA', 'GenomeB', 'GenomeC'],
        )


if __name__ == '__main__':
    unittest.main(verbosity=2)
