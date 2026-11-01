import gzip
import tempfile
import unittest
from pathlib import Path

try:
    from PyamilySeq.incremental import prepare_incremental_inputs, IncrementalInputError
except ImportError:
    from incremental import prepare_incremental_inputs, IncrementalInputError


class TestPrepareIncrementalInputs(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.out = self.root / 'out'
        self.out.mkdir()
        self.original_fasta = self.root / 'old.fasta'
        self.original_fasta.write_text(
            '>g1|gene1\nATGAAA\n'
            '>g1|gene1_alt\nATGCCC\n'
            '>g2|gene2\nATGGGG\n', encoding='utf-8'
        )
        self.cluster_file = self.root / 'old.clstr'
        self.cluster_file.write_text(
            '>Cluster 0\n'
            '0 6nt, >g1|gene1... *\n'
            '1 6nt, >g1|gene1_alt... at 80.00%\n'
            '>Cluster 1\n'
            '0 6nt, >g2|gene2... *\n', encoding='utf-8'
        )
        self.new_fasta = self.root / 'new.fasta'
        self.new_fasta.write_text(
            '>g3|gene3 note ignored\nATGTTT\n'
            '>g3|gene4\nATGCCC\n', encoding='utf-8'
        )

    def tearDown(self):
        self.temp.cleanup()

    def test_builds_full_and_representative_fastas(self):
        result = prepare_incremental_inputs(
            str(self.cluster_file), str(self.original_fasta), str(self.new_fasta), str(self.out)
        )
        full = (self.out / 'incremental_combined_sequences_dna.fasta').read_text()
        reps = (self.out / 'incremental_representatives_dna.fasta').read_text()
        self.assertEqual(result['original_sequence_count'], 3)
        self.assertEqual(result['original_cluster_count'], 2)
        self.assertEqual(result['representative_count'], 2)
        self.assertEqual(result['new_gene_count'], 2)
        self.assertEqual(result['new_genomes'], ['g3'])
        self.assertIn('>g1|gene1_alt\nATGCCC\n', full)
        self.assertIn('>g3|gene4\nATGCCC\n', full)
        self.assertIn('>g1|gene1\nATGAAA\n', reps)
        self.assertIn('>g2|gene2\nATGGGG\n', reps)
        self.assertNotIn('>g1|gene1_alt\n', reps)
        self.assertIn('>g3|gene3\nATGTTT\n', reps)

    def test_duplicate_new_id_raises(self):
        self.new_fasta.write_text('>g3|gene3\nATGAAA\n>g3|gene3\nATGCCC\n', encoding='utf-8')
        with self.assertRaisesRegex(IncrementalInputError, 'Duplicate sequence identifier'):
            prepare_incremental_inputs(
                str(self.cluster_file), str(self.original_fasta), str(self.new_fasta), str(self.out)
            )

    def test_old_new_id_collision_raises(self):
        self.new_fasta.write_text('>g1|gene1\nATGAAA\n', encoding='utf-8')
        with self.assertRaisesRegex(IncrementalInputError, 'reuses'):
            prepare_incremental_inputs(
                str(self.cluster_file), str(self.original_fasta), str(self.new_fasta), str(self.out)
            )

    def test_missing_original_representative_raises(self):
        self.original_fasta.write_text('>g1|gene1_alt\nATGCCC\n', encoding='utf-8')
        with self.assertRaisesRegex(IncrementalInputError, 'absent from'):
            prepare_incremental_inputs(
                str(self.cluster_file), str(self.original_fasta), str(self.new_fasta), str(self.out)
            )

    def test_generated_delimiter_can_be_rewritten(self):
        self.new_fasta.write_text('>g3|gene3\nATGTTT\n', encoding='utf-8')
        result = prepare_incremental_inputs(
            str(self.cluster_file), str(self.original_fasta), str(self.new_fasta), str(self.out),
            sequence_delimiter='::', generated_ids_use_pipe=True,
        )
        normalised = Path(result['new_genes_dna_fasta']).read_text(encoding='utf-8')
        self.assertIn('>g3::gene3\nATGTTT\n', normalised)
        self.assertEqual(result['new_genomes'], ['g3'])

    def test_gzip_inputs(self):
        cluster_gz = self.root / 'old.clstr.gz'
        with gzip.open(cluster_gz, 'wt', encoding='utf-8') as handle:
            handle.write(self.cluster_file.read_text())
        fasta_gz = self.root / 'old.fasta.gz'
        with gzip.open(fasta_gz, 'wt', encoding='utf-8') as handle:
            handle.write(self.original_fasta.read_text())
        result = prepare_incremental_inputs(
            str(cluster_gz), str(fasta_gz), str(self.new_fasta), str(self.out)
        )
        self.assertEqual(result['representative_count'], 2)


if __name__ == '__main__':
    unittest.main()