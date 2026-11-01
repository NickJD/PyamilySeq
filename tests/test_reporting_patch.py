import ast, csv, math, tempfile
from collections import Counter, OrderedDict
from pathlib import Path
from types import SimpleNamespace

source = Path('/mnt/data/workflow_reporting_fixed.py').read_text()
tree = ast.parse(source)
func_names = {'_prevalence_bin', '_prevalence_bin_label', '_write_statistics', '_write_family_summary'}
namespace = {'math': math, 'Counter': Counter, 'csv': csv}
for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name in func_names:
        exec(compile(ast.Module(body=[node], type_ignores=[]), '<extracted>', 'exec'), namespace)

class Family:
    def __init__(self, family_id, stage, prevalence, group):
        self.family_id = family_id
        self.stage = stage
        self.prevalence = prevalence
        self.category = f'{stage}_core_{group}'
        self.copy_number = {'G0': 1}
        self.members = ['G0|gene1']
        self.representative = 'G0|gene1'
        self.source_families = ['source_1']
        self.source_clusters = ['1']

counts = {
    'First': {'99':319, '95':184, '15':1082, '0':45147},
    'extended': {'99':615, '95':112, '15':250, '0':816},
    'combined': {'99':0, '95':0, '15':0, '0':3},
    'only_Second': {'99':0, '95':0, '15':0, '0':1565},
}
N=5892
thresholds={'99':math.ceil(.99*N), '95':math.ceil(.95*N), '15':math.ceil(.15*N)}
families={}; categories=OrderedDict(); fid=0
for stage, bins in counts.items():
    for group, n in bins.items():
        categories[f'{stage}_core_{group}']=[]
        prevalence = thresholds[group] if group in thresholds else thresholds['15']-1
        for _ in range(n):
            key=f'F{fid}'; fid+=1
            fam=Family(key,stage,prevalence,group)
            families[key]=fam
            categories[f'{stage}_core_{group}'].append(key)
ledger=SimpleNamespace(families=families, taxa=[f'G{i}' for i in range(N)], genomes=[f'G{i}' for i in range(N)], assignment={}, sorted_families=lambda ids: [families[x] for x in ids])
with tempfile.TemporaryDirectory() as tmp:
    path=Path(tmp)/'summary_statistics.txt'
    namespace['_write_statistics'](str(path),ledger,list(families),categories,'Species',12054451,[99,95,15,0])
    out=path.read_text()
    assert 'ALL_STAGES\t934\t296\t1332\t47531\t50093' in out, out[-1200:]
    assert 'Families present in >= 99% of genomes (at least 5834/5892): 934' in out
    assert 'Families present in >= 95% of genomes (at least 5598/5892): 1230' in out
    assert 'Families present in >= 15% of genomes (at least 884/5892): 2562' in out
    family_path=Path(tmp)/'family_summary.tsv'
    namespace['_write_family_summary'](str(family_path), ledger, ['F0'], 'Species', [99,95,15,0])
    with family_path.open() as handle:
        row=next(csv.DictReader(handle,delimiter='\t'))
    assert row['prevalence_bin']=='>=99%'
    assert row['prevalence_denominator']=='5892'
    assert row['prevalence_percent']=='99.0156'
    print(out)