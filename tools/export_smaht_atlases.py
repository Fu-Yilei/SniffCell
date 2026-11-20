"""Package the frozen 2026-10-16 SMaHT catalogs without read-level data.

Requires the authorized local SMaHT atlas/build artifacts. Existing destinations
are never overwritten. Tissue TSV contents are preserved byte for byte; blood
changes are restricted to three non-scoring provenance columns.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path


TISSUES = {
    '3E': ('colon_ascending', ['3E']),
    '3G': ('colon_descending', ['3G']),
    '3I': ('liver', ['3I']),
    '3O': ('aorta', ['3O']),
    '3Q': ('lung', ['3Q']),
    '3S': ('heart', ['3S']),
    '3AD': ('skin', ['3AD', '3AF']),
}
METADATA_CHANGES = {'atlas_source', 'paired_stability_filter', 'paired_donors'}


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            digest.update(block)
    return digest.hexdigest()


def stats(path):
    counts = collections.Counter()
    orders = set()
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            counts[row.get('modification') or 'modifiedC'] += 1
            orders.add(row['code_order'])
    if len(orders) != 1:
        raise ValueError(f'Expected a single shared class schema: {path.name}')
    return dict(rows=sum(counts.values()), rows_by_modification=dict(counts),
                classes=next(iter(orders)).split('|'))


def compress_exact(source, destination, expected_sha):
    digest = hashlib.sha256()
    with source.open('rb') as inp, destination.open('xb') as raw:
        with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as output:
            for block in iter(lambda: inp.read(2**20), b''):
                digest.update(block)
                output.write(block)
    if digest.hexdigest() != expected_sha:
        raise ValueError(f'Source differs from frozen run manifest: {source.name}')
    return dict(compressed_sha256=sha(destination), uncompressed_sha256=digest.hexdigest(),
                source_filename=source.name, source_sha256=expected_sha,
                source_content_preserved=True, **stats(destination))


def export_blood(source, destination, expected_sha):
    if sha(source) != expected_sha:
        raise ValueError('Blood source differs from frozen production catalog')
    original_scoring = hashlib.sha256()
    exported_scoring = hashlib.sha256()
    with gzip.open(source, 'rt', newline='') as inp, destination.open('xb') as raw:
        reader = csv.DictReader(inp, delimiter='\t')
        fields = reader.fieldnames
        scoring_fields = [field for field in fields if field not in METADATA_CHANGES]
        with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as zipped:
            with io.TextIOWrapper(zipped, newline='') as output:
                writer = csv.DictWriter(output, fieldnames=fields, delimiter='\t', lineterminator='\n')
                writer.writeheader()
                for row in reader:
                    original_scoring.update(json.dumps([row[k] for k in scoring_fields], separators=(',', ':')).encode()+b'\n')
                    row['atlas_source'] = Path(row['atlas_source']).name
                    # The original find export carried stale false/zero metadata
                    # despite the subsequent two-donor stability filter.
                    if row['modification'] in {'5mC', '5hmC'}:
                        row['paired_stability_filter'] = 'True'
                        row['paired_donors'] = '2'
                    if any(v.startswith('/') for v in row.values()):
                        raise ValueError('Unexpected local path in blood export')
                    exported_scoring.update(json.dumps([row[k] for k in scoring_fields], separators=(',', ':')).encode()+b'\n')
                    writer.writerow(row)
    if original_scoring.digest() != exported_scoring.digest():
        raise ValueError('Export altered a scoring field')
    return dict(source_filename=source.name, source_compressed_sha256=expected_sha,
                compressed_sha256=sha(destination), scoring_fields=scoring_fields,
                scoring_fields_sha256=original_scoring.hexdigest(),
                scoring_fields_preserved=True,
                metadata_changes={'atlas_source': 'basename only',
                                  'paired_stability_filter': 'True on FANS m/h rows',
                                  'paired_donors': '2 on FANS m/h rows'}, **stats(destination))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smaht-root', type=Path, required=True)
    parser.add_argument('--atlas-root', type=Path, default=Path('atlases'))
    args = parser.parse_args()
    root = args.smaht_root/'sniffcell'
    release = args.atlas_root/'smaht_20261016'
    blood_dir = args.atlas_root/'blood_5hmc'
    release.mkdir(parents=True, exist_ok=True)
    blood_dir.mkdir(parents=True, exist_ok=True)
    audit = root/'plans/all_smaht_20261007_v2/nonblood_submission_20261007/atlas_audit.tsv'
    records = {Path(r['atlas']).parent.name: r for r in csv.DictReader(audit.open(), delimiter='\t')}
    entries = []
    for code, (name, codes) in TISSUES.items():
        record = records[code]
        source = root/f'ctDMR/sniffcell_lite_find/{code}/{code}.ctdmr.tsv'
        dest = release/f'{name}.ctdmr.tsv.gz'
        info = compress_exact(source, dest, record['sha256'])
        if info['rows'] != int(record['rows']):
            raise ValueError(f'Unexpected row count for {code}')
        if code == '3AD' and records['3AF']['sha256'] != record['sha256']:
            raise ValueError('Skin catalogs no longer identical')
        entries.append(dict(tissue=name, tissue_codes=codes, file=dest.name, **info))
        print(f'Packaged {name}: {info["rows"]} rows', flush=True)

    build = root/'atlases/whole_blood_20261016/release'
    qc = json.loads((build/'integrated_atlas_QC.json').read_text())
    source = build/'blood_universal.modifiedC_plus_5mc_5hmc.ctdmr.tsv.gz'
    dest = blood_dir/source.name
    info = export_blood(source, dest, qc['atlas_sha256'])
    if info['rows_by_modification'] != qc['rows_by_modification']:
        raise ValueError('Blood assay counts changed')
    backbone = json.loads((build/'atlas_QC.json').read_text())
    gate = json.loads((root/'plans/whole_blood_integrated_20261016/pilot_gate.json').read_text())
    if gate['technical_gate'] != 'pass':
        raise ValueError('Production validation gate did not pass')
    fractions = [r['expected_class_compatible_fraction'] for r in gate['reference_diagnostics']]
    blood_manifest = dict(reference_build='GRCh38', file=dest.name, **info,
        source_backbone_sha256=qc['source_backbone_sha256'], source_fans_sha256=qc['source_fans_sha256'],
        discovery_parameters=backbone['parameters'],
        reference_profile_counts=backbone['reference_profile_counts'],
        minimum_finite_backbone_references_per_class=2,
        source_files={k: v['sha256'] for k, v in backbone['source_files'].items()},
        fans_reference=dict(donors=2, sorted_samples=8, classes=['B-cell','Monocyte','NK-cell','T-cell'],
            all_blocks_observed=True, same_effect_direction_both_donors=True,
            minimum_effect_each_donor=0.30, minimum_median_effect=0.40,
            candidate_rows={'5mC':73676,'5hmC':55}, stable_rows={'5mC':71846,'5hmC':48}),
        integration='Preserve modifiedC backbone; add stable m/h rows outside its extraction intervals',
        granulocyte_mh_profiles_available=False,
        validation=dict(platform_smoke_checks=len(gate['pilots']), reciprocal_donor_checks=len(fractions),
            expected_class_compatible_fraction_range=[min(fractions),max(fractions)],
            scope='Bounded known-origin read checks; compatible multiple labels count; no granulocyte FANS control; not genome-wide accuracy'))
    (blood_dir/'manifest.json').write_text(json.dumps(blood_manifest, indent=2)+'\n')
    entries.append(dict(tissue='whole_blood', tissue_codes=['3A'], file='../blood_5hmc/'+dest.name, **info))
    brain = args.atlas_root/'brain_5hmc/brain_universal.paper_backbone_plus_5mc_5hmc.four_donors.ctdmr.tsv.gz'
    if sha(brain) != records['brain_5hmc']['sha256']:
        raise ValueError('Branch brain catalog differs from production input')
    entries.append(dict(tissue='brain', tissue_codes=['3AK','3AL','3AM','3AN','3AO'],
        file='../brain_5hmc/'+brain.name, compressed_sha256=sha(brain), source_content_preserved=True, **stats(brain)))
    manifest = dict(release='SMaHT 2026-10-16 production catalogs', reference_build='GRCh38',
        original_code_commit='26e0ad9f6b5ebfd853910feacc9fa50037872861', catalogs=entries,
        tissue_discovery_parameters=dict(diff_threshold=0.4,min_rows=1,min_cpgs=3,max_gap_bp=2000),
        required_scoring_fixes=['observed-CpG means for all assays','explicit missing granulocyte m/h profiles'],
        additional_reporting_correction='Per-read variant readable labels now use consensus codes; historical readable outputs were not regenerated',
        scope=dict(original_inventory=560,completed=407,blood=34,nonblood=373,
                   existing_callset_annotations=405,without_callset=2),
        not_run=dict(missing_tissue_catalogs={'esophagus':44,'muscle':35,'testis':37,'adrenal':22,'ovary':8},
                     fibroblast_controls=7))
    (release/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print('Verified branch brain catalog and exported blood scoring fields', flush=True)


if __name__ == '__main__':
    main()
