#!/usr/bin/env python3
"""Independent integrity checks on the actual exported assignment results."""
import csv
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def table(filename, delimiter='\t'):
    with (ROOT/'results'/filename).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def main():
    vcf = [line.split('\t') for line in (ROOT/'input/student_06.vcf').read_text().splitlines()
           if line and not line.startswith('#')]
    with zipfile.ZipFile(ROOT/'input/variants.zip') as archive:
        assert (ROOT/'input/student_06.vcf').read_bytes() == archive.read('student_06.vcf')
    annotations = table('annotations.tsv')
    assert len(annotations) == len(vcf) == 50
    assert annotations == table('annotations.csv', ',')
    assert len(table('clinvar.tsv')) == len(table('gnomad.tsv')) == 50
    assert [r['original_id'] for r in annotations] == [r[2] for r in vcf]
    assert len({r['variant_id'] for r in annotations}) == 50
    required = ['original_id', 'gene', 'most_severe_consequence', 'selected_transcript',
                'transcript_selection_rule', 'rsid', 'clinvar_present', 'gnomad_status']
    assert all(r[k] not in ['', 'NA'] for r in annotations for k in required)
    transitions = sum(a[3]+a[4] in ['AG','GA','CT','TC'] for a in vcf)
    assert transitions == 41 and len(vcf)-transitions == 9
    assert sum(r[0] == '21' for r in vcf) == 11
    assert sum(r[0] == '22' for r in vcf) == 39
    assert sum(r['most_severe_consequence'] == 'missense_variant' for r in annotations) == 32
    assert sum(r['most_severe_consequence'] == 'synonymous_variant' for r in annotations) == 18
    tx = table('vep_transcripts.tsv')
    for r in annotations:
        selected = [t for t in tx if t['original_id'] == r['original_id'] and t['is_selected'] == '1']
        assert len(selected) == 1
        assert r['most_severe_consequence'] in selected[0]['consequences'].split('|')
        raw = json.loads((ROOT/'raw/gnomad'/(r['variant_id']+'.json')).read_text())
        variant = raw.get('data', {}).get('variant')
        if r['gnomad_af'] != 'NA':
            assert variant and r['gnomad_status'] == 'present'
            data = variant[r['gnomad_af_source']]
            assert data['filters'] == [] and data['an'] > 0
            assert 0 <= float(r['gnomad_af']) <= 1
            assert math.isclose(float(r['gnomad_af']), data['ac']/data['an'], abs_tol=1e-12)
        if r['clinvar_present'] == 'no':
            assert r['clinvar_classification'] == r['clinvar_variation_id'] == 'NA'
    # Verify classifications against exact alleles in original NCBI XML, not rsID alone.
    archives = ET.parse(ROOT/'raw/clinvar.xml').getroot().findall('.//VariationArchive')
    for r in annotations:
        matches = []
        for v in archives:
            loc = v.find('./ClassifiedRecord/SimpleAllele/Location/SequenceLocation[@Assembly="GRCh38"]')
            if loc is not None and (loc.get('Chr'), loc.get('positionVCF'),
                    loc.get('referenceAlleleVCF'), loc.get('alternateAlleleVCF')) == (
                    r['chrom'], r['pos'], r['ref'], r['alt']):
                matches.append(v)
        assert (r['clinvar_present'] == 'yes') == bool(matches)
        for v in matches:
            assert v.get('VariationID') in r['clinvar_variation_id'].split(' || ')
            assert v.findtext('./ClassifiedRecord/Classifications/GermlineClassification/Description') in r['clinvar_classification']
    by_id = {r['original_id']: r for r in annotations}
    assert by_id['rs116555717']['clinvar_variation_id'] == '3305447'
    assert by_id['rs116555717']['clinvar_classification'] == 'Conflicting classifications of pathogenicity'
    assert by_id['rs114053718']['clinvar_classification'] == 'Uncertain significance'
    assert by_id['rs1135638']['clinvar_present'] == 'no'
    assert by_id['rs1135638']['selected_biotype'] == 'nonsense_mediated_decay'
    assert by_id['rs1135638']['mane_consequences'] == 'synonymous_variant'
    assert by_id['rs16994704']['clinvar_classification'] == 'Benign'
    assert by_id['rs16994704']['sift_prediction'] == 'deleterious'
    assert by_id['rs61735928']['rsid'] == 'rs20550'
    assert by_id['rs61735928']['rsid_status'] == 'merged'
    merged = json.loads((ROOT/'raw/dbsnp/rs61735928.json').read_text())
    assert merged['merged_snapshot_data']['merged_into'] == ['20550']
    for request in json.loads((ROOT/'raw/requests.json').read_text()):
        actual = hashlib.sha256((ROOT/'raw'/request['file']).read_bytes()).hexdigest()
        assert actual == request['sha256'], request['file']
    assert sum(r['clinvar_present'] == 'yes' for r in annotations) == 19
    assert sum(r['matches_input_variant'] == 'no' for r in table('clinvar_records_fetched.tsv')) == 5
    print('PASS: 50 unique variants; input order; VEP transcripts; exact ClinVar alleles;')
    print('PASS: cohort-specific gnomAD AF=AC/AN; missing-value semantics; source SHA-256 hashes.')


if __name__ == '__main__':
    main()
