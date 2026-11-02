#!/usr/bin/env python3
"""Build annotation tables offline from the saved, unmodified API responses."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET
from fetch_annotations import VEP_OPTIONS

ROOT = Path(__file__).resolve().parents[1]
RULE = 'VEP_most_severe_then_MANE_Select_then_canonical_then_protein_coding_then_transcript_id'


def read_json(path):
    return json.loads((ROOT / path).read_text())


def read_vcf():
    rows = []
    for line in (ROOT / 'input/student_06.vcf').read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        chrom, pos, ident, ref, alt, qual, filt, info = line.split('\t')
        assert len(ref) == len(alt) == 1 and ref in 'ACGT' and alt in 'ACGT' and ref != alt
        rows.append(dict(chrom=chrom, pos=int(pos), original_id=ident, ref=ref, alt=alt,
                         variant_id=f'{chrom}-{pos}-{ref}-{alt}',
                         substitution_type='transition' if {ref, alt} in [{'A', 'G'}, {'C', 'T'}]
                         else 'transversion'))
    assert len(rows) == len({r['variant_id'] for r in rows}) == 50
    return rows


def write_table(name, rows):
    assert rows
    path = ROOT / 'results' / name
    path.parent.mkdir(exist_ok=True)
    with path.open('w', newline='', encoding='utf-8-sig' if path.suffix == '.csv' else 'utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter=',' if path.suffix == '.csv' else '\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: 'NA' if v is None or v == '' else v for k, v in row.items()})


def select_transcript(vep, alt):
    """Preserve VEP's global most-severe term, then choose a supporting transcript.

    MANE/canonical are tie-breakers only, never a reason to silently replace the
    global term with a less severe term. All transcript rows are exported.
    """
    worst = vep['most_severe_consequence']
    candidates = [t for t in vep.get('transcript_consequences', [])
                  if worst in t['consequence_terms'] and t['variant_allele'] == alt]
    if not candidates:
        # This assignment has only transcript-level worst consequences. Abort if
        # this changes, rather than assigning an unrelated transcript or gene.
        raise ValueError('No transcript supports worst consequence: '+vep['id'])
    return min(candidates, key=lambda t: (not bool(t.get('mane_select')),
               not bool(t.get('canonical')), t.get('biotype') != 'protein_coding', t['transcript_id']))


def clinvar_index():
    """Index only classified SINGLE alleles by exact assembly+chrom+POS+REF+ALT."""
    root = ET.parse(ROOT / 'raw/clinvar.xml').getroot()
    index = defaultdict(list)
    audit = []
    for archive in root.findall('.//VariationArchive'):
        allele = archive.find('./ClassifiedRecord/SimpleAllele')
        if allele is None:
            raise ValueError('Review unexpected non-simple ClinVar record '+archive.get('VariationID'))
        loc = allele.find('./Location/SequenceLocation[@Assembly="GRCh38"]')
        if loc is None:
            raise ValueError('No GRCh38 location '+archive.get('VariationID'))
        key = (loc.get('Chr'), int(loc.get('positionVCF')), loc.get('referenceAlleleVCF'), loc.get('alternateAlleleVCF'))
        cs = archive.find('./ClassifiedRecord/Classifications/GermlineClassification')
        classification = cs.findtext('Description') if cs is not None else None
        record = dict(clinvar_variation_id=archive.get('VariationID'),
                      clinvar_accession=archive.get('Accession')+'.'+archive.get('Version'),
                      clinvar_classification=classification,
                      clinvar_review_status=cs.findtext('ReviewStatus') if cs is not None else None,
                      clinvar_last_evaluated=cs.get('DateLastEvaluated') if cs is not None else None,
                      clinvar_conflict_explanation=cs.findtext('Explanation') if cs is not None else None,
                      clinvar_conditions=' | '.join(sorted(set(
                          e.text for e in archive.findall('./ClassifiedRecord/RCVList/RCVAccession/ClassifiedConditionList/ClassifiedCondition')
                          if e.text))),
                      clinvar_record_updated=archive.get('DateLastUpdated'),
                      clinvar_url='https://www.ncbi.nlm.nih.gov/clinvar/variation/'+archive.get('VariationID')+'/')
        index[key].append(record)
        audit.append(dict(chrom=key[0], pos=key[1], ref=key[2], alt=key[3], **record))
    return index, audit


def gnomad_result(row):
    response = read_json('raw/gnomad/'+row['variant_id']+'.json')
    if response.get('errors'):
        assert all('Variant not found' in e['message'] for e in response['errors']), response
    v = (response.get('data') or {}).get('variant')
    result = dict(gnomad_status='not_found', gnomad_dataset='gnomad_r4', gnomad_af=None,
                  gnomad_af_source=None, gnomad_exome_af=None, gnomad_exome_ac=None,
                  gnomad_exome_an=None, gnomad_exome_filters=None, gnomad_genome_af=None,
                  gnomad_genome_ac=None, gnomad_genome_an=None, gnomad_genome_filters=None,
                  gnomad_url='https://gnomad.broadinstitute.org/variant/'+row['variant_id']+'?dataset=gnomad_r4')
    if v is None:
        if not response.get('errors') and 'variant' not in response.get('data', {}):
            raise ValueError('Malformed gnomAD response: '+row['variant_id'])
        return result
    assert v['variant_id'] == row['variant_id'] and v['reference_genome'] == 'GRCh38'
    result['gnomad_status'] = 'present_no_unfiltered_af'
    for source in ['exome', 'genome']:
        counts = v.get(source)
        if counts is None:
            continue
        ac, an, af = counts['ac'], counts['an'], counts['af']
        assert 0 <= ac <= an
        if an > 0:
            assert af is not None and abs(af-ac/an) < 1e-10
        result.update({f'gnomad_{source}_af': af, f'gnomad_{source}_ac': ac,
                       f'gnomad_{source}_an': an,
                       f'gnomad_{source}_filters': 'PASS' if counts['filters'] == [] else '|'.join(counts['filters'])})
    # Keep cohorts separate: sample overlap prevents naive AC/AN summation.
    for source in ['genome', 'exome']:
        counts = v.get(source)
        if counts and counts['an'] > 0 and counts['af'] is not None and counts['filters'] == []:
            result.update(gnomad_status='present', gnomad_af=counts['af'], gnomad_af_source=source)
            break
    return result


def main():
    rows = read_vcf()
    vep_list = read_json('raw/vep.json')
    assert len(vep_list) == 50
    vep_by_key = {}
    for v in vep_list:
        assert v['assembly_name'] == 'GRCh38' and len(v['allele_string'].split('/')) == 2
        ref, alt = v['allele_string'].split('/')
        key = (v['seq_region_name'], v['start'], ref, alt)
        assert key not in vep_by_key
        vep_by_key[key] = v
    cv_index, cv_audit = clinvar_index()
    input_ids = {(r['chrom'], r['pos'], r['ref'], r['alt']): r['original_id'] for r in rows}
    for record in cv_audit:
        matched_id = input_ids.get((record['chrom'], record['pos'], record['ref'], record['alt']))
        record['matches_input_variant'] = 'yes' if matched_id else 'no'
        record['matched_original_id'] = matched_id
    annotations, transcript_rows = [], []
    for row in rows:
        key = (row['chrom'], row['pos'], row['ref'], row['alt'])
        v = vep_by_key[key]
        assert v['id'] == row['original_id']
        t = select_transcript(v, row['alt'])
        colocated = [c for c in v.get('colocated_variants', [])
                     if c['start'] == row['pos'] and c.get('strand', 1) == 1
                     and row['ref'] in c.get('allele_string', '').split('/')
                     and row['alt'] in c.get('allele_string', '').split('/')]
        rsids = sorted(set(c['id'] for c in colocated if c['id'].startswith('rs')))
        assert rsids, 'No allele-matching current rsID for '+row['original_id']
        rsid_status = 'unchanged'
        if row['original_id'] not in rsids:
            resolution = read_json('raw/dbsnp/'+row['original_id']+'.json')
            targets = ['rs'+ident for ident in resolution['merged_snapshot_data']['merged_into']]
            assert set(targets) & set(rsids)
            rsid_status = 'merged'
        mane = sorted([x for x in v.get('transcript_consequences', []) if x.get('mane_select')],
                      key=lambda x: x['transcript_id'])
        annotation = dict(**row, gene=t.get('gene_symbol') or t.get('gene_id'),
            all_genes='|'.join(sorted(set(x.get('gene_symbol') or x.get('gene_id')
                               for x in v.get('transcript_consequences', [])
                               if x.get('gene_symbol') or x.get('gene_id')))),
            most_severe_consequence=v['most_severe_consequence'], selected_transcript=t['transcript_id'],
            selected_biotype=t['biotype'], selected_transcript_consequences='|'.join(t['consequence_terms']),
            transcript_selection_rule=RULE, selected_mane_select=t.get('mane_select'),
            selected_canonical=t.get('canonical', 0), impact=t.get('impact'), hgvsc=t.get('hgvsc'),
            hgvsp=t.get('hgvsp'), sift_prediction=t.get('sift_prediction'), sift_score=t.get('sift_score'),
            polyphen_prediction=t.get('polyphen_prediction'), polyphen_score=t.get('polyphen_score'),
            mane_transcripts='|'.join(x['transcript_id'] for x in mane),
            mane_consequences='|'.join(sorted(set(c for x in mane for c in x['consequence_terms']))),
            mane_hgvsp='|'.join(x.get('hgvsp', 'NA') for x in mane),
            rsid='|'.join(rsids), rsid_status=rsid_status,
            transcript_count=len(v.get('transcript_consequences', [])),
            **gnomad_result(row))
        hits = cv_index.get(key, [])
        annotation['clinvar_present'] = 'yes' if hits else 'no'
        cvfields = ['clinvar_variation_id', 'clinvar_accession', 'clinvar_classification',
                    'clinvar_review_status', 'clinvar_last_evaluated', 'clinvar_conflict_explanation',
                    'clinvar_conditions', 'clinvar_record_updated', 'clinvar_url']
        for field in cvfields:
            values = [h[field] for h in hits if h[field]]
            annotation[field] = ' || '.join(values) if values else None
        annotations.append(annotation)
        for tc in v.get('transcript_consequences', []):
            transcript_rows.append(dict(original_id=row['original_id'], variant_id=row['variant_id'],
                gene=tc.get('gene_symbol') or tc.get('gene_id'),
                transcript_id=tc['transcript_id'], biotype=tc.get('biotype'),
                consequences='|'.join(tc['consequence_terms']),
                is_selected=int(tc['transcript_id'] == t['transcript_id']), mane_select=tc.get('mane_select'),
                canonical=tc.get('canonical', 0), hgvsc=tc.get('hgvsc'), hgvsp=tc.get('hgvsp'),
                sift_prediction=tc.get('sift_prediction'), sift_score=tc.get('sift_score'),
                polyphen_prediction=tc.get('polyphen_prediction'), polyphen_score=tc.get('polyphen_score')))
    write_table('annotations.tsv', annotations)
    write_table('annotations.csv', annotations)
    write_table('vep_transcripts.tsv', transcript_rows)
    write_table('clinvar_records_fetched.tsv', cv_audit)
    write_table('clinvar.tsv', [{k: v for k, v in a.items()
                if k in ['original_id', 'variant_id'] or k.startswith('clinvar_')} for a in annotations])
    write_table('gnomad.tsv', [{k: v for k, v in a.items()
                if k in ['original_id', 'variant_id'] or k.startswith('gnomad_')} for a in annotations])
    summary = dict(variant_count=len(rows), chromosome_counts=dict(Counter(r['chrom'] for r in rows)),
                   substitution_counts=dict(Counter(r['substitution_type'] for r in rows)),
                   ti_tv_ratio=sum(r['substitution_type'] == 'transition' for r in rows) /
                               sum(r['substitution_type'] == 'transversion' for r in rows),
                   consequences=dict(Counter(a['most_severe_consequence'] for a in annotations)),
                   clinvar_present=sum(a['clinvar_present'] == 'yes' for a in annotations),
                   clinvar_absent=sum(a['clinvar_present'] == 'no' for a in annotations),
                   clinvar_classifications=dict(Counter(a['clinvar_classification'] for a in annotations
                                                       if a['clinvar_present'] == 'yes')),
                   gnomad_status_counts=dict(Counter(a['gnomad_status'] for a in annotations)),
                   gnomad_primary_af_source=dict(Counter(a['gnomad_af_source'] for a in annotations)),
                   transcript_row_count=len(transcript_rows),
                   selected_non_mane=[a['original_id'] for a in annotations if not a['selected_mane_select']],
                   worst_differs_from_mane=[a['original_id'] for a in annotations
                       if a['most_severe_consequence'] not in a['mane_consequences'].split('|')],
                   merged_rsids={a['original_id']: a['rsid'] for a in annotations if a['rsid_status']=='merged'})
    (ROOT / 'results/statistics.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n')
    metadata = dict(assembly='GRCh38', species='homo_sapiens', variant=6,
                    input_sha256=hashlib.sha256((ROOT/'input/student_06.vcf').read_bytes()).hexdigest(),
                    ensembl_release=read_json('raw/ensembl_release.json'),
                    ensembl_rest=read_json('raw/ensembl_rest.json'),
                    vep_options=VEP_OPTIONS,
                    transcripts='Ensembl, all returned transcripts; no pick/per_gene filtering',
                    selection_rule=RULE, gnomad_dataset='gnomad_r4',
                    gnomad_af_rule='Unfiltered genome AF if AN>0, else unfiltered exome AF; never combine cohorts',
                    clinvar_search='Union of rsID (All Fields) and chromosome+GRCh38 position (CHR+CPOS/chrpos38)',
                    clinvar_match='ClassifiedRecord/SimpleAllele/Location/SequenceLocation: exact GRCh38, Chr, positionVCF, referenceAlleleVCF, alternateAlleleVCF',
                    clinvar_classification='VCV aggregate GermlineClassification; retain conflicts verbatim',
                    rsid_resolution='Original IDs preserved; current allele-matching rsIDs from VEP; merges independently verified using NCBI RefSNP raw/dbsnp/',
                    missing_values='NA means unavailable/not applicable; never replaces missing AF with zero',
                    request_provenance='raw/requests.json')
    (ROOT / 'results/annotation_settings.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
