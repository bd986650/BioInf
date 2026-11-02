#!/usr/bin/env python3
"""Download public GRCh38 annotations; Python 3 standard library only.

Responses are cached verbatim. Remove only the intended cached response to
refresh it. A failed request never becomes a "variant absent" annotation.
"""
import argparse
import hashlib
import http.client
import json
from pathlib import Path
import time
from datetime import datetime, timezone
import urllib.error
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'raw'
VEP_OPTIONS = dict(canonical=1, mane=1, hgvs=1, protein=1, numbers=1,
                   transcript_version=1, variant_class=1)
GNOMAD_QUERY = '''query Variant($id: String!) {
  variant(variantId: $id, dataset: gnomad_r4) {
    variant_id reference_genome
    exome { ac an af filters }
    genome { ac an af filters }
  }
}'''


class RetryableServiceError(Exception):
    pass


def variants():
    return [line.split('\t') for line in
            (ROOT / 'input/student_06.vcf').read_text().splitlines()
            if line and not line.startswith('#')]


def fetch(filename, url, payload=None, binary=False):
    """Save successful response bytes and a non-secret reproducibility log."""
    path = RAW / filename
    if path.exists():
        body = path.read_bytes()
        return body if binary else json.loads(body)
    data = json.dumps(payload).encode() if payload is not None else None
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, data=data, headers={
                'Content-Type': 'application/json',
                'Accept': 'application/xml' if binary else 'application/json',
                'User-Agent': 'BioInf-Assignment4-variant6/1.0'})
            with urllib.request.urlopen(req, timeout=120) as response:
                body = response.read()
                response_headers = dict(response.headers)
            if not binary:
                result = json.loads(body)
                if isinstance(result, dict) and result.get('error'):
                    raise RuntimeError(result['error'])
                if isinstance(result, dict) and result.get('errors'):
                    messages = [e['message'] for e in result['errors']]
                    if any('overloaded' in m.lower() or 'rate limit' in m.lower() for m in messages):
                        raise RetryableServiceError('; '.join(messages))
                    if not all('Variant not found' in m for m in messages):
                        raise RuntimeError(messages)
            else:
                result = body
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(body)
            logpath = RAW / 'requests.json'
            log = json.loads(logpath.read_text()) if logpath.exists() else []
            log.append(dict(file=filename, url=url, method='POST' if data else 'GET',
                            payload=payload, retrieved_at_utc=datetime.now(timezone.utc).isoformat(),
                            response_headers=response_headers,
                            sha256=hashlib.sha256(body).hexdigest()))
            logpath.write_text(json.dumps(log, indent=2, ensure_ascii=False)+'\n')
            print('Saved', filename, len(body), 'bytes', flush=True)
            return result
        except (urllib.error.URLError, TimeoutError, ConnectionError,
                http.client.RemoteDisconnected, http.client.IncompleteRead, RetryableServiceError) as e:
            if attempt == 4:
                raise
            print('Retry', attempt+1, filename, str(e), flush=True)
            time.sleep(min(2 ** attempt, 15))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', choices=['all', 'vep', 'clinvar', 'gnomad'], default='all')
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    records = variants()
    assert len(records) == 50
    if args.source in ('all', 'vep'):
        fetch('ensembl_release.json', 'https://rest.ensembl.org/info/data?content-type=application/json')
        fetch('ensembl_rest.json', 'https://rest.ensembl.org/info/rest?content-type=application/json')
        vep = fetch('vep.json', 'https://rest.ensembl.org/vep/homo_sapiens/region?'+
              urllib.parse.urlencode(VEP_OPTIONS), {'variants': [' '.join(r) for r in records]})
        for v in vep:
            if v['id'] not in [c['id'] for c in v.get('colocated_variants', [])]:
                fetch('dbsnp/'+v['id']+'.json',
                      'https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/'+v['id'][2:])
    if args.source in ('all', 'clinvar'):
        # Search by original rsIDs AND independently by GRCh38 chromosome/position.
        # Exact REF/ALT matching is done against the downloaded VCV XML later.
        base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        fetch('clinvar_database_info.json', base+'einfo.fcgi?db=clinvar&retmode=json')
        terms = {
            'rs': ' OR '.join(r[2]+'[All Fields]' for r in records),
            'coordinates': ' OR '.join('('+r[0]+'[CHR] AND '+r[1]+'[CPOS])' for r in records)
        }
        ids = set()
        for mode, term in terms.items():
            result = fetch('clinvar_search_'+mode+'.json', base+'esearch.fcgi?'+
                           urllib.parse.urlencode(dict(db='clinvar', term=term, retmode='json', retmax=1000)))
            search = result['esearchresult']
            # phrasesnotfound lists absent rsIDs/positions in an OR query;
            # this is expected negative evidence, not a transport/query error.
            errors = search.get('errorlist', {})
            if errors.get('fieldsnotfound') or search.get('warninglist') or search.get('ERROR'):
                raise RuntimeError('ClinVar query warning/error: '+json.dumps(search))
            if int(search['count']) > len(search['idlist']):
                raise RuntimeError('ClinVar result is truncated')
            ids.update(search['idlist'])
            time.sleep(0.4)  # NCBI: <=3 requests/s without API key.
        if ids:
            idstr = ','.join(sorted(ids, key=int))
            fetch('clinvar_summary.json', base+'esummary.fcgi?'+
                  urllib.parse.urlencode(dict(db='clinvar', id=idstr, retmode='json')))
            time.sleep(0.4)
            xml = fetch('clinvar.xml', base+'efetch.fcgi?'+urllib.parse.urlencode(dict(
                db='clinvar', id=idstr, rettype='vcv', is_variationid='true', retmode='xml')), binary=True)
            if b'VariationArchive' not in xml:
                raise RuntimeError('ClinVar XML contains no VariationArchive records')
        else:
            raise RuntimeError('Unexpectedly no ClinVar records for the entire input')
    if args.source in ('all', 'gnomad'):
        for r in records:
            vid = '-'.join([r[0], r[1], r[3], r[4]])
            cached = (RAW / ('gnomad/'+vid+'.json')).exists()
            result = fetch('gnomad/'+vid+'.json', 'https://gnomad.broadinstitute.org/api',
                           {'query': GNOMAD_QUERY, 'variables': {'id': vid}})
            if result.get('errors'):
                # A known variant-not-found message is a valid negative response;
                # every other API error must be resolved rather than hidden.
                messages = [e['message'] for e in result['errors']]
                if not all('Variant not found' in m for m in messages):
                    raise RuntimeError(messages)
            if not cached:
                time.sleep(1.0)  # Public gnomAD service: sequential, limited request rate.


if __name__ == '__main__':
    main()
