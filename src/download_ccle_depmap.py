#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download the CCLE baseline expression and cell-line metadata of a DepMap
Public release from its figshare article (the DepMap portal API answers
automated requests with a verification page; figshare is DepMap's own
public distribution of each release).

Release: DepMap 24Q4 Public (figshare article 27993248), the most recent
full release on figshare.

Files (into conf.CCLE_DIR / DEPMAP_RELEASE):
- OmicsExpressionProteinCodingGenesTPMLogp1.csv
    rows = DepMap models (ACH-...), columns = 'SYMBOL (EntrezID)',
    values = log2(TPM + 1), protein-coding genes.
- Model.csv
    model metadata (CCLE name, stripped cell line name, lineage, ...).

Each file is checked against the md5 published by figshare; a file already
present with the right md5 is not downloaded again.

Usage
-----
    python download_ccle_depmap.py
"""

import hashlib
import json
import urllib.request

from conf import CCLE_DIR

DEPMAP_RELEASE = 'DepMap_24Q4'
FIGSHARE_ARTICLE = 27993248
FILES = ['OmicsExpressionProteinCodingGenesTPMLogp1.csv', 'Model.csv']


def md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def article_files(article_id):
    url = f'https://api.figshare.com/v2/articles/{article_id}/files?page_size=200'
    with urllib.request.urlopen(url, timeout=60) as r:
        return {f['name']: f for f in json.load(r)}


if __name__ == '__main__':
    out_dir = CCLE_DIR / DEPMAP_RELEASE
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = article_files(FIGSHARE_ARTICLE)
    for name in FILES:
        f = meta[name]
        expected = f.get('computed_md5') or f.get('supplied_md5')
        dest = out_dir / name
        if dest.exists() and md5(dest) == expected:
            print(f'{name}: present, md5 ok')
            continue
        print(f'{name}: downloading {f["size"] / 1e6:.0f} MB')
        urllib.request.urlretrieve(f['download_url'], dest)
        got = md5(dest)
        assert got == expected, f'{name}: md5 {got} != {expected}'
        print(f'{name}: md5 ok -> {dest}')
