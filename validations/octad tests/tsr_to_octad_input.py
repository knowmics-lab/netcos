# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 12:19:52 2025

@author: L-F-S
"""

import pandas as pd

signature_filename = 'ipf_signature_gene_id.csv'
signature_octad_format_filename = 'ipf.oct'

df = pd.read_csv(signature_filename, sep=';')
df.rename(columns={
    'DE_log2_FC':'log2FoldChange', 'adj.p.value':'padj',
    'p.value':'pval','gene':'symbol','gene_id':'gene'}
    , inplace=True)

df[['log2FoldChange','pval','padj','gene','symbol']].to_csv(signature_octad_format_filename, index=False, sep=',')
