import pandas as pd
import numpy as np
import mygene

ogee_essentiality =pd.read_csv('/proj/yunligrp/users/tianyou/gRNA/OGEE/9606_all.txt',sep='\t')
ogee_essentiality = ogee_essentiality[['locus','essential','dataType']]
ogee_essentiality = ogee_essentiality.groupby(['locus','essential']).count()
ogee_essentiality['sum'] = ogee_essentiality['dataType'].groupby(['locus']).transform('sum')
ogee_essentiality['OGEE_prop_Essential'] = ogee_essentiality['dataType']/ogee_essentiality['sum']
ogee_essentiality = ogee_essentiality[['dataType','OGEE_prop_Essential']].unstack(fill_value=0)
ogee_essentiality = pd.DataFrame(ogee_essentiality.to_records()).iloc[:,[0,3]]
ogee_essentiality.columns = ['locus', 'OGEE_prop_Essential']

## change gene IDs to names
ens2entrez = ogee_essentiality['locus'].unique()
mg = mygene.MyGeneInfo()
out = mg.querymany(ens2entrez, scopes='ensembl.gene', fields='entrezgene', species='human', as_dataframe=True)
out = out['entrezgene'].reset_index()
out.columns = ['locus', 'geneId']
out = out.drop_duplicates(subset=['geneId'],keep=False)
ogee_essentiality = pd.merge(ogee_essentiality,out,how='left',on='locus')
ogee_essentiality.to_csv("/proj/yunligrp/users/tianyou/gRNA/OGEE/OGEE_essentiality.txt",sep="\t")

