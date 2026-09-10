"""Export a derived scoring catalog with exact backbone text and path-free metadata."""
import argparse,csv,gzip,hashlib,io,json
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()

def main():
 p=argparse.ArgumentParser();p.add_argument('--integrated',type=Path,required=True);p.add_argument('--backbone',type=Path,required=True);p.add_argument('--five-mc',type=Path,required=True);p.add_argument('--five-hmc',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
 with a.backbone.open() as f:
  reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;legacy=list(reader)
 for src in [a.five_mc,a.five_hmc]:
  with src.open() as f:
   for row in csv.DictReader(f,delimiter='\t'):
    assert row['paired_donors']=='4' and row['paired_min_support'] in ['3','3.0'], 'Expected four-donor paired source'
 dest=a.output_dir/'brain_cereb.paper_backbone_plus_5mc_5hmc.four_donors.ctdmr.tsv.gz'
 assert not dest.exists(),dest
 counts={};legacy_i=0;seen=set()
 with gzip.open(a.integrated,'rt') as inp,dest.open('wb') as raw:
  with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as zipped:
   with io.TextIOWrapper(zipped,newline='') as out:
    writer=csv.DictWriter(out,fieldnames=fields+['modification','evidence_id'],delimiter='\t',lineterminator='\n');writer.writeheader()
    for row in csv.DictReader(inp,delimiter='\t'):
     mod=row['modification'];assert mod in ['modifiedC','5mC','5hmC'];record={k:row.get(k,'') for k in fields}
     if mod=='modifiedC':
      original=legacy[legacy_i];assert all(row[k]==original[k] for k in ['chr','start','end']);record=original;legacy_i+=1
     key=(record['chr'],record['start'],record['end'],mod);assert key not in seen;seen.add(key)
     record=record|{'modification':mod,'evidence_id':row['evidence_id']}
     assert not any(v.startswith('/') for v in record.values())
     writer.writerow(record);counts[mod]=counts.get(mod,0)+1
 assert legacy_i==len(legacy)==152961
 manifest={'catalog':dest.name,'sha256':sha(dest),'rows':sum(counts.values()),'rows_by_modification':counts,'donors':['bcontrol1','bcontrol2','bmsa','bparkinsons'],'heldout_donor':None,'reference_build':'GRCh38','policy':'Paper backbone preserved; add ONT markers only outside backbone; separate m and h scoring rows','inputs':[{'filename':x.name,'sha256':sha(x)} for x in [a.backbone,a.five_mc,a.five_hmc]],'paired_filter':{'donors':4,'min_support':3,'support_effect':0.30,'min_effect':0.15,'median_effect':0.40},'validation_scope':'Catalog structure and exact original backbone fields checked; benchmark uses a separate bcontrol1-held-out three-donor atlas'}
 (a.output_dir/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
