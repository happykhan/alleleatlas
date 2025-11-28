#!/usr/bin/env python3
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p='out/hnsw_sweep.json'
with open(p,'r') as f:
    data=json.load(f)

# compute total_time = build_time + rerank_time + minhash_time (if present)
def total_time(entry):
    return sum(entry.get(k,0) for k in ('build_time','rerank_time','minhash_time'))

for e in data:
    e['total_time']=total_time(e)
    if 'recall_at_k' not in e:
        e['recall_at_k']=float('nan')

# pick entries with candidate_multiplier==1 and candidate_k==5 (fast configs)
fast=[e for e in data if e.get('candidate_multiplier',None)==1 and e.get('candidate_k',None)==5]
fast_sorted=sorted(fast, key=lambda x:(-x['recall_at_k'], x['total_time']))

# pick overall top by recall/time (recall divided by total_time), avoid zero time
ranked=sorted([e for e in data if e.get('total_time',0)>0], key=lambda x: (- (x.get('recall_at_k',0)/(x.get('total_time',1e-9))), -x.get('recall_at_k',0)))

# select top 6 configs
top_by_recall=fast_sorted[:6]
top_by_recalltime=ranked[:6]

os.makedirs('out',exist_ok=True)
with open('out/hnsw_sweep_top_configs.json','w') as f:
    json.dump({'top_by_recall':[{k:v for k,v in e.items() if k in ['dim','M','ef_construction','ef_search','candidate_multiplier','candidate_k','method','recall_at_k','mst_edge_overlap','total_time']} for e in top_by_recall]}, f,indent=2)

with open('out/hnsw_sweep_top_by_recalltime.json','w') as f:
    json.dump({'top_by_recalltime':[{k:v for k,v in e.items() if k in ['dim','M','ef_construction','ef_search','candidate_multiplier','candidate_k','method','recall_at_k','mst_edge_overlap','total_time']} for e in top_by_recalltime]}, f,indent=2)

# plotting recall vs total_time colored by dim
dims=sorted(set(e['dim'] for e in data))
colors={d:plt.cm.viridis(i/len(dims)) for i,d in enumerate(dims)}

plt.figure(figsize=(7,4))
for e in data:
    plt.scatter(e['total_time'], e['recall_at_k'], color=colors[e['dim']], s=18, alpha=0.6)
for i,e in enumerate(top_by_recalltime):
    plt.annotate(f"{i+1}:{e['dim']},ef{e['ef_search']},M{e['M']}", (e['total_time'], e['recall_at_k']))
plt.xlabel('total_time (s)')
plt.ylabel('recall@k')
plt.title('HNSW sweep: recall vs time')
plt.grid(alpha=0.2)
plt.savefig('out/hnsw_sweep_recall_vs_time.png', bbox_inches='tight', dpi=150)
plt.close()

# plotting mst vs total_time
plt.figure(figsize=(7,4))
for e in data:
    plt.scatter(e['total_time'], e.get('mst_edge_overlap',float('nan')), color=colors[e['dim']], s=18, alpha=0.6)
for i,e in enumerate(top_by_recalltime):
    plt.annotate(f"{i+1}:{e['dim']},ef{e['ef_search']},M{e['M']}", (e['total_time'], e.get('mst_edge_overlap',0)))
plt.xlabel('total_time (s)')
plt.ylabel('mst_edge_overlap')
plt.title('HNSW sweep: MST overlap vs time')
plt.grid(alpha=0.2)
plt.savefig('out/hnsw_sweep_mst_vs_time.png', bbox_inches='tight', dpi=150)
plt.close()

print('Wrote: out/hnsw_sweep_top_configs.json, out/hnsw_sweep_top_by_recalltime.json,')
print('       out/hnsw_sweep_recall_vs_time.png, out/hnsw_sweep_mst_vs_time.png')
print('\nTop configs (by recall/time):')
for i,e in enumerate(top_by_recalltime[:8]):
    print(i+1, e.get('method','hnsw'), 'dim', e['dim'], 'M', e['M'], 'ef', e['ef_search'], 'recall', e['recall_at_k'], 'time', round(e['total_time'],4))
