"""
Evaluation pipeline for the thesis run.
Per pair: SpatialSimBench E1-E3, BEASTsim Dice/DSC + DE biological signal,
BP components+summaries, perturbation robustness, interpretability.
Writes each pair to its own JSON shard the moment it finishes, so a worker crash or
pool hang does not lose finished pairs. The final CSVs assemble from the shards on disk.
"""
from __future__ import annotations
import argparse, json, os, platform, sys, time, traceback, glob, re
for _v in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_v, '2')
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError as FTimeout
from pathlib import Path
import numpy as np, pandas as pd
SD = Path(__file__).resolve().parent; sys.path.insert(0, str(SD))
import ssb_common as C, spatialsimbench_metrics as M, perturbations as P, bp_module as BP

def _mem_gb():
    try:
        import resource; return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2),2)
    except Exception: return float("nan")
def _safe(s): return re.sub(r"[^A-Za-z0-9._-]","_",str(s))

def process_pair(row, data_root, corr_sample, max_obs, pert_max_obs, replicates, seed, fig_dir, shard_dir, gene_cap, representative, inner_threads=4):
    t0=time.time(); tag=f'{row.get("platform")}/{row.get("dataset_variant")}/{row.get("simulator")}'
    rng=np.random.default_rng(seed)
    res={"row":row,"tag":tag,"metrics":[],"summaries":[],"components":[],"perturbations":[],"status":"ok"}
    try:
        R=C.read_h5ad(os.path.join(data_root,row["real_h5ad"])); S=C.read_h5ad(os.path.join(data_root,row["simulated_h5ad"]))
        R,S,genes=C.common_genes(R,S)
        def cap(a,n): return a[rng.choice(a.shape[0],n,replace=False)] if (n and a.shape[0]>n) else a
        Rb,Sb=cap(R,max_obs),cap(S,max_obs); XR,XS=C.dense(Rb),C.dense(Sb)
        if gene_cap and XR.shape[1]>gene_cap:
            gi=rng.choice(XR.shape[1],gene_cap,replace=False); XR=XR[:,gi]; XS=XS[:,gi]
        xyR=C.get_coords(Rb,row.get("coordinate_key_real","obsm/spatial")); xyS=C.get_coords(Sb,row.get("coordinate_key_sim","obsm/spatial"))
        labR=C.get_labels(Rb,row.get("label_key_real","obs/Cell_Type")); labS=C.get_labels(Sb,row.get("label_key_sim","obs/cell_type"))
        pl,dv,sm=row.get("platform"),row.get("dataset_variant"),row.get("simulator")
        def add(fam,met,val): res["metrics"].append(dict(platform=pl,dataset_variant=dv,simulator=sm,seed=seed,family=fam,metric=met,value=val,evidence_state="computed"))
        for fam,fn in [("E1_data_properties",lambda:M.data_property_metrics(XR,XS,n_corr_sample=corr_sample,seed=seed)),
                       ("E2_SVG_moran",lambda:M.svg_and_moran_metrics(XR,XS,xyR,xyS)),
                       ("E2_clustering",lambda:M.clustering_metrics(XR,XS,labR,seed=seed)),
                       ("E2_ct_proportion",lambda:M.ct_proportion_metrics(labR,labS)),
                       ("BEASTsim_similarity_DSC",lambda:M.beastsim_similarity(XR,XS,xyR,xyS,labR,labS)),
                       ("BEASTsim_biological_signal_DE",lambda:M.de_signal_preservation(XR,XS,labR,labS))]:
            try:
                for k,v in fn().items(): add(fam,k,v)
            except Exception as e: add(fam,"__error__",str(e))
        comp,copy_sim=M.fast_components(XR,XS,xyR,xyS,labR,labS)
        for k,v in comp.items(): res["components"].append(dict(platform=pl,simulator=sm,dataset_variant=dv,seed=seed,component=k,score=v,evidence_state="computed"))
        bs=BP.named_summaries(comp,copy_sim)
        for k,v in bs.items():
            add("biological_plausibility",k,v); res["summaries"].append(dict(platform=pl,simulator=sm,dataset_variant=dv,seed=seed,summary=k,score=v))
        Rp,Sp=cap(R,pert_max_obs),cap(S,pert_max_obs); XRp,XSp=C.dense(Rp),C.dense(Sp)
        if gene_cap and XRp.shape[1]>gene_cap:
            gi=rng.choice(XRp.shape[1],gene_cap,replace=False); XRp=XRp[:,gi]; XSp=XSp[:,gi]
        xyRp=C.get_coords(Rp,row.get("coordinate_key_real","obsm/spatial")); xySp=C.get_coords(Sp,row.get("coordinate_key_sim","obsm/spatial"))
        labRp=C.get_labels(Rp,row.get("label_key_real","obs/Cell_Type")); labSp=C.get_labels(Sp,row.get("label_key_sim","obs/cell_type"))
        c0,cs0=M.fast_components(XRp,XSp,xyRp,xySp,labRp,labSp); s0=BP.named_summaries(c0,cs0)
        from concurrent.futures import ThreadPoolExecutor as _TPE
        def _ptask(pr):
            pert,rep=pr; base=P._sev(pert)[1]; exp=P.EXPECTED.get(base,""); out=[]
            try:
                Xp,xyp,lp=P.apply_perturbation(XSp,xySp,labSp,pert,np.random.default_rng(seed+rep),real_X=XRp)
                cp,csp=M.fast_components(XRp,Xp,xyRp,xyp,labRp,lp); sp=BP.named_summaries(cp,csp)
                for s in ["Balanced Realism","Spatial Tissue Fidelity","Novelty-Gated Plausibility"]:
                    out.append(dict(platform=pl,simulator=sm,dataset_variant=dv,seed=seed,perturbation=pert,replicate=rep,summary=s,score=sp[s],delta=sp[s]-s0[s],expected_damaged_layer=exp,evidence_state="computed"))
            except Exception as e:
                out.append(dict(platform=pl,simulator=sm,dataset_variant=dv,seed=seed,perturbation=pert,replicate=rep,summary="error",score=float("nan"),delta=float("nan"),expected_damaged_layer=exp,evidence_state=f"failed:{e}"))
            return out
        _tasks=[(pert,rep) for pert in P.BATTERY if pert!="baseline" for rep in range(replicates)]
        with _TPE(max_workers=max(1,inner_threads)) as _tp:
            for rows in _tp.map(_ptask,_tasks): res["perturbations"].extend(rows)
        if representative and labR is not None:
            try:
                fd=Path(fig_dir); fd.mkdir(parents=True,exist_ok=True)
                MR,cats=M._ct_cooccurrence(XR,xyR,labR); MS,_=M._ct_cooccurrence(XS,xyS,labS if labS is not None else labR)
                pd.DataFrame(MR,index=cats,columns=cats).to_csv(fd/f"{_safe(tag)}_neighbourhood_real.csv")
                pd.DataFrame(MS,index=cats,columns=cats).to_csv(fd/f"{_safe(tag)}_neighbourhood_sim.csv")
            except Exception: pass
        res["compute_seconds"]=round(time.time()-t0,2)
        add("E3_scalability","compute_seconds",res["compute_seconds"]); add("E3_scalability","peak_memory_gb",_mem_gb())
    except Exception as e:
        res["status"]="failed"; res["error"]=str(e); res["traceback"]=traceback.format_exc()
    try: json.dump(res,open(os.path.join(shard_dir,_safe(tag)+f"__seed{seed}"+".json"),"w"))
    except Exception: pass
    return {"tag":tag,"status":res["status"],"compute_seconds":res.get("compute_seconds"),"n_metrics":len(res["metrics"]),"n_pert":len(res["perturbations"]),"error":res.get("error")}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest",required=True); ap.add_argument("--data-root",default=".")
    ap.add_argument("--output-root",default="outputs"); ap.add_argument("--workers",type=int,default=os.cpu_count())
    ap.add_argument("--corr-sample",type=int,default=2000); ap.add_argument("--max-obs",type=int,default=8000)
    ap.add_argument("--pert-max-obs",type=int,default=3000); ap.add_argument("--replicates",type=int,default=3)
    ap.add_argument("--gene-cap",type=int,default=2500); ap.add_argument("--global-timeout",type=int,default=1800)
    ap.add_argument("--seed",type=int,default=0)
    ap.add_argument("--seeds",type=str,default="",help="comma-separated seeds for multi-seed replication of the main run; overrides --seed when set")
    a=ap.parse_args()
    seeds=[int(x) for x in a.seeds.split(",") if x.strip()!=""] or [a.seed]
    out=Path(a.output_root); out.mkdir(parents=True,exist_ok=True)
    fig_dir=out/"figure_sources"; shard_dir=out/"_shards"; shard_dir.mkdir(exist_ok=True)
    man=pd.read_csv(a.manifest).to_dict("records"); seen=set()
    for r in man: r["_rep"]=r.get("platform") not in seen; seen.add(r.get("platform"))
    print(f"[thesis] {len(man)} pairs | seeds={seeds} | replicates={a.replicates} | max_obs={a.max_obs} | gene_cap={a.gene_cap}",flush=True)
    t0=time.time(); cw=min(a.workers,len(man),32)
    inner=max(2,min(8,(os.cpu_count() or 8)//max(1,cw)))
    print(f'[thesis] outer_workers={cw} inner_threads={inner} (~{cw*inner} cores engaged)',flush=True)
    import multiprocessing as mp
    ctx=mp.get_context('spawn')
    for si,seed in enumerate(seeds):
        if len(seeds)>1: print(f"[thesis] === seed {seed} ({si+1}/{len(seeds)}) ===",flush=True)
        try:
            with ProcessPoolExecutor(max_workers=cw, mp_context=ctx) as ex:
                futs={ex.submit(process_pair,r,a.data_root,a.corr_sample,a.max_obs,a.pert_max_obs,a.replicates,seed,str(fig_dir),str(shard_dir),a.gene_cap,(r["_rep"] and si==0),inner):r for r in man}
                try:
                    for f in as_completed(futs,timeout=a.global_timeout):
                        s=f.result(); print(f"[{s['status']:4}] seed{seed} {s['tag']}: {s['n_metrics']} metrics, {s['n_pert']} pert rows, {s.get('compute_seconds')}s"+(f"  ERR={s['error']}" if s['status']!='ok' else ""),flush=True)
                except FTimeout: print(f"[thesis] global timeout on seed {seed}; assembling shards",flush=True)
                ex.shutdown(wait=False,cancel_futures=True)
        except Exception as e: print(f"[thesis] pool error on seed {seed} ({e}); assembling shards",flush=True)
    results=[]
    for p in glob.glob(str(shard_dir/"*.json")):
        try: results.append(json.load(open(p)))
        except Exception: pass
    def col(k): return pd.DataFrame([x for r in results for x in r.get(k,[])])
    col("metrics").to_csv(out/"metric_long.csv",index=False)
    col("summaries").to_csv(out/"named_plausibility_baseline.csv",index=False)
    col("components").to_csv(out/"component_evidence.csv",index=False)
    pr=col("perturbations"); pr.to_csv(out/"perturbation_response.csv",index=False)
    if len(pr): pr[pr.summary=="Novelty-Gated Plausibility"].groupby(["perturbation","expected_damaged_layer"])["delta"].median().reset_index().to_csv(out/"expected_vs_observed.csv",index=False)
    json.dump({"python":platform.python_version(),"numpy":np.__version__,"cpu_count":os.cpu_count(),"workers":cw},open(out/"environment.json","w"),indent=2)
    json.dump({"pairs":len(man),"shards_written":len(results),"ok":sum(r.get("status")=="ok" for r in results),"failed":sum(r.get("status")!="ok" for r in results),"elapsed_seconds":round(time.time()-t0,1),"failures":[{"tag":r["tag"],"error":r.get("error")} for r in results if r.get("status")!="ok"]},open(out/"run_summary.json","w"),indent=2)
    json.dump({"deliverable_1_spatialsimbench":["metric_long.csv (E1_data_properties,E2_*,E3_scalability)"],"beastsim_native_measures":["metric_long.csv (BEASTsim_similarity_DSC, BEASTsim_biological_signal_DE)"],"deliverable_3_bp_prototype":["named_plausibility_baseline.csv","component_evidence.csv"],"deliverable_3_robustness":["perturbation_response.csv","expected_vs_observed.csv"],"deliverable_3_interpretability":["figure_sources/*_neighbourhood_{real,sim}.csv"]},open(out/"contract_coverage.json","w"),indent=2)
    print(f"[thesis] DONE in {time.time()-t0:.1f}s | shards={len(results)} -> {out}",flush=True); return 0

if __name__=="__main__":
    raise SystemExit(main())
