"""Magnitude-stratified candidate density diagnostic; not survey completeness."""
import argparse,csv,json,math
from pathlib import Path
RED_AREA=0.31462445722188803
TOTAL_AREA=0.528823684998808
EDGES=(-math.inf,23.,24.,25.,26.,26.2,math.inf)

def rate_ratio(red,outside,red_area=RED_AREA,total_area=TOTAL_AREA):
    if not 0<red_area<total_area or min(red,outside)<0: raise ValueError("invalid areas/counts")
    if red==0 or outside==0: return {"ratio":None,"poisson_wald_95pct":None}
    rr=(red/red_area)/(outside/(total_area-red_area))
    h=1.96*math.sqrt(1/red+1/outside)
    return {"ratio":rr,"poisson_wald_95pct":[rr*math.exp(-h),rr*math.exp(h)]}

def fisher_two_sided(table):
    a,b=table[0];c,d=table[1]
    if any(v<0 or int(v)!=v for v in (a,b,c,d)):raise ValueError("invalid table")
    n=a+b+c+d;row=a+b;col=a+c;den=math.comb(n,row)
    prob=lambda k:math.comb(col,k)*math.comb(n-col,row-k)/den
    observed=prob(a)
    return sum(prob(k) for k in range(max(0,row-(n-col)),min(row,col)+1)
               if prob(k)<=observed*(1+1e-12))

def proximity_groups(rows,max_separation_arcsec=0.5):
    parent=list(range(len(rows)))
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]]
            i=parent[i]
        return i
    for i in range(len(rows)):
        for j in range(i+1,len(rows)):
            dx=(rows[i]["ra"]-rows[j]["ra"])*3600*math.cos(math.radians((rows[i]["dec"]+rows[j]["dec"])/2))
            dy=(rows[i]["dec"]-rows[j]["dec"])*3600
            if dx*dx+dy*dy<=max_separation_arcsec**2:parent[root(i)]=root(j)
    groups={}
    for i,row in enumerate(rows):groups.setdefault(root(i),[]).append(row)
    if any(len({r["red"] for r in group})>1 for group in groups.values()):
        raise ValueError("mixed red-epoch coverage in proximity group")
    return [{"mag":sum(r["mag"] for r in group)/len(group),"red":group[0]["red"]}
            for group in groups.values()]

def stratification(data,red_area,total_area):
    bright=[z for z in data if z["mag"]<25.]
    faint=[z for z in data if z["mag"]>=25.]
    a=sum(z["red"] for z in bright);b=sum(z["red"] for z in faint)
    c=len(bright)-a;d=len(faint)-b
    return {
      "bright_m444_lt25":{"red":a,"nonred":c,"density_ratio":rate_ratio(a,c,red_area,total_area)},
      "faint_m444_ge25":{"red":b,"nonred":d,"density_ratio":rate_ratio(b,d,red_area,total_area)},
      "red_bright_vs_faint_odds_ratio_relative_to_nonred":a*d/(b*c) if b*c else None,
      "fisher_two_sided_exploratory_p":fisher_two_sided([[a,b],[c,d]])}

def analyze(rows,red_area=RED_AREA,total_area=TOTAL_AREA):
    data=[];ids=set()
    for row in rows:
        key=int(row["candidate_id"])
        if key in ids:raise ValueError("duplicate ID")
        ids.add(key)
        mag=float(row["mag_F444W"]);red=int(row["red365"])
        if red not in (0,1) or not math.isfinite(mag):raise ValueError("invalid source")
        data.append({"mag":mag,"red":red,"ra":float(row["ra_deg"]),"dec":float(row["dec_deg"])})
    bins=[]
    for lo,hi in zip(EDGES[:-1],EDGES[1:]):
        sub=[z for z in data if lo<=z["mag"]<hi]
        nr=sum(z["red"] for z in sub);no=len(sub)-nr
        bins.append({"lower_mag_inclusive":None if lo==-math.inf else lo,
                     "upper_mag_exclusive":None if hi==math.inf else hi,
                     "red365_count":nr,"nonred_count":no,
                     "density_ratio":rate_ratio(nr,no,red_area,total_area)})
    groups=proximity_groups(data)
    return {
      "status":"EXPLORATORY_MAGNITUDE_STRATIFIED_COUNTS_NOT_COMPLETENESS",
      "source":"frozen 125-row Zhen-Kai v13 candidate photometry fixture",
      "areas_deg2":{"red365_WHT_positive":red_area,"nonred_WHT_positive":total_area-red_area},
      "n_rows":len(data),"n_unique_ids":len(ids),"mag_bins":bins,
      **stratification(data,red_area,total_area),
      "proximity_groups_0p5arcsec":{"n_groups":len(groups),**stratification(groups,red_area,total_area)},
      "limitations":[
        "25 AB is a descriptive exploratory threshold, not preregistered; Fisher p is not confirmatory.",
        "Poisson Wald CIs assume independent counts and fixed area, neglect spatial correlation and depth.",
        "Close-pair groups are proximity components, not confirmed distinct astrophysical sources.",
        "F444W magnitude stratification does not control F277W local noise, flags, crowding or F356W SNR.",
        "Public 60mas WHT grid loses about 0.96% per-tile native area; geometric, not selection matched.",
        "No claim of intrinsic underdensity or PM recovery completeness."]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--fixture",default="pm_footprint/results/v13_epoch_candidate_photometry_fixture_20261031.csv")
    p.add_argument("--out",default="pm_footprint/results/magnitude_stratified_density_20261101.json")
    a=p.parse_args()
    with open(a.fixture,newline="") as f:r=analyze(csv.DictReader(f))
    if r["n_rows"]!=125 or r["proximity_groups_0p5arcsec"]["n_groups"]!=115:
        raise ValueError("frozen fixture unexpected counts")
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    Path(a.out).write_text(json.dumps(r,indent=2)+"\n")
    print(json.dumps({"bright":r["bright_m444_lt25"],"faint":r["faint_m444_ge25"]},indent=2))

if __name__=="__main__":main()
