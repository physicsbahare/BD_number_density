"""Package completed depth-limited runs and regenerate compact paper tables."""
import json,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.cosmosweb_counts import TEFF_DENSITY_2024


def run():
    repo=Path(__file__).resolve().parents[1];root=Path('/home/bahareh/results')
    out=repo/'data/ali_depth_limited_2026-11-28';out.mkdir(parents=True,exist_ok=True)
    base=json.loads((root/'ali_final_counts/count_summary.json').read_text())
    bins=pd.read_csv(root/'ali_final_counts/temperature_component_counts.csv')
    bins.to_csv(out/'central_temperature_component_counts.csv',index=False)
    shutil.copyfile(root/'ali_final_counts/count_summary.json',out/'reused_central_execution.json')
    parts=[]
    for name in ['sensitivity','convergence','native_correlation']:
        folder=root/('ali_depth_limited_'+name)
        parts.append(pd.read_csv(folder/'sensitivity_bins.csv'))
        shutil.copyfile(folder/'execution.json',out/(name+'_execution.json'))
        shutil.copyfile(folder/'distance_boundaries.csv',out/(name+'_distance_boundaries.csv'))
    data=pd.concat(parts,ignore_index=True).drop_duplicates(['scenario','teff_lo_K','teff_hi_K'],keep='last')
    data.to_csv(out/'sensitivity_temperature_bins.csv',index=False)
    rows=[]
    for scenario,d in data.groupby('scenario',sort=False):
        row=dict(scenario=scenario,total=float(d.total.sum()),thin=float(d.thin.sum()),thick=float(d.thick.sum()),halo=float(d.halo.sum()),
                 measured_density_sigma=float(np.sqrt(np.nansum(d.measured_density_sigma**2))),
                 lower_limit_bin_count=float(d.loc[d.density_status=='lower_limit','total'].sum()),quadrature_error=float(d.quadrature_error.sum()))
        rows.append(row)
    pd.DataFrame(rows).to_csv(out/'sensitivity_summary.csv',index=False)
    cdf=json.loads((root/'ali_depth_limited_cdf_uncertainty.json').read_text())
    shutil.copyfile(root/'ali_depth_limited_cdf_uncertainty.json',out/'cdf_sampling_uncertainty.json')
    main=base['results']['five_empirical'];four=base['results']['four_matched'];area=base['five_area_deg2']
    selected=bins[bins.selection=='five_empirical'].copy()
    selected['effective_selected_volume_pc3']=selected.apply(lambda r:r.total/TEFF_DENSITY_2024[(int(r.teff_lo_K),int(r.teff_hi_K))][0],axis=1)
    selected.to_csv(out/'five_band_bins_and_effective_volumes.csv',index=False)
    n=main['total'];old=23.843714
    result=dict(status='ACCEPTED_EXPLICITLY_DEPTH_LIMITED_POPULATION_PREDICTION',
        base_commit='d431c5efec8b740e47323480c744d96a654a6aff',
        estimand='Expected signed native-SNR photometric passes within intrinsic F444W aperture signal >= 5 empirical background sigma. Not an exact catalog-selection prediction.',
        final_count=n,area_deg2=area,surface_density_arcmin2=n/(3600*area),
        thin=main['thin'],thick=main['thick'],halo=main['halo'],
        measured_density_sigma=main['measured_density_uncertainty'],
        conditional_cdf_sampling_sigma=cdf['conditional_cdf_sampling_standard_error'],
        step_cdf_diagnostic_count=cdf['expected_count_ecdf_step_quadrature'],
        interpolated_minus_step_cdf_count=n-cdf['expected_count_ecdf_step_quadrature'],
        partial_density_and_cdf_sigma=float(np.hypot(main['measured_density_uncertainty'],cdf['conditional_cdf_sampling_standard_error'])),
        lower_limit_bin_contribution=main['incomplete_density_bin_contribution'],
        four_original_count=base['results']['four_original']['total'],four_original_area_deg2=base['four_area_deg2'],
        four_matched_count=four['total'],f814w_reduction=four['total']-n,f814w_reduction_fraction=1-n/four['total'],
        previous_conditional_count=old,change_from_previous=n-old,fractional_change_from_previous=n/old-1,
        spatial_512_to_1024_delta=float(data.loc[data.scenario=='convergence','total'].sum()-n),
        no_fixed_distance_cap=True,all_distance_result_still_divergent=True,
        exact_auto_positivity_joint_distribution_modelled=False,
        uncertainty_limitations='Partial measured-density and finite-CDF uncertainty only; lower-limit bins, boundary-definition sensitivity, estimator calibration, correlations, PSF and population-model systematics are not a combined confidence interval.',
        central_result_provenance='/home/bahareh/results/ali_final_counts/count_summary.json; 512 cells/tile, depth domain, catalog SNR, source-noise nodes 15; old rmax_pc metadata is unused in depth mode.')
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    plt.rcParams.update({'svg.fonttype':'none','svg.hashsalt':'ali-depth-limited'})
    fig,ax=plt.subplots(figsize=(6.4,4))
    vals=[next(r['total'] for r in rows if r['scenario']=='boundary4'),n,next(r['total'] for r in rows if r['scenario']=='boundary6')]
    ax.plot([4,5,6],vals,'o-');ax.set(xlabel='Intrinsic F444W depth threshold / empirical sigma',ylabel='Expected five-band count',title='Population-boundary sensitivity (not confidence bounds)');ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'boundary_sensitivity.svg',metadata={'Date':None});plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));x=np.arange(len(selected));bottom=np.zeros(len(selected))
    for name in ['thin','thick','halo']:
        ax.bar(x,selected[name],bottom=bottom,label=name);bottom+=selected[name].to_numpy()
    ax.set_xticks(x,[f'{int(r.teff_lo_K)}–{int(r.teff_hi_K)}' for _,r in selected.iterrows()],rotation=55,ha='right',fontsize=8)
    ax.set(ylabel='Depth-limited expected count',xlabel='Effective temperature (K)');ax.legend();fig.tight_layout()
    fig.savefig(out/'temperature_components.svg',metadata={'Date':None});plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
