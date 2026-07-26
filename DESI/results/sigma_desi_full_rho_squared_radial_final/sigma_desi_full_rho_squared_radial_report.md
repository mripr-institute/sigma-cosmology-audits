# DESI Full rho r^2 Audit

\mu=0.082912607552,
\qquad
\gamma=0.38603416,

r_\ast^2=9.470447610693826,

r^2=\|Z\hat n-C\|^2,

N_{\rm retained}=13,097,304,

\mathrm{pooled\ retained\ r^2\ coverage}
=
0.9708459316265691.

## Fixed Identity

\sigma(r^2)=\log\!\left(\frac{\mu^2(1+\gamma r^2)}{\mu+\gamma}\right)-\mu r^2

r_\ast^2=\frac1\mu-\frac1\gamma

r^2=\|Z\hat n-C\|^2

- mu = 0.082912607552
- gamma = 0.38603416
- r_star_squared = 9.470447610693826
- r_tail_squared = 21.531339548905194

## Corrected r^2 Geometry

- completed closure manifest SHA-256 = 6479c51b15ed642d1870c4e48f3d56590f6f4651068ebf6d20295ff2ae8ba74a
- best_scan_offset = -2.8000000000000007
- scan_coordinate = [0.19517549336413728, 1.8538582506826466, -2.0892860295244273]

## Corpus Readout

- retained objects = 13097304
- pooled retained r^2 coverage = 0.9708459316265691
- pooled retained r^2 uncovered = 0.02915406837343093

## Retained r^2 Coverage

- BGS: 0.048767788150996694 to 50.8793377097994; support mass = 0.9334414062448388
- LRG: 0.006512160949640133 to 83.55451052247751; support mass = 0.9933340908424543
- ELG: 0.0041593731511966325 to 83.81182491734174; support mass = 0.9934904763565024
- QSO: 0.0010850585494512899 to 70.70884690283145; support mass = 0.9834164971018396

## Fast/Direct Route Equivalence

- r^2 execution route = compare
- max abs = 2.842170943040401e-14
- pass = True
- full retained corpus compared = True

## Source Identity Checks

- r_star_squared_formula_delta = 0.00000000000000051595955919645144576470473423687570734575656339789867117351794511858186583843836
- d_sigma_d_r_squared_at_r_star_squared = 0.0000000000000000035469642430692814659591215443660375053324280564380496150977876937646978590632283
- d2_sigma_d_r_squared2_plus_mu_squared = -0.00000000000000000058817610857316015284834712537741966255793298612543541682930062552039566706519814
- rho_r_squared_equals_exp_sigma_max_abs = 1.1985855186283733264973703418407940551402602312002927470616517106228423184595114e-4098
- F_r_squared_zero = 0.0
- F_r_squared_far = 1.0
- fast_vs_direct_r_squared_max_abs = 2.842170943040401e-14
- direct_coordinate_residual_max_abs = 2.842170943040401e-14
- linear_radius_path_pass = True
- forbidden_calls = []

## Input and Output Custody

- source SHA-256 = f1e5ce2f915a7cccf42c807a35f4bbba8e42fe6deaa041ce30c2282a32ed130e
- cache_schema_version = desi_full_rho_squared_radial_cache_v2

- BGS data SHA-256 = 6f53fa8833dc207e36070857c857a1bf7747bd4fe29f8c58f2f563d9024aa179
- BGS random SHA-256 = 52b3e8c2528824349214d49034cb3c752b6703876027ac1c5c781cc78ca6f15a
- LRG data SHA-256 = 4b8ad646de35c7ac17d7fc19dc2432914dd61ceef016a3eed9aabb19af0a4dff
- LRG random SHA-256 = 12065a89957189993b3592ed80ba243f19bb2375fd47ed17197fa6b559871376
- ELG data SHA-256 = 80f3365934a91d31c54af2741cc21953d5f52f265d2296a9666b1094bf7f198f
- ELG random SHA-256 = 12065a89957189993b3592ed80ba243f19bb2375fd47ed17197fa6b559871376
- QSO data SHA-256 = c0b1ffe67794c14664fb5a8f205cf44b2631ee7bfcce9c2c318ca473846c6491
- QSO random SHA-256 = 12065a89957189993b3592ed80ba243f19bb2375fd47ed17197fa6b559871376

## Detailed Per-Tracer Readouts

- BGS: retained_count = 4549176; CDF_L1 = 0.07431060171510609; KS_max = 0.2615243980570453; r^2 density RMS = 0.024500863397479585
- BGS: d log density / d(r^2) L1 = 0.3705214523126395; RMS = 0.545127698847981; r_star_squared residual = -0.0059003546982921585
- LRG: retained_count = 2831821; CDF_L1 = 0.22154339499185363; KS_max = 0.3718038704190798; r^2 density RMS = 0.020980126512983283
- LRG: d log density / d(r^2) L1 = 0.39677744594926545; RMS = 0.6398966790286399; r_star_squared residual = -0.008180011997022701
- ELG: retained_count = 3436588; CDF_L1 = 0.22590037788658857; KS_max = 0.39676679591494435; r^2 density RMS = 0.01634713412448543
- ELG: d log density / d(r^2) L1 = 0.41413382307429997; RMS = 0.6473363348734755; r_star_squared residual = -0.00886756113910094
- QSO: retained_count = 2279719; CDF_L1 = 0.2174125245111732; KS_max = 0.3850937667521875; r^2 density RMS = 0.025845044703233098
- QSO: d log density / d(r^2) L1 = 0.16892135087929488; RMS = 0.35931296230143006; r_star_squared residual = -0.008332423129378354
