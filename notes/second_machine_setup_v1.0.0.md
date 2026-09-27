# Setting up a second machine with tudatpy 1.0.0 (conda), identical to the first

Goal: a second machine that gives the same results as the first one (tudatpy 1.0.0 from the
`tudat-team` conda channel, Python 3.12).

## 1. Best option: clone the first machine's environment exactly

On the **first** machine (the one that works), with the tudat environment active:

```bash
conda list --explicit > tudat100_explicit.txt      # exact package URLs + builds (same OS only)
conda env export --no-builds > tudat100.yml         # all versions, no build strings
conda env export --from-history > tudat100_hist.yml # only the packages you asked for (most portable)
conda list tudatpy                                  # note the build string, e.g. py312h3558349_3
```

On the **second** machine:

- **Same OS as the first** (e.g. both Windows): `conda create -n tudat100 --file tudat100_explicit.txt`.
  This gives bit-identical packages.
- **Different OS:** `conda env create -n tudat100 -f tudat100.yml`. This gives the same versions but a
  different build of tudatpy. See §4 for what differs numerically.
  - A full `env export` from Windows lists Windows-only packages (`ucrt`, `vc`, `vs2015_runtime`, ...)
    and one from Linux lists Linux-only ones (`libgcc-ng`, `_openmp_mutex`, ...). If the solver
    reports them as not found, delete those lines, or use `tudat100_hist.yml` and pin versions from
    the table in §2.

## 2. From scratch (if the first machine's export is not available)

```bash
conda create -n tudat100 -c tudat-team -c conda-forge python=3.12 tudatpy=1.0.0 numpy=1.26 scipy matplotlib pandas spiceypy pyerfa mpmath
conda activate tudat100
```

Versions in the working environment used for this work (Linux):

| package | version |
|---|---|
| tudatpy | 1.0.0 (tudat-team, py312h3558349_3) |
| tudat-resources | 2.4 |
| python | 3.12 |
| numpy | 1.26.4 |
| scipy | 1.17.1 |
| spiceypy | 8.2.0 |
| pyerfa | 2.0.1.5 |
| cspice | 67 |
| boost-cpp | 1.84.0 |

Notes:
- Pin `tudatpy=1.0.0`. Newer tudatpy versions renamed several APIs, e.g. `ancilliary` → `ancillary`.
- The solver can be slow; `mamba` or `micromamba` is much faster than `conda` here.
- tudat-resources downloads SPICE kernels and other data on first use and needs internet access once.

## 3. Check the install

```bash
python -c "import tudatpy; print(tudatpy.__version__)"          # 1.0.0
python -c "from tudatpy.interface import spice; spice.load_standard_kernels(); print('spice ok')"
conda list tudatpy                                               # compare the build string with machine 1
```

Then run the same small script on both machines and compare the outputs, for example the
standard averaged-Doppler test (Cassini 1-s harness, `2026-09-24_avg_doppler_floor.py <tag> --csv`
with its `data/cassini/` inputs). On v1.0.0, metric A (averaged − instantaneous at 1 s) should give a
robust scatter of about 57 mHz, and the n-way range point noise about 0.4 mm.

## 4. What will differ if the OS differs (Windows vs Linux)

- **Averaged Doppler (`dsn_n_way_averaged_doppler`):** the per-point numerical noise is about
  1 mHz at Ka for 60-s counts. Its exact values differ between platforms, because MSVC's
  `long double` is a plain double. Measured Linux vs Windows on the Rhea R1 package: means agree to
  ≤ 0.05 mHz, point differences 0.4–0.6 mHz white. Full OD estimates agree within ~0.2σ.
  Compare prefits, not single-run estimates.
- **Everything without that differencing** (instantaneous Doppler, range, propagation) agrees
  closely across platforms.
- The fix for the averaged-Doppler noise ([atogni/tudatpy#1](https://github.com/atogni/tudatpy/pull/1),
  against v1.0.0, still a draft) is **not** in the conda package. Using it means building v1.0.0 from
  source with the patch (branch `claude/dsn-averaged-doppler-fix-v1.0.0`).

## 5. Data and kernels

Copy these from the first machine; they are not part of the conda package:
- Cassini SPK `060111R_SCPSE_05320_05348.bsp`, `de442s.bsp`, `sat441.bsp`, `lsk_981207.tls`,
  `cas_v44.tf`, and the Cassini SCLK/CK kernels used by the Rhea driver;
- the Rhea package (`rhea_r1_pc60_package`) and the `data/cassini/export/` CSVs for the 1-s harness.

Keep the same directory layout, or adjust the paths at the top of the scripts
(`CASSINI_ROOT`, `SAT_SPK`, `PLANET_SPK` for the Rhea driver).
