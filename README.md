This was originally a fork of Marco's original NuMIFlux [repo](https://github.com/marcodeltutto/NuMIFlux) for flugg files, combined with Krishan's [modifications](https://github.com/kvjmistry/NuMIFlux) for dk2nu files. However, it turns out that as I progressively keep making severe changes (including overhauling the code structure), this is worth detaching and maintaining separately. All credit to them and other predecessors for the main meat of the code.

- Various simplifications mainly
- Plan to improve some more areas to ensure easy comparisons
- Currently, `FluggFlux.cc` looks at FLUGG ntuples and `Dk2NuFlux.cc` looks at Dk2Nu ntuples
- Check `<flux_type>_example.py` in `examples/` directory for example usage

## Instructions

After setting up `uboonecode` or optionally only the dependencies (see below)
```
source setup_numiana.sh
make all
```
If you just want to work with individual `<flux_type>`s where `<flux_type>` is either `dk2nu` or `flugg`, one can do
```
make base
make <flux_type>
```
To start afresh, as usual one can do
```
make clean
```

### Dependencies
- `flugg` just needs `root` basically
- `dk2nu` pulls in additionally `ppfx`, `boost`, `dk2nu`

### File Locations
- `flugg` : originally available on bluearc
    - FHC : `/nusoft/data/flux/blackbird-numix/flugg_mn000z200i_rp11_lowth_pnut_f11f093bbird_target/root`
    - RHC : `/nusoft/data/flux/blackbird-numix/flugg_mn000z-200i_rp11_lowth_pnut_f11f093bbird_target/root`
    - Also copied to
        - `/exp/uboone/data/users/bnayak/ppfx/flugg_studies/flugg_files`
        - `/pnfs/uboone/persistent/users/bnayak/flux_files/flugg_files` (dcache)
    - See [wiki](https://cdcvs.fnal.gov/redmine/projects/numi-beam-sim/wiki/Locations_of_shared_files) for more information

- `dk2nu` :
    - old (bugged) : `/cvmfs/uboone.osgstorage.org/stash/uboonebeam/numi_dk2nu_zero_threshold`
    - new (various bugfixed) : `/pnfs/uboone/persistent/users/bnayak/flux_files`

#### Picking a dk2nu production under `/pnfs/.../flux_files`

The `/pnfs/uboone/persistent/users/bnayak/flux_files` directory fans out into the
different dk2nu campaigns that have been staged for MicroBooNE.  The main
subdirectories you will see are:

| Directory | What it contains | When to use it |
| --- | --- | --- |
| `flugg_files/` | Legacy FLUGG ntuples mirrored from BlueArc. | Only if you explicitly need the FLUGG production described above. |
| `uboone_beamsim_g4.10.3.p01c/` and `uboone_beamsim_g4.10.3.p03e/` | Early Geant4 10.3 dk2nu re-productions.  They predate the geometry and focusing bug fixes. | Historical comparisons; **not** recommended for new analyses. |
| `uboone_beamsim_g4.10.4/` | The standard MicroBooNE bug-fixed dk2nu sample produced with Geant4 10.4 and the corrected NuMI beam line. | Use this for most current work—it is the dataset the repository examples assume. |
| `uboone_beamsim_g4.10.4_nothresh/` | Same as above but with the neutrino energy filter removed (“zero threshold”).  Files are larger and were mainly for validation. | Only if you need the no-threshold production; otherwise stick with the standard 10.4 set. |
| `uboone_geometrybugfix/` | Transitional samples used while validating the geometry fix. | Ignore for normal workflows. |
| `end/`, `nova/` | Non-MicroBooNE productions staged in the same area. | Do not use for MicroBooNE flux work. |

If you are unsure, default to `uboone_beamsim_g4.10.4/` for MicroBooNE dk2nu
studies and to `flugg_files/` only when reproducing the older FLUGG-based
analyses.

#### Anatomy of `uboone_beamsim_g4.10.4/`

Once you descend into `uboone_beamsim_g4.10.4/` you will encounter three
top-level folders:

| Subdirectory | What it holds |
| --- | --- |
| `me000z200i/` | Forward horn current (FHC) dk2nu production.  The "200i" tag is the horn current sign/convention used by NuMI jobs. |
| `me000z-200i/` | Reverse horn current (RHC) dk2nu production. |
| `from_isafa/` | A convenience mirror produced by Isafá that bundles a handful of merged FHC chunks (`g4numi_10.4_FHC_*.root`).  You can read them directly, but most workflows instead use the canonical `me000z±200i` trees. |

Each horn-current folder is organised by production run: `run0/`, `run1/`, …
Inside a run directory you will find three subfolders:

* `files/` – the actual dk2nu ntuples you want to process.  For example,
  ```
  /pnfs/uboone/persistent/users/bnayak/flux_files/uboone_beamsim_g4.10.4/me000z200i/run0/files/g4numi_minervame_me000z200i_0*.root
  ```
  chains every FHC chunk from `run0`.  Point `Dk2NuFlux` at that glob (or write
  it to a file list) and the reader will iterate through the entire sample.
* `CACHE/` – the job-submission payload: `g4numi.mac`, the wrapped job script,
  and the archived UPS products.  These are useful only if you need to inspect or
  reproduce the original beamline simulation.
* `log/` – per-job stdout/stderr.  When debugging failed grid jobs or sanity
  checking the production history, inspect the corresponding log file (they are
  named `g4numi_me000z±200i_<job>_<try>.log`).

The RHC directory mirrors the same structure with
`g4numi_minervame_me000z-200i_*` files.  Mixing horn polarities is as simple as
concatenating both globs in your file list before instantiating `Dk2NuFlux`.

#### Loading both horn polarities at once

1. Build a small text file that lists the FHC and RHC dk2nu globs you want to
   process.  For the Geant4 10.4 production the relevant patterns are:

   ```bash
   cat <<'EOF' > dk2nu_g4104_fhc_rhc.files
   /pnfs/uboone/persistent/users/bnayak/flux_files/uboone_beamsim_g4.10.4/me000z200i/run*/files/g4numi_minervame_me000z200i_*.root
   /pnfs/uboone/persistent/users/bnayak/flux_files/uboone_beamsim_g4.10.4/me000z-200i/run*/files/g4numi_minervame_me000z-200i_*.root
   EOF
   ```

   `TChain::Add` understands wildcards, so you can keep the globs as-is.  If you
   prefer to freeze the exact file list (for reproducibility or grid splitting),
   expand the glob with `ifdh` or `samweb` and append the concrete file names to
   the same text file instead.

2. Feed that list into the `Dk2NuFlux` constructor that accepts an explicit file
   list:

   ```python
   import ROOT

   ROOT.gSystem.Load("$NUMIANA_DIR/lib/Dk2NuFlux_cc.so")
   flux = ROOT.Dk2NuFlux(True, "dk2nu_g4104_fhc_rhc.files", "dk2nu_g4104_fhc_rhc.root")
   flux.CalculateFlux()
   ```

   The reader will chain both horn polarities, track the combined POT exposure,
   and write a single output file.  You can extend the list with additional
   production runs or filtered subsets by adding more lines to the file.

### Choosing between FLUGG and dk2nu samples

- Use the **FLUGG** samples when you want the legacy FLUKA+Geant4 flux ntuples that
  were historically used in MicroBooNE cross-checks.  They expose the FLUGG `h10`
  tree, so they are quick to read but do **not** embed dk2nu-style decay metadata
  or per-decay PPFX reweighting hooks.  Analyses that just need the pre-weighted
  flux spectra out of FLUGG can stay within this pathway.
- Use the **dk2nu** samples when you need the modern Geant4 flux ntuples with the
  `dk2nuTree`/`dkmetaTree` pair.  Those ntuples carry the per-decay information
  required by PPFX and the common NuMI reweighting infrastructure, so they are the
  right choice for any production analysis or when you intend to run universes.
- Both readers in this repository produce the same flattened output interface, so
  you can switch between them with the same downstream code.  The main trade-off is
  whether you need PPFX-capable metadata (dk2nu) or just want the faster FLUGG
  spectra.

### Working with the FLUGG file lists

The FLUGG directories on dCache/BlueArc contain one ROOT file per FLUGG job.  The
file stem encodes the beam configuration (e.g. horn polarity, thin-target setting),
and the trailing `_target_####.root` number simply increments with each job.  To use
the full data set you usually glob the files, for example

```
/pnfs/uboone/persistent/users/bnayak/flux_files/flugg_files/rhc/flugg_mn000z-200i_rp11_bs1.1_pnut_lowth_f11f093bbird_target_7*.root
```

Pass that glob to `FluggFlux` (or place the expanded file list into a text file and
use the `isfilelist` constructor) and ROOT will chain every chunk automatically.
When submitting to the grid, you can also split the list into subsets if you only
need a limited POT exposure.

### Job Submission
- Usage for `jobsub/submit_grid.sh`. Arguments
    - `-p` : tar and ship local ppfx installation (`$PPFX_DIR` must be set and pointing to local ppfx installation)
    - `-s` : set seed for PPFX universes, pulls PPFX config from `dk2nu/ppfx/inputs_ubnumi_multisim.xml`
        - multisim xml file runs 100 universes currently, modify to run your own set
    - `-n` : number of jobs to run (max 999 jobs for grid script reasons in `jobsub/numiana_job.sh`)
    - `-f` : number of files per job to run (default 1)
    - `-i` : folder containing input `dk2nu` files. The submission takes the first `n*f` files and runs `f` file(s) per job
    - `-o` : job output directory
- Runs `jobsub/numiana_job.sh` on the grid which in turn runs the python macro `jobsub/run_ppfxunivs.py`

### Changes

- ~Common output interface, currently its replicated across both `<flux_type>` classes~
- Parallel Processing
    - ~Atleast for FLUGG, this should be easy~
        - Check `examples/flugg_parallel.sh -h`
    - ~Can do for Dk2Nu as well but PPFX has a bunch of static instances that is likely not thread safe~
        - just submit jobs lol..
    - ~Parallelism for analyzing output TTrees from this code~
        - Check `examples/analyzer_example.C` for analyzing the outputs
        - This uses some newer experimental ROOT features like `RDataFrame` and `RVecOps`.
        - I tested this on `root 6.26`, your mileage may vary based on your version
        - With a good SSD and 15 cores, I'm now able to produce 100 flux universes from 50M neutrinos in about 90s
