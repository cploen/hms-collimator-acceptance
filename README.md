# HMS collimator acceptance study

Determine whether reconstructed production tracks intersect the physical HMS
collimator aperture. Start with good electrons; inspect projected positions and
all six pairs of xptar, yptar, ytar, delta. Record outside events and signed edge
margins. Do not impose angular acceptance cuts before diagnosing the aperture.

## Established inputs

- Example: NPS pass2 run 2189, tree T, 396680 entries.
- Scalar Double_t leaves: H.gtr.{x,y,th,ph,dp,ok},
  H.extcor.{x,y,th,ph,dp,ok,xsieve,ysieve}, H.react.{x,y,z,ok}.
- Event identifiers: fEvtHdr.fRun and fEvtHdr.fEvtNum (UInt_t).
- Previous optics electron selection: Cherenkov >2, calorimeter E/p >0.65.
  The previous optics delta and foil cuts are not baseline cuts for this study.
- Saved HCANA source at commit prefix fcef8e23e162, THcExtTarCor.cxx,
  sets target x,y in cm and projects HMS to the sieve as
  xsieve = xtar + 168*xptar; ysieve = ytar + 168*yptar (cm).
  Verify this against the production replay; the source is a reference, not
  confirmed provenance for NPS pass2.

## First-pass results: run 2189

Reported from interactive ROOT checks on ifarm:

- Electron branches H.cer.npeSum and H.cal.etottracknorm confirmed present.
- Baseline: H.gtr.ok, Cherenkov >2, track-normalized calorimeter E/p >0.65,
  at least one NPS cluster with time in [147,153] ns. Count each event once.
- Selected events: 132627. No additional cluster energy cut was applied.
- Both gtr and extcor target-to-sieve residual plots displayed zero mean and
  standard deviation for all selected events. This checks internal projection
  consistency, not independent optics accuracy.
- At provisional z=166.4 cm, nominal large-collimator front aperture:
  7229/132627 events outside (5.45%). Geometry and alignment remain unconfirmed.
- Nominal centered front aperture in transport coordinates (cm):
  abs(x)<=11.646, abs(y)<=4.575,
  abs(x)/11.646+abs(y)/4.575<=1.5.

## Agreed implementation scope

Keep histograms and canvases in ROOT output files. Produce entrance and exit
projections, physical aperture overlays, all six pairs of xptar/yptar/ytar/delta,
inside/outside comparisons, signed boundary margins, vertex dependence, and
outside-event identifiers. Save the selection and geometry with the output.
Develop on one representative run, then apply compact run-level checks to all
production runs within a setting before combining compatible runs.

Hall A R-functions supply context for correlated boundaries and margins; the
physical HMS aperture defines this study. Epsilon-dependent analysis is deferred
until Christine discusses the desired scope with Charles.

## Pending checks

- Confirm installed collimator, entrance/exit z positions, dimensions, alignment,
  and the run's replay geometry/mispointing parameters from recorded sources.
  Sieve at 168 cm must not be equated to the collimator entrance.
- Implement the reusable ROOT macro and validate it on ifarm.
- Record run-setting metadata before combining runs.

## References

- https://hallcweb.jlab.org/wiki/index.php/Main_Page#HMS_Collimator/Sieve
- https://hallaweb.jlab.org/publications/Technotes/files/2001/01-055.pdf

Outside fractions among recorded events are reconstruction consistency
diagnostics, not absolute efficiencies.

## Running the first pass

From the repository root (ROOT 6):

```sh
root -l -b -q 'macros/hms_acceptance.C("/cache/hallc/c-nps/analysis/pass2/replays/production/nps_hms_coin_2189_0_1_-1.root","results/run2189_v1")'
```

The macro reads **config/study.conf** (the authoritative runtime configuration).
The JSON selection file is an earlier baseline record. To use another geometry
or selection, copy the .conf and supply its path as the third macro argument.

Outputs: acceptance.root (histograms, canvases, configuration and counts),
outside_events.csv, summary.txt. Open the ROOT file in a TBrowser to view the
saved canvases. Results are ignored by git. Use a new output directory on rerun;
the macro refuses to overwrite an existing acceptance.root.

Outside means outside either face; front and back counts are also saved
separately. CSV failed_face_mask is 1=front, 2=back, 3=both. The event is counted
once regardless of how many NPS clusters satisfy the timing selection. Nonfinite
target coordinates are counted separately and excluded from diagnostic histograms.
Histograms have fixed display ranges; inspect underflow/overflow before interpreting
fractions. Vertex-fraction bins with no events have no acceptance information.

The provisional back plane is front+6.35 cm, using documented back dimensions.
Both faces share a center and orientation. For a straight track and a linearly
flared convex octagon, passing both faces implies passage through the interior.
These assumptions remain subject to geometry confirmation.

Local validation: interpreted ROOT 6.36 synthetic test checked timing rejection,
central passage, horizontal/vertical/diagonal failures, an exit-only failure, saved
counts, and canvas persistence. Production validation on ifarm is pending. Local
ACLiC compilation was blocked by a macOS C++ header-path issue; the documented
command uses interpreted ROOT.

## Performance

The macro enables only leaves required by the target variables, identifiers,
and configured selection (including Sum$ cluster arrays and their counters).
Selection expressions should name leaves directly, rather than use external aliases.
It prints an initial message and progress every 10000 entries. A delay before
"Starting" can occur during file opening or formula setup; a delay after a progress
message can occur during reading. Farm I/O performance has not been measured locally.
The optimized path was checked with a variable-length cluster array and multiple
in-window clusters per event; selected-event counts remained unchanged.

## KinC_x36_5 run-list context

User-supplied list: 145 runs, 1369 segments. HMS angles are -16.43 or -16.44 deg;
beam energies span 10.54089453–10.5437901 GeV. Targets include LD2 (75 runs), LH2
(49), dummy (19), and carbon (2). NPS distances include 307, 407, and 357 cm;
NPS angles and trigger metadata also vary. HMS central momentum is not supplied.
Do not pool targets or assume identical selection conditions across this list.
The NPS cluster timing requirement can couple the selected HMS population to NPS
configuration. Retain target, trigger, and NPS configuration labels in run comparisons.
The initial 2189 check uses only segment 0 of 12 listed segments for that run.
