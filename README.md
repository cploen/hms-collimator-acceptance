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
