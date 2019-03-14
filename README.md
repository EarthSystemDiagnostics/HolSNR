**HolSNR** These R scripts (and processed data) belong to the manuscript:
Reschke, M., Rehfeld, K., Laepple, T.: Empirical estimate of the signal content of Holocene temperature proxy records, Climate of the Past, 2019
and allow reproducing the main analyses presented in the paper.

*EstimatingCorOverDistance.R* allows estimating the spatial correlation of proxy records and model data (MPI6k and T21k) and *EstimatingUnSNR.R* allows estimating the time uncertainty and signal-to-noise ratios of the proxy records by estimating the mismatch between the correlation of neighboured modified model data and proxy records. Both R scripts are only an example of the analyses applied to one dataset (here: M13) and necessitate some additional functions summarised in *additionalFunctions.R*.
The datasets M13, LH14 and R18 contain the following components (here M13 as an example):
- a list *M13_meta* containing information on the metadata (1) name of the record, (2) latitude, (3) longitude, (4) mean resolution, (5) proxy type
- a list *M13_Prx* of proxy records restricted to 6 to 0ky BP and saved as zoo-objects (requires the R-package *zoo*)
- a list *M13_MPI6k* of model data extracted from MPI6k at proxy position, subsampled to proxy resolution and saved as zoo-objects (used for estimating correlations over distance)
- a list *M13_T21k* of model data extracted from T21k at proxy position, subsampled to proxy resolution and saved as zoo-objects (used for estimating correlations over distance)
- a list *M13_rawMPI6k* of annual model data extracted from MPI6k at proxy position (used for estimating time uncertainty and signal-to-noise ratios)
- a list *M13_rawT21k* of annual model data extracted from T21k at proxy position and restricted to 6 to 0ky BP (used for estimating time uncertainty and signal-to-noise ratios)
For the availability of the original datasets, please see the section *Data availability* in the paper.

Please contact Maria Reschke (mreschke@awi.de) at the Alfred-Wegener-Institute, Helmholtz Centre for Polar and Marine Research, Germany, for further information.
