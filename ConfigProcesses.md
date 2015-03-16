# Configuration Processes #

The processes configuration key is a list of process names as defined by the model.  These names are model specific and, in some models, can be included or excluded by configuration.  Below are process names and definitions for CMAQ and CAMx.  If you don't know which processes you should use, I have made a recommendation without any warranty.  The recommendation is over inclusive and you will get warnings for processes that are missing.

CMAQ:
`  processes: [INIT, XADV, YADV, ADJC, HADV, ZADV, HDIF, VDIF, EMIS, DDEP, CHEM, AERO, CLDS, PING, FCONC]`

Defintions:
  * XADV: Advection in the E-W direction (not compatible with YAMO advection)
  * YADV: Advection in the N-S direction (not compatible with YAMO advection)
  * ADV2: XADV + YADV (not compatible with YAMO advection)
  * HADV: Advection in the E-W-N-S direction (only compatible with YAMO advection)
  * ZADV: Vertical advection
  * ADV3: Sum of XADV, YADV, and ZADV (not compatible with YAMO advection)
  * MADV: Sum of HADV and ZADV (only compatible with YAMO advection)
  * ADJC: Mass adjustment (not compatible with YAMO advection)
  * TADV: Sum of XADV, YADV, ZADV, ADJC (not compatible with YAMO advection)
  * HDIF: Horizontal diffusion
  * VDIF: Vertical diffusion
  * TDIF: Sum of HDIF and VDIF
  * TRAN: Sum of XADV, YADV, ZADV, ADJC, HDIF, and VDIF (not compatible with YAMO advection)
  * TRNM: Sum of HADV, ZADV, ADJC, HDIF, and VDIF (only compatible with YAMO advection)
  * EMIS: Emissions
  * DDEP: Dry deposition
  * CHEM: Chemistry
  * AERO: Aerosols
  * CLDS: Cloud processes and aqueous chemistry
  * PING: Plume-in-grid

CAMx:
`  processes: [INIT, CHEM, EMIS, PTEMIS, PIG, A_W, A_E, A_S, A_N, A_B, A_T, DIL, D_W, D_E, D_S, D_N, D_B, D_T, DDEP, WDEP, INORGACHEM, ORGACHEM, AQACHEM, ACHEM, FCONC]`

Definitions:
  * INIT: Initial concentration
  * CHEM: Gas-phase chemistry
  * EMIS: Gridded emissions
  * PTEMIS: Point emissions
  * PIG: Plume-In-Grid
  * A\_W: Advection from West
  * A\_E: Advection from East
  * A\_S: Advection from South
  * A\_N: Advection from North
  * A\_B: Advection from bottom
  * A\_T: Advection from top
  * DIL: Dilution
  * D\_W: Diffusion from West
  * D\_E: Diffusion from East
  * D\_S: Diffusion from South
  * D\_N: Diffusion from North
  * D\_B: Diffusion from bottom
  * D\_T: Diffusion from top
  * DDEP: Dry deposition
  * WDEP: Wet deposition
  * ACHEM: in CAMx v4.4 and lower all aerosol chemistry is lumped
  * INORGACHEM: Inorganic aerosol chemistry
  * ORGACHEM: organic aerosol chemistry
  * AQACHEM: aqueous aerosol chemistry
  * FCONC: Final concentration