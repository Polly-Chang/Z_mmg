# Zmumugamma studies
## Contents
There are three codes I wrote to complete the efficiency & SF measurement. \
The SF are provided for different cut based and mVAIDs with R9 bins. Try to check how to provide them and write them accordingly. \
The environment I use to run "Z_mmg_coffea.py" is stored in "environment.yml".
### "environment.yml"
Create an environment by:\
`conda env create -n <new_env_name> -f environment.yml`\
new_env_name: name the environment as you want

### "Z_mmg_coffea.py"
Run this code to select the Z->mumugamma events and store them in "root file".

Example command:\
`python Z_mmg_coffea.py 69200 6688000 try_NanoAODv12 "/data5/NanoAOD/Run3Summer22EENanoAODv12/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8" data_EFG.root DY_postEE.root`

### "SF_eff.C"
Run it to calculate efficiency and scale factor and create the root files.

### "pu_syserr.C"
Run it to include the background systematic uncertainty, pile up systematic uncertainty and statistical uncertainty then get the final summary plot with uncertainties.