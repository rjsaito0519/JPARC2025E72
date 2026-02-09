#!/bin/bash

# Remove old links/files if they exist
unlink set_tdc_window.py 2>/dev/null
unlink update_param.py 2>/dev/null
unlink create_runlist.py 2>/dev/null
unlink run_hdprm.py 2>/dev/null
unlink run_hdphc.py 2>/dev/null

# Link new scripts
ln -s calibration/scripts/set_tdc_window.py
ln -s calibration/scripts/update_param.py 
ln -s calibration/scripts/create_runlist.py
ln -s calibration/scripts/run_hdprm.py
ln -s calibration/scripts/run_hdphc.py

chmod +x calibration/scripts/*.py
