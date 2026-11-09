#!/bin/bash

# Remove old links/files if they exist
for name in set_tdc_window.py update_param.py create_runlist.py dst_create_runlist.py \
            run_hodo.py run_dcprm.py run_d5.py run_tpc.py \
            iterate_dc_resi.py iterate_d5_align.py iterate_hodo_chain.py calib_dispatch.py
do
  unlink "$name" 2>/dev/null
done

# Link front ends to the repository root
ln -s calibration/scripts/set_tdc_window.py
# ln -s calibration/scripts/update_param.py
ln -s calibration/scripts/create_runlist.py
ln -s calibration/scripts/dst_create_runlist.py
ln -s calibration/scripts/run_hodo.py
ln -s calibration/scripts/run_dcprm.py
ln -s calibration/scripts/run_d5.py
ln -s calibration/scripts/run_tpc.py
ln -s calibration/scripts/iterate_dc_resi.py
ln -s calibration/scripts/iterate_d5_align.py
ln -s calibration/scripts/iterate_hodo_chain.py
ln -s calibration/scripts/calib_dispatch.py

chmod +x calibration/scripts/*.py
