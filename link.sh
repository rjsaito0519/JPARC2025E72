#!/bin/bash

unlink set_tdc_window.py
ln -s scripts/set_tdc_window.py
unlink update_param.py 
ln -s scripts/update_param.py 
unlink create_runlist.py 
ln -s scripts/create_runlist.py

chmod +x *.py
chmod +x *.sh
