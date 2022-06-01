#!/bin/bash
# setup the env
# Note - using poi="" in pyhf so use the dev version

if [[ ! -d env ]]
then
    python3 -m venv env
    source env/bin/activate
    python3 -m pip install --upgrade pip
    python3 -m pip install \
        --upgrade "git+https://github.com/scikit-hep/pyhf.git#egg=pyhf" \
        cabinetry=='0.4.1'
else
    source env/bin/activate
fi

