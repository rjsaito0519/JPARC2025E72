#!/bin/sh

VAR=`date`
git add .
git commit -m "$VAR"
git push origin main