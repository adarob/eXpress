#!/bin/bash

scp -r html/* login.math.berkeley.edu:/web/bio/eXpress
ssh login.math.berkeley.edu 'chmod -R 755 /web/bio/eXpress/downloads'