#!/bin/bash

scp -r html/* login.math.berkeley.edu:/web/bio/docs/eXpress
ssh login.math.berkeley.edu 'chmod -R 755 /web/bio/docs/eXpress/downloads'