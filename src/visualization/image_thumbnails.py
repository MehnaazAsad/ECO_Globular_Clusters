#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 23:31:42 2017

@author: asadm2
"""

from PIL import Image
from glob import glob
from fpdf import FPDF

size = 1600,1600
files = ['hlsp_coma_hst_acs-wfc_v24_f475w_v1_ivm-drz-cl_test_acs_wfc_f475w.fits',
'hst_10861_18_acs_wfc_f475w_test_acs_wfc_f475w.fits',
'hst_10861_b4_acs_wfc_f475w_test_acs_wfc_f475w.fits',
'hst_11711_02_acs_wfc_f475w_test_acs_wfc_f475w.fits',
'hst_11711_03_acs_wfc_f475w_test_acs_wfc_f475w.fits',
'hlsp_coma_hst_acs-wfc_v24_f814w_v1_ivm-drz-cl_test_acs_wfc_f814w.fits',
'hst_10861_18_acs_wfc_f814w_test_acs_wfc_f814w.fits',
'hst_10861_b4_acs_wfc_f814w_test_acs_wfc_f814w.fits',
'hst_11711_02_acs_wfc_f814w_test_acs_wfc_f814w.fits']

for image in files:
    Image.open(image).thumbnail(size).save("thumbnail_%s_%s"%(image,"_".join(size)))

pdf = FPDF()
arr = glob('thumbnail*')
for image in arr:
    pdf.add_page()
    pdf.image(image,0,0,1600,1600)
pdf.output("ECO00026.pdf","F")
