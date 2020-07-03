#!/usr/bin/python
##
# \file         bonemeasurement.py
# \author       Bill Hill
# \date         October 2018
# \version      $Id$
# \par
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be
# useful but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
# \brief        Computes long bone endpoints from DXA scan images.
#               Image coordinates are assumed to have increasing line
#               values from head to toes and increasing column values
#               from right to left (ie person scanned lying on their
#               back and projection through table to sensor underneath).
##

from __future__ import print_function
import io
import os
import re
import csv
import sys
import json
import argparse
import subprocess
import tempfile
import traceback
import numpy as np
import ctypes as c
import math as m
import Wlz as w
import copy as cp
import pydicom
import matplotlib.image as pig

libc = c.CDLL("libc.so.6")

libc.fopen.restype = c.POINTER(w.FILE)

wlz_rgba_channel = {
  'red':        '-C 1',
  'green':      '-C 2',
  'blue':       '-C 3',
  'hue':        '-C 4',
  'saturation': '-C 5',
  'brightness': '-C 6',
  'cyan':       '-C 7',
  'magenta':    '-C 8',
  'yellow':     '-C 9',
  'modulus':    '-m'}

#======================================================================#

prog            = 'bonemeasurement.py'
version         = '1.0.0'
args            = []
assaylist       = []
def_base_dir    = os.environ['HOME'] + '/Projects/Biobank/BoneMeasurement'
def_scan_dir    = def_base_dir + '/InitialScans'
def_work_dir    = def_base_dir + '/Measurement/Working'
vrb_level       = 0                 # Level of verbosity, 0 = none, 9 = most
def_model       = '1848871'         # base file name of the model
def_format      = 'dcm'             # default image format
img_size_min    = 64                # minimum size of image
roi_obj_size    = 160               # size of roi images
min_bgd_area    = 180               # Minimum area for filling in background
min_bgd_value   = 250               # Minimum value for removing background
min_cnf_value   = 0.55              # Default minimum confidence value

#======================================================================#
# Tuning parameters for region of interest search.
#
# Control of filters when extracting ROI for registration
prm_filter_sobel  = False           # Use sobel filter
prm_filter_window = True            # Use window function
# Peak detection in image profiles
prm_peak_sigma  = 11                # Gaussian sigma for smoothed image
prm_peak_minval = 20                # Minimum peak value in smoothed image
# Centre/top of head
prm_head_s      = 0.25              # Top of head is in this fraction of the
                                    # image from top (line 0) toward bottom
prm_head_w      = 10.0              # Minimum width of head
# Tip of toes
prm_toes_s      = 0.25              # Tip of itoes is in this fraction of the
                                    # image from bottom toward top (line 0)
prm_toes_w      = 10.0              # Minimum width of toes
# Neck and shoulder search
prm_head_cen    = 0.06              # Fraction of body length from top of head
                                    # to centre of head
prm_shoulder_e  = 0.25              # Fraction of body length from start of
                                    # head in which to search for shoulders
# Centre of pelvis search
prm_addomen_s   = 0.35              # Fraction of body length from top of
                                    # head at which to start pelvis search
prm_pelvis_e    = 0.55              # Fraction of body length from top of
                                    # head at which to end pelvis search
# Hip regions
prm_hip_s       = 0.01              # Fraction of body length from
                                    # base of spine at which to start searching
                                    # for hips
prm_hip_e       = 0.15              # Fraction of body length from
                                    # base of spine at which to end searching
                                    # for hips
# Elbow region search is based on fractions of head - neck distance
prm_elbow_i     = 0.35              # Length down arm from shoulder where
                                    # arm tracking begins - discontinuities
                                    # are checked from here down
prm_elbow_s     = 0.60              # Length down arm from shoulder where
                                    # search for elbow begins
prm_elbow_e     = 0.70              # Length down arm from shoulder where
                                    # search for elbow ends
# Wrist region search is based on fractions of upper arm length
prm_wrist_s     = 0.50              # Length down arm from elbow where
                                    # search for wrist begins
prm_wrist_e     = 1.10              # Length down arm from elbow where
                                    # search for wrist ends
prm_wrist_l     = 0.06              # Fraction of found wrist - elbow
                                    # vector to add to found wrist location
# Knee region search is based on fractions of head - neck distance
prm_knee_s     = 0.60               # Length down arm from shoulder where
                                    # search for knee begins
prm_knee_e     = 1.20               # Length down arm from shoulder where
                                    # search for knee ends

#======================================================================#
# Debug stuff
#
debug_mask = 0x0000                 # Don't set this unless debugging!

debug_lut = {
    'none':       0x0000,
    'exit':       0x0001,
    'combine':    0x0002,
    'confidence': 0x0004,
    'map':        0x0008,
    'proc_assay': 0x0010,
    'proc_model': 0x0020,
    'register':   0x0030
    }

def debugValue(s): #{
  v = 0
  if(bool(debug_mask) and bool(s) and (s in debug_lut)): #{
    v = debug_mask & debug_lut[s]
  #}
  return v
#}

def debugFile(s): #{
  return prog.split('.')[0] + '_debug_' + s
#}

#======================================================================#

def isPython3(): #{
  rtn = (sys.version_info[0] == 3)
  return rtn
#}

def isVerbose(lvl): #{
  return((vrb_level > 0) and (lvl <= vrb_level))
#}

def errMsg(msg): #{
  print(prog + ' ERROR : ' + msg, file=sys.stderr)
#}

def vrbMsg(lvl, msg): #{
  if(isVerbose(lvl)): #{
    print(prog + ': ' + msg, file=sys.stderr)
  #}
#}

def runCmd(cmd): #{
  vrbMsg(5, 'runCmd() cmd = ' + cmd)
  status = 0
  try: #{
    rtn = subprocess.check_output(cmd, shell=True)
    rtn = rtn.decode('ascii')
    rtn = rtn.replace('\n','')
  except: #}{
    rtn = ''
    status = 1
  #}
  vrbMsg(5, 'runCmd() status = ' + str(status) + ', rtn = ' + rtn)
  return(status, rtn)
#}

def roiCentresFile(image): #{
  roc_file = args.work_dir + '/roc-' + str(image) + '.jsn'
  return roc_file
#}

def tmpFileBase(): #{
  tmp_file_base = args.work_dir + '/tmp'
  return tmp_file_base
#}

def hstFileBase(image): #{
  hst_file_base = args.work_dir + '/hst-' + str(image)
  return hst_file_base
#}

def bfiFileBase(image): #{
  bfi_file_base = args.work_dir + '/bfi-' + str(image)
  return bfi_file_base
#}

def smtFileBase(image): #{
  smt_file_base = args.work_dir + '/smt-' + str(image)
  return smt_file_base
#}

def regFileBase(model, assay, roc): #{
  reg_file_base = args.work_dir + '/reg-' + model +'-' + assay + '-' + roc
  return reg_file_base
#}

def roiFileBase(image): #{
  roi_file = args.work_dir + '/roi-' + str(image)
  return roi_file
#}

def trxFileBase(model, assay, roc): #{
  trx_base = (args.work_dir + '/trx-' + str(model) +'-' + str(assay) + '-' +
              roc + '-')
  return trx_base
#}

def wlzObjToNP(obj): #{
  vrbMsg(5, 'wlzObjToNP() obj = ' + str(obj))
  gtbl = {w.WLZ_GREY_INT: c.c_int,
          w.WLZ_GREY_SHORT: c.c_short,
          w.WLZ_GREY_UBYTE: c.c_ubyte,
          w.WLZ_GREY_FLOAT: c.c_float,
          w.WLZ_GREY_DOUBLE: c.c_double}
  vt = None
  ary = None
  aryc = None
  org = w.WlzIVertex2()
  err_num = c.c_int(w.WLZ_ERR_FILE_OPEN)
  bb = w.WlzBoundingBox2I(obj, c.byref(err_num))
  if(not bool(err_num)): #{
    sz = w.WlzIVertex2()
    sz.vtX = bb.xMax - bb.xMin + 1
    sz.vtY = bb.yMax - bb.yMin + 1
    org.vtX = bb.xMin
    org.vtY = bb.yMin
    gt = w.WlzGreyTypeFromObj(obj, c.byref(err_num))
    if(not (gt in gtbl)): #{
      err_num = w.WLZ_GREY_ERROR
    else: #}{
      vt = gtbl[gt]
    #}
    UPP = c.POINTER(c.POINTER(vt))
    UPV = c.POINTER(c.c_void_p)
    aryc = c.cast(0, UPV)
    err_num = w.WlzToArray2D(c.byref(aryc), obj, sz, org, 0, c.c_int(gt))
  #}
  if(not bool(err_num)): #{
    aryc = c.cast(aryc, UPP)
    ary = np.ctypeslib.as_array(aryc.contents, (sz.vtY, sz.vtX))
    w.AlcFree(aryc)
  #}
  vrbMsg(5, 'wlzObjToNP() err_num = ' + str(err_num) +
      ', org = (' + str(org.vtX) + ', ' + str(org.vtY) + ')' +
      ', ary = ' + str(ary))
  return(err_num, org, ary)
#}

def readWoolzObj(obj_file): #{
  vrbMsg(5, 'readWoolzObj() obj_file = ' + obj_file)
  err_num = c.c_int(w.WLZ_ERR_FILE_OPEN)
  obj = None
  fn = obj_file.encode('ascii') if (isPython3) else obj_file
  fp = libc.fopen(fn, 'rb')
  if(not bool(fp)): #{
    errMsg('Failed to open Woolz object file ' + obj_file + ' for reading')
  else: #}{
    err_num = c.c_int(w.WLZ_ERR_NONE)
    obj = w.WlzAssignObject(w.WlzReadObj(fp, c.byref(err_num)), None)
    libc.fclose(fp)
    if(bool(err_num)): #{
      errMsg('Failed to read Woolz object from file ' + obj_file)
    #}
  #}
  vrbMsg(5, 'readWoolzObj() err_num = ' + str(err_num) + ', obj = ' + str(obj))
  return err_num, obj
#}

def writeWoolzObj(obj_file, obj): #{
  vrbMsg(5, 'writeWoolzObj() obj_file = ' + obj_file + ', obj = ' + str(obj))
  err_num = c.c_int(w.WLZ_ERR_FILE_OPEN)
  fn = obj_file.encode('ascii') if (isPython3) else obj_file
  fp = libc.fopen(fn, 'wb')
  if(not bool(fp)): #{
    errMsg('Failed to open Woolz object file ' + obj_file + ' for writing')
  else: #}{
    err_num = w.WlzWriteObj(fp, obj)
    libc.fclose(fp)
    if(bool(err_num)): #{
      errMsg('Failed to write Woolz object to file ' + obj_file)
    #}
  #}
  vrbMsg(5, 'writeWoolzObj() err_num = ' + str(err_num))
  return err_num
#}

def readJsnFile(filename): #{
  vrbMsg(5, 'readJsnFile() filename = ' + filename)
  status = 0
  jsn = None
  try: #{
    f = io.open(filename, 'r')
    jsn = json.load(f)
    f.close()
  except: #}{
    status = 1
  #}
  return status, jsn
#}

def getPixelSz(assay, scan_dir, fmt): #{
  pix_sz = [1.0, 1.0]
  if(fmt == 'dcm'): #{
    file_path = scan_dir + '/' + assay + '.' + fmt
    img = pydicom.read_file(file_path)
    pix_sz[0] = float(img.ExposedArea[0]) / float(img.Columns)
    pix_sz[1] = float(img.ExposedArea[1]) / float(img.Rows)
  #}
  return pix_sz
#}

def meanSurroundValue(obj_tst, obj_ref): #{
  mean_value = 0.0
  err_num = c.c_int(w.WLZ_ERR_NONE)
  tmp_obj = [None, None, None]
  tmp_obj[0] = w.WlzAssignObject(
      w.WlzDilation(obj_tst, c.c_int(w.WLZ_8_CONNECTED),
          c.byref(err_num)), None)
  if(not bool(err_num)): #{
    err_num = c.c_int(w.WLZ_ERR_NONE)
    tmp_obj[1] = w.WlzAssignObject(
        w.WlzDiffDomain(tmp_obj[0], obj_tst, c.byref(err_num)), None)
  #}
  dummy = w.WlzFreeObj(tmp_obj[0])
  vrbMsg
  if(not bool(err_num)): #{
    err_num = c.c_int(w.WLZ_ERR_NONE)
    tmp_obj[2] = w.WlzAssignObject(
        w.WlzGreyTransfer(tmp_obj[1], obj_ref, 0, c.byref(err_num)), None)
  #}
  dummy = w.WlzFreeObj(tmp_obj[1])
  if(not bool(err_num)): #{
    err_num = c.c_int(w.WLZ_ERR_NONE)
    mean = c.c_double(0.0)
    area = w.WlzGreyStats(tmp_obj[2], None, None, None, None, None,
        c.byref(mean), None, c.byref(err_num))
    mean_value = mean.value
  #}
  dummy = w.WlzFreeObj(tmp_obj[2])
  return mean_value, err_num
#}

def fillBackground(tmp_file, in_idx, out_idx): #{
  vrbMsg(5, 'fillBackground() tmp_file = ' + tmp_file + 
          ', in_idx = ' + str(in_idx) +
          ', out_idx = ' + str(out_idx))
  # Read in Woolz object
  in_obj = None
  out_obj = None
  thr_obj = None
  n_rgn = None
  lbl_ary = None
  in_obj_file = tmp_file + str(in_idx) + '.wlz'
  out_obj_file = tmp_file + str(out_idx) + '.wlz'
  err_num, in_obj = readWoolzObj(in_obj_file)
  # Threshold to get regions of background, but these will also include
  # some dense bone foreground.
  if(not bool(err_num)): #{
    hilo = w.WLZ_THRESH_HIGH
    tv = w.WlzPixelV()
    tv = c.c_int(min_bgd_value)
    err_num = c.c_int(w.WLZ_ERR_NONE)
    vrbMsg(5, 'fillBackground() thresholding object read from = ' + 
        in_obj_file + ' at value ' + str(min_bgd_value))
    thr_obj = w.WlzAssignObject(
        w.WlzThresholdI(in_obj, tv, hilo, c.byref(err_num)), None)
    if(bool(err_num)): #{
      errMsg('Failed to threshold input Woolz object')
    #}
  #}
  # Label to get connected regions.
  if(not bool(err_num)): #{
    vrbMsg(5, 'fillBackground() running connected component labeling')
    n_rgn = c.c_int(0)
    max_rgn = 1000
    lbl_ary = c.cast(w.AlcCalloc(max_rgn, c.sizeof(w.WlzObject)),
        c.POINTER(c.POINTER(w.WlzObject)))
    err_num = w.WlzLabel(thr_obj, c.byref(n_rgn), c.byref(lbl_ary),
        max_rgn, 1, c.c_int(w.WLZ_8_CONNECTED))
                       
  #}
  if(not bool(err_num)): #{
    vrbMsg(5, 'fillBackground() found ' + str(n_rgn.value) +
        ' possible background regions')
    bb = w.WlzBoundingBox2I(in_obj, None)
    max_head_ln = (bb.yMax - bb.yMin + 1) / 5
    vrbMsg(5, 'fillBackground() selecting regions by size and position')
    # Only want the regions not surrounded by high valued voxels in the
    # original image. Squeeze all other regions out of the array freeing
    # them as we go.
    idx0 = 0
    for idx1 in range(0, n_rgn.value): #{
      mean, err_num = meanSurroundValue(lbl_ary[idx1], in_obj)
      vrbMsg(9, 'fillBackground() idx1 = ' + str(idx1) +
          ', mean = ' + str(mean) + ', err_num = ' + str(err_num))
      if(mean < 100.0): #{
        lbl_ary[idx0] = lbl_ary[idx1]
        idx0 = idx0 + 1
      else: #}{
        dummy = w.WlzFreeObj(lbl_ary[idx1])
        lbl_ary[idx1] = None
      #}
    #}
    n_selected = idx0
    # Form the union of the selected regions
    vrbMsg(5, 'fillBackground() forming union of ' + str(n_selected) + 
        ' selected regions')
    err_num = c.c_int(w.WLZ_ERR_NONE)
    rgn_obj = w.WlzAssignObject(
        w.WlzUnionN(n_selected, lbl_ary, 0, c.byref(err_num)), None)
    for idx0 in range(n_selected): #{
      dummy = w.WlzFreeObj(lbl_ary[idx0])
    #}
    w.AlcFree(lbl_ary)
    vrbMsg(9, 'fillBackground() union of selected regions rgn_obj = ' +
        str(rgn_obj) + ',errNum = ' + str(err_num))
    if(not bool(err_num)): #{
      vrbMsg(5, 'fillBackground() masking regions')
      # Set value 0 using mask union of selected regions
      out_obj = w.WlzAssignObject(
          w.WlzGreyMaskI(in_obj, rgn_obj, 0, c.byref(err_num)), None)
      dummy = w.WlzFreeObj(rgn_obj)
    #}
  #}
  if(not bool(err_num)): #{
    vrbMsg(5, 'fillBackground() write object to ' + out_obj_file)
    errNum = writeWoolzObj(out_obj_file, out_obj)
  #}
  if(bool(out_obj)): #{
    dummy = w.WlzFreeObj(out_obj)
  #}
  vrbMsg(5, 'fillBackground() returning status = ' + str(err_num))
  return err_num 
#}

def dist(p0, p1): #{
  vrbMsg(9, 'dist() p0 = ' + str(p0) + ' p1 = ' + str(p1))
  l = -1.0
  if((type(p0) is list) and (len(p0) >= 2) and
     (type(p1) is list) and (len(p1) >= 2)): #{
    d = [p0[0] - p1[0], p0[1] - p1[1]]
    l = m.sqrt((d[0] * d[0]) + (d[1] * d[1]))
  #}
  vrbMsg(9, 'dist() l = ' + str(l))
  return l 
#}

def findPeaks(x, sigma, min_hgt, show = False, troughs = False): #{
  import matplotlib.pyplot as plt
  vrbMsg(5, 'findPeaks() sigma = ' + str(sigma) +
         ', min_hgt = ' + str(min_hgt) +
         ', show = ' + str(show) +
         ', troughs = ' + str(troughs))
  sigma = max(sigma, 1.0)
  x = np.atleast_1d(x).astype('float64')
  if(troughs): #{
    x = 255 - x
  #}
  off = int(np.floor(m.fabs(sigma / 2.0)) + 1)
  # Subtract offset 1D values (this is DoG because alread convolved with
  # Gaussian)
  dx = x[off:] - x[:-off]
  pkp = np.where((np.hstack((dx, 0)) <= 0)  & (np.hstack((0, dx)) > 0))[0]
  # remove peaks < minimum peak height
  pkp = pkp[x[pkp] >= min_hgt]
  # Remove offset from subtraction
  pkp += int(m.floor(off / 2) - 1)
  vrbMsg(9, 'findPeaks() pkp = ' + str(pkp))
  # Find peak widths where widths are width at half height or next minimum
  # which ever is closest to peak position. Hopefully this is enough without
  # fitting Gaussians.
  pkw = np.zeros(len(pkp))
  for idp in range(0, len(pkp)): #{
    p = pkp[idp]
    y0 = x[p]
    y02 = y0 / 2
    y1 = y0
    for idx in range(p + 1, len(x), 1): #{
      y2 = x[idx]
      if((y2 < y02) or (y2 > y1)): #{
        break
      #}
      y1 = y2
    #}
    pkw[idp] = idx - p
    y1 = y0
    for idx in range(p - 1, 0, -1): #{
      y2 = x[idx]
      if((y2 < y02) or (y2 > y1)): #{
        break
      #}
      y1 = y2
    #}
    pkw[idp] += p - idx
  #}
  vrbMsg(9, 'findPeaks() pkw = ' + str(pkw))
  if(show): #{
    if(troughs): #{
      x = 255 - x
    #}
    dummy, ax = plt.subplots(1, 1, figsize=(8, 4))
    ax.plot(x, 'b', lw=1)
    ax.plot(pkp, x[pkp], '+', mec='r')
    plt.show()
  #}
  vrbMsg(5, 'findPeaks() pkp = ' + str(pkp) + ', pkw = ' + str(pkw))
  return(pkp,pkw)
#}

def matchProfilePeak(pkp, p, w): #{
  """
  Match the given peak by choosing the closest within a given range.
  Matched peak must be in range p +/- w

  Args:
  pkp:          Given array of peak positions
  p:            Position of given peak to match
  w:            Width

  Returns:      Index into the given array of the matched peak or None if
                no match found
  """
  vrbMsg(5, 'matchProfilePeak() pkp = ' + str(pkp) + 
         ', p = ' + str(p) + ',w = ' + str(w))
  pi = None
  off = 2 * w
  p0 = p - w
  p1 = p + w
  for i in range(0, len(pkp)): #{
    if ((pkp[i] > p0) and (pkp[i] < p1)): #{
      o2 = abs(pkp[i] - p)
      if(o2 < off): #{
        pi = i
        off = o2
      #}
      break
    #}
  #}
  vrbMsg(5, 'matchProfilePeak() pi = ' + str(pi))
  return pi 
#}

def matchProfilePeakPair(pkp, cpp): #{
  """
  Match a pair of peaks finding the closest pair to the given peak.

  Args:
  pkp:          Given array of peak positions
  p:            Position of given peak to match

  Returns:      An array of two indices. It's possible that the indices are
                None if there's no peak. The first peak index will always
                be lower valued and non identical to the second.
  """
  vrbMsg(5, 'matchProfilePeakPair() pkp = ' + str(pkp) + ', cpp = ' + str(cpp))
  npp = [None, None]
  for i in range(0, len(pkp)): #{
    if(npp[0] is None): #{
      npp[0] = i
    elif(npp[1] is None): #}{
      npp[1] = i
    else: #{
      vrbMsg(5, 'matchProfilePeakPair() i = ' + str(i) + ', npp = ' + str(npp))
      d0 = abs(pkp[i] - cpp)
      d1 = abs(pkp[npp[0]] - cpp)
      if(d0 < d1): #{
        j = npp[0]
        npp[0] = i
        d0 = abs(pkp[j] - cpp)
        d1 = abs(pkp[npp[1]] - cpp)
        if(d0 < d1): #{
          npp[1] = j
        #}
      #}
    #}
  #}
  if((not (npp[1] is None)) and (npp[1] < npp[0])): #{
    i = npp[0]
    npp[0] = npp[1]
    npp[1] = i
  #}
  vrbMsg(5, 'matchProfilePeakPair() npp = ' + str(npp))
  return npp 
#}

def findHighestPkPair(x, pkp): #{
  """
  Find the highest two peaks of the given peak array using the given data
  array to find heights. If there's only one peak then an array with a single
  element is returned.
  """
  vrbMsg(5, 'findHighestPkPair() x = [...], pkp = ' + str(pkp))
  mi = [0,0]
  mv = x[pkp[0]]
  for i in range(1, len(pkp)): #{
    nv = x[pkp[i]]
    if(nv > mv): #{
      mi[0] = i
      mv = nv
    #}
  #}
  mv = None
  for i in range(0, len(pkp)): #{
    if(i != mi[0]): #{
      nv = x[pkp[i]]
      if((mv is None) or (nv > mv)): #{
        mi[1] = i
        mv = nv
      #}
    #}
  #}
  if(mi[0] == mi[1]): #{
    mi = [mi[0]]
  elif(mi[0] > mi[1]): #{
    mi = [mi[1], mi[0]]
  #}
  vrbMsg(5, 'findHighestPkPair() mi = ' + str(mi))
  return mi 
#}

def findProfileHeadCT(rois, ary): #{
  # Find first substantial single peak
  vrbMsg(3, 'findProfileHeadCT()')
  sz = np.shape(ary)
  max_head_ln = int(sz[0] * prm_head_s)
  for ln in range(0, max_head_ln, 1): #{
    pkp,pkw = findPeaks(ary[ln], prm_peak_sigma, prm_peak_minval, show = False)
    if((len(pkp) == 1) and (pkw[0] > prm_head_w)): #{
      rois['head_ct'] = [pkp[0],ln]
      break
    #}
  #}
  vrbMsg(3, 'findProfileHeadCT(rois[\'head_ct\'] = ' + str(rois['head_ct']))
#}

def findProfileToes(rois, ary): #{
  # Find last substantial peak (may be either left or right)
  vrbMsg(3, 'findProfileToes()')
  sz = np.shape(ary)
  min_toes_ln = int(sz[0] - (sz[0] / prm_toes_s))
  for ln in range(sz[0] - 1, min_toes_ln, -1): #{
    pkp,pkw = findPeaks(ary[ln], prm_peak_sigma, prm_peak_minval, show = False)
    if((len(pkp) >= 1) and (pkw[0] > prm_toes_w)): #{
      rois['toes_e'] = [pkp[0],ln]
      break
    #}
  #}
  vrbMsg(3, 'findProfileHeadCT(rois[\'toes_e\'] = ' + str(rois['toes_e']))
#}


def findProfileShoulders(rois, ary): #{
  # Have found top of head and tip of toes so can estimate head size as
  # approximately 1/8th this distance, use this to help in searching for
  # the shoulders. First find neck below top of head then first triple peak
  # after start of neck
  vrbMsg(3, 'findProfileShoulders()')
  if(bool(rois['toes_e']) and bool(rois['head_ct'])): #{
    sz = np.shape(ary)
    head_s = rois['head_ct'][1]
    body_len = rois['toes_e'][1] - head_s
    min_head_centre_ln = int(head_s + (body_len * prm_head_cen))
    max_shoulders_ln = int(head_s + (body_len * prm_shoulder_e))
    # Find widest part of head
    w = [0.0, 0.0]
    for ln in range(head_s, min_head_centre_ln, 1): #{
      pkp,pkw = findPeaks(ary[ln], prm_peak_sigma, prm_peak_minval,
          show = False)
      n_p = len(pkp)
      if(n_p > 0): #{
        w[1] = (pkp[n_p - 1] - pkp[0]) + ((pkw[0] + pkw[n_p - 1]) / 2.0)
        if(w[1] > w[0]): #{
          w[0] = w[1]
        #}
      #}
    #}
    vrbMsg(5, 'findProfileShoulders() ln ' + str(ln) +
        ' head width >= ' + str(w[0]))
    neck = False
    for ln in range(min_head_centre_ln, max_shoulders_ln, 1): #{
      pkp,pkw = findPeaks(ary[ln], prm_peak_sigma, prm_peak_minval,
          show = False)
      if((len(pkp) == 1) and ((4 * pkw[0]) < (3 * w[0]))): #{
        neck = True
        break
      #}
    #}
    vrbMsg(5, 'findProfileShoulders() ln = ' + str(ln) + ' neck = ' + str(neck))
    cen = rois['head_ct'][0]
    match_cnt = 0
    # for ln in range(ln, max_shoulders_ln, 1): #{
    for ln in range(ln, 200, 1): #{
      pkp,pkw = findPeaks(ary[ln], prm_peak_sigma, prm_peak_minval,
          show = False)
      if(len(pkp) == 3): #{
        l0 = pkp[1] - pkp[0]
        l2 = pkp[2] - pkp[1]
        r = float(l0) / float(l2)
        if((r > 0.3) and (r < 1.5)): #{
          match_cnt = match_cnt + 1
          if(match_cnt > 4): #{
            # Require 3 peaks at 4 consecutive lines to avoid noise
            rois['shoulder_r'] = [pkp[0], ln]
            rois['neck_cb']    = [pkp[1], ln]
            rois['shoulder_l'] = [pkp[2], ln]
            break
          #}
        #}
      else: #}{
        match_cnt = 0
      #}
    #}
  #}
  vrbMsg(3, 'findProfileShoulders(rois[\'shoulder_l\'] = ' +
      str(rois['shoulder_l']))
  vrbMsg(3, 'findProfileShoulders(rois[\'shoulder_r\'] = ' +
      str(rois['shoulder_r']))
  vrbMsg(3, 'findProfileShoulders(rois[\'neck_cb\'] = ' +
      str(rois['neck_cb']))
#}

def findProfilePelvis(rois, ary): #{
  # Have found top of head and tip of toes. Now search for the centre
  # of the pelvis. Start at abdomen and find highest peak (of spine).
  vrbMsg(3, 'findProfilePelvis()')
  if(bool(rois['toes_e']) and bool(rois['head_ct'])): #{
    last_ln = None
    pki = None
    pkp = None
    head_s = rois['head_ct'][1]
    body_len = rois['toes_e'][1] - head_s
    abdomen_s = int(head_s + (body_len * prm_addomen_s))
    pelvis_e = int(head_s + (body_len * prm_pelvis_e))
    # Over 4 lines find average position of the highest peak (spine).
    p_mean = 0
    w_mean = 0.0
    for ln in range(abdomen_s, abdomen_s + 4, 1): #{
      x = ary[ln]
      pkp,pkw = findPeaks(x, prm_peak_sigma, prm_peak_minval, show = False)
      i0 = 0
      for i1 in range(1, len(pkp)): #{
        if(x[pkp[i1]] > x[pkp[i0]]): #{
         i0 = i1
        #}
      #}
      p_mean += pkp[i0]
    #}
    p_mean = p_mean * 0.25
    last_spp = [int(p_mean), abdomen_s + 4, 0]
    vrbMsg(5, 'findProfilePelvis() p_mean = ' + str(p_mean))
    # Continue working down until the spine peak is gone.
    for ln in range(abdomen_s + 5, pelvis_e, 1): #{
      x = ary[ln]
      pkp,pkw = findPeaks(x, prm_peak_sigma, prm_peak_minval, show = False)
      # Is there a peak at the spine position?
      mpi = matchProfilePeak(pkp, last_spp[0], prm_peak_sigma)
      vrbMsg(7, 'findProfilePelvis() ln = ' + str(ln) + ', mpi = ' + str(mpi))
      if(mpi is None): #{
        break
      else: #}{
        last_spp = [pkp[mpi], ln, pkw[mpi]]
      #}
    #}
    if((ln > abdomen_s + 5) and (ln < pelvis_e - 1)): #{
      rois['pelvis_c'] = last_spp[:2]
    #}
  #}
  vrbMsg(3, 'findProfilePelvis() rois[\'pelvis_c\'] = ' +
      str(rois['pelvis_c']))
#}

def findProfileHip(rois, ary): #{
  # Have found top of head, centre of pelvis and tip of toes. Now search for
  # the head of the femur (hip) by first finding the peak either side of the
  # centre of the pelvis, then tracking this down until it reaches a maximum
  # width.
  vrbMsg(3, 'findProfileHip()')
  if(bool(rois['toes_e']) and bool(rois['pelvis_c'])): #{
    body_len = rois['toes_e'][1] - rois['head_ct'][1]
    pelvis_c = rois['pelvis_c'][0]
    hip_s = rois['pelvis_c'][1]  
    hip_e = int(hip_s + (body_len * prm_hip_e))
    hip_s = int(hip_s + (body_len * prm_hip_s))
    tp = rois['pelvis_c'][0]
    last_ln = None
    last_pkp = [None, None]
    last_pkh = [None, None]
    for ln in range(hip_s, hip_e, 1): #{
      x = ary[ln]
      ln_from_start = ln - hip_s
      vrbMsg(5, 'findProfileHip() ln = ' + str(ln))
      pkp,pkw = findPeaks(x, prm_peak_sigma, prm_peak_minval, show = False)
      if(last_ln is None): #{
        # Find central peaks, ie either side of pelvis (using pelvis centre).
        pkpi = matchProfilePeakPair(pkp, rois['pelvis_c'][0]) 
        if((pkpi[0] is None) or (pkpi[1] is None)): #{
          last_ln = None
          break
        else: #}{
          last_pkp = [pkp[pkpi[0]], pkp[pkpi[1]]]
          last_pkh = [x[last_pkp[0]], x[last_pkp[1]]]
        #}
      else: #}{
        # Track central peaks down to find maximum peak height
        pkpi = [None, None]
        pkpi[0] = matchProfilePeak(pkp, last_pkp[0], 3 * prm_peak_sigma)
        pkpi[1] = matchProfilePeak(pkp, last_pkp[1], 3 * prm_peak_sigma)
        if((pkpi[0] is None) or (pkpi[1] is None)): #{
          last_ln = None
          break
        else: #}{
          this_pkp = [pkp[pkpi[0]], pkp[pkpi[1]]]
          this_pkh = [x[this_pkp[0]], x[this_pkp[1]]]
          if((this_pkh[0] < last_pkh[0]) and (this_pkh[1] < last_pkh[1])): #{
            break
          else: #}{
            last_pkp = this_pkp
            last_pkh = this_pkh
          #}
        #}
      #}
      last_ln = ln
      vrbMsg(5, 'findProfileHip() ln = ' + str(ln) +
             ', last_pkp = [' + str(last_pkp[0]) + ', ' +
             str(last_pkp[1]) + '], ' +
             'last_pkh  = [' + str(last_pkh[0]) + ', ' +
             str(last_pkh[1]) + ']')
    #}
    if(not (last_ln is None)): #{
      rois['hip_r'] = [last_pkp[0], last_ln]
      rois['hip_l'] = [last_pkp[1], last_ln]
    #}
  #}
  vrbMsg(3, 'findProfileHip( rois[\'hip_l\'] = ' + str(rois['hip_l']))
  vrbMsg(3, 'findProfileHip( rois[\'hip_r\'] = ' + str(rois['hip_r']))
#}

def findProfileElbow(rois, ary): #{
  # Have found at least one shoulder and the centre of the pelvis. Find
  # elbows by tracking the first and last peaks down from the shoulder
  # finding the position of the maximum height of these peaks between
  # set fractions of the neck - pelvis distance down the arms.
  vrbMsg(3, 'findProfileElbow()')
  if(bool(rois['neck_cb']) and bool(rois['pelvis_c']) and
     (bool(rois['shoulder_l']) or bool(rois['shoulder_r']))): #{
    neck_ln = rois['neck_cb'][1]
    pelvis_ln = rois['pelvis_c'][1]
    neck_pelvis_dst = (pelvis_ln - neck_ln)
    elbow_s = int(neck_pelvis_dst * prm_elbow_s)
    elbow_e = int(neck_pelvis_dst * prm_elbow_e)
    discon  = [rois['shoulder_r'] is None, rois['shoulder_l'] is None]
    arm_max = [None, None]
    shoulder = [None if discon[0] else rois['shoulder_r'][:2],
                None if discon[1] else rois['shoulder_l'][:2]]
    arm_pkp = [None, None]
    finished = [False, False]
    d_to_s = [None, None]
    vrbMsg(5, 'findProfileElbow() neck_pelvis_dst = ' + str(neck_pelvis_dst) +
        ', neck_ln = ' + str(neck_ln) +
        ', elbow_e = ' + str(elbow_e) + 
        ', elbow_s = ' + str(elbow_s))
    srt_ln = neck_ln + int(neck_pelvis_dst * prm_elbow_i)
    end_ln = neck_ln + elbow_e
    for ln in range(srt_ln, end_ln, 1): #{
      vrbMsg(5, 'findProfileElbow() ln = ' + str(ln) +
          ', arm_pkp = ' + str(arm_pkp))
      x = ary[ln]
      pkp,pkw = findPeaks(x, prm_peak_sigma, prm_peak_minval, show = False)
      if(len(pkp) < 2): #{
        discon = [True, True]
        break
      #}
      # Check position of the first and last peak is contiguous with
      # the previous
      for i in range(0, 2): #{
        if(not (discon[i] or finished[i])): #{
          j = i * (len(pkp) - 1)
          if(arm_pkp[i] is None): #{
            arm_pkp[i] = pkp[j]
            continue
          #}
          d0 = abs(pkp[j] - arm_pkp[i])
          if((ln >= elbow_s) and (d0 > 3 * prm_peak_sigma)): #{
            # Discontinuity
            discon[i] = True
            vrbMsg(1, 'findProfileElbow() Arm ' + str(i) +
                  ' discontinuity at [' + str(pkp[j]) + ',' + str(ln) + '].')
            continue
          #}
          arm_pkp[i] = pkp[j]
          d_to_s[i] = dist(shoulder[i], [arm_pkp[i], ln])
          # Look for maximum intensity value
          # TODO Elbow may not be at maximum intensity. Look for elbow bend too.
          if(d_to_s[i] > elbow_e): #{
            finished[i] = True
          elif(d_to_s[i] > elbow_s): #}{
            p = x[arm_pkp[i]]
            if((arm_max[i] is None) or (p > arm_max[i][2])): #{
              arm_max[i] = [arm_pkp[i], ln, p]
            #}
          #}
        #}
        vrbMsg(5, 'findProfileElbow() ln = ' + str(ln) + 
            ' arm_max = ' + str(arm_max))
      #}
      if((discon[0] and discon[1]) or (finished[0] and finished[1])): #{
        break
      #}
    #}
    if(not (discon[0] or (arm_max[0] is None))): #{
      rois['elbow_r'] = [arm_max[0][0], arm_max[0][1]]
    #}
    if(not (discon[1] or (arm_max[1] is None))): #{
      rois['elbow_l'] = [arm_max[1][0], arm_max[1][1]]
    #}
  #}
  vrbMsg(3, 'findProfileElbow( rois[\'elbow_l\'] = ' + str(rois['elbow_l']))
  vrbMsg(3, 'findProfileElbow( rois[\'elbow_r\'] = ' + str(rois['elbow_r']))
#}

def findProfileWrist(rois, ary): #{
  # Have found at least one shoulder and elbow, the centre of the neck and
  # pelvis. # Find wrists by tracking down the arms from the elbows checking
  # for discontinuities and then looking for the minimum height peak that's
  # between set fractions of the uper arm length. The minimum will be above the
  # wrist so correct for this using an additional fraction of the elbow -
  # found wrist. The wrist minimum is a local minimum which implies need to
  # ensure this isn't over run, once mimimum has been held for set
  # number of lines then define it to be the local minimum and end search.
  # Discontinuity if there are less than two peaks of first and last peaks
  # run close to edge of image.
  vrbMsg(3, 'findProfileWrist()')
  if(bool(rois['shoulder_l']) and bool(rois['shoulder_r']) and
     (bool(rois['elbow_l']) or bool(rois['elbow_r']))): #{
    elbow_min = 0
    wrist_max = 0
    inset_lft = prm_peak_sigma # insets determine if peaks close to edge
    inset_rgt = len(ary[0]) - (inset_lft + 1)
    up_arm_len = [0, 0]
    elbow_s = [None, None]
    wrist_e = [None, None]
    arm_pkp = [None, None]
    arm_min = [None, None]
    wrist_min_cnt = [0, 0]
    discon  = [False, False]
    finished = [False, False]
    if(bool(rois['elbow_r']) and bool(rois['shoulder_r'])): #{
      elbow_s[0] = rois['elbow_r'][:2]
      elbow_min = elbow_s[0][1]
      up_arm_len[0] = dist(rois['shoulder_r'][:2], elbow_s[0])
      wrist_e[0] = rois['elbow_r'][1] + up_arm_len[0]
      wrist_max = wrist_e[0]
    #}
    if(bool(rois['elbow_l']) and bool(rois['shoulder_l'])): #{
      elbow_s[1] = rois['elbow_l'][:2]
      up_arm_len[1] = dist(rois['shoulder_l'][:2], elbow_s[1])
      wrist_e[1] = rois['elbow_l'][1] + up_arm_len[1]
      if(bool(rois['elbow_r'])): #{
        elbow_min = min(elbow_min, elbow_s[1][1])
        wrist_max = max(wrist_max, wrist_e[1])
      else: #}{
        elbow_min = elbow_s[1][1]
        wrist_max = wrist_e[1]
      #}
    #}
    elbow_min = int(elbow_min)
    wrist_max = int(wrist_max)
    last_arm_pkp = elbow_s[:]
    vrbMsg(7, 'findProfileWrist() elbow_min = ' + str(elbow_min) +
        ', wrist_max = ' + str(wrist_max))
    for ln in range(elbow_min, wrist_max, 1): #{
      vrbMsg(5, 'findProfileWrist() ln = ' + str(ln))
      x = ary[ln]
      pkp,pkw = findPeaks(x, prm_peak_sigma, prm_peak_minval, show = False)
      # If less than two peaks of peaks are at left or right edge of
      # image have discontinuity
      if(len(pkp) < 2): #{
        discon = [True, True]
      else: #}{
        if(pkp[0] < inset_lft): #{
          discon[0] = True
        #}
        if(pkp[1] > inset_rgt): #{
          discon[1] = True
        #}
      #}
      if(discon[0] and discon[1]): #{
        break
      #}
      # Check position of first and last peaks are contiguous with
      # the previous
      for i in range(0, 2): #{
        vrbMsg(7, 'findProfileWrist() discon = ' + str(discon) +
            ', finished = ' + str(finished))
        if(not ((elbow_s[i] is None) or discon[i] or finished[i])): #{
          j = i * (len(pkp) - 1)
          arm_pkp[i] = [pkp[j], ln]
          d0 = dist(last_arm_pkp[i], arm_pkp[i])
          vrbMsg(5, 'findProfileWrist() arm ' + str(i) + 
              ' arm_pkp ' + str(arm_pkp[i]) + ' d0 ' + str(d0))
          if(d0 > prm_peak_sigma): #{
            # Discontinuity
            vrbMsg(1, 'findProfileWrist() Arm ' + str(i) +
                ' discontinuity at [' + str(pkp[j]) + ',' + str(ln) + '].')
            discon[i] = True
            continue
          #}
          d1 = dist(elbow_s[i], arm_pkp[i])
          vrbMsg(5, 'findProfileWrist() arm ' + str(i) +
              ' d1 ' + str(d1) + ' up_arm_len ' + str(up_arm_len[i]))
          if(d1 > up_arm_len[i] * prm_wrist_e): #{
            finished[i] = True
          elif(d1 > up_arm_len[i] * prm_wrist_s): #}{
            vrbMsg(7, 'findProfileWrist() arm_min = ' + str(arm_min) +
              ', wrist_min_cnt = ' + str(wrist_min_cnt))
            if(arm_min[i] is None): #{
              arm_min[i] = [arm_pkp[i][0], arm_pkp[i][1],  x[arm_pkp[i][0]]]
            else: #}{
              m = x[arm_pkp[i][0]]
              if(m < arm_min[i][2]): #{
                arm_min[i] = [arm_pkp[i][0], arm_pkp[i][1], m]
                wrist_min_cnt[i] = 0
              else: #}{
                wrist_min_cnt[i] = wrist_min_cnt[i] + 1
                if(wrist_min_cnt[i] > prm_peak_sigma): #{
                  finished[i] = True
                #}
              #}
            #}
          #}
          last_arm_pkp[i] = arm_pkp[i]
        #}
        if(finished[0] and finished[1]): #{
          break
        #}
      #}
      vrbMsg(5, 'findProfileWrist() ln = ' + str(ln) +
          ' arm_min = ' + str(arm_min))
    #}
    for i in range(0, 2): #{
      if((not discon[i]) and (not (arm_min[i] is None))): #{
        elr = 'elbow_r' if (i == 0) else 'elbow_l'
        wlr = 'wrist_r' if (i == 0) else 'wrist_l'
        w = arm_min[i][:2]
        e = rois[elr]
        for j in range(0, 2): #{
          w[j] = int(w[j] + ((w[j] - e[j]) * prm_wrist_l))
        #}
        rois[wlr] = w
      #}
    #}
  #}  
  vrbMsg(3, 'findProfileWrist( rois[\'wrist_l\'] = ' + str(rois['wrist_l']))
  vrbMsg(3, 'findProfileWrist( rois[\'wrist_r\'] = ' + str(rois['wrist_r']))
#}

def findProfileKnee(rois, ary): #{
  # Have found at least one hip, the neck and the centre of the pelvis.
  # Find knees by tracking the hip peaks down from the hip finding
  # finding the position of the maximum height of these peaks between
  # set fractions of the neck - pelvis distance down the legs.
  vrbMsg(3, 'findProfileKnee()')
  if(bool(rois['neck_cb']) and bool(rois['pelvis_c']) and
     (bool(rois['hip_l']) or bool(rois['hip_r']))): #{
    neck_ln = rois['neck_cb'][1]
    pelvis_ln = rois['pelvis_c'][1]
    neck_pelvis_dst = (pelvis_ln - neck_ln)
    knee_s = int(neck_pelvis_dst * prm_knee_s)
    knee_e = int(neck_pelvis_dst * prm_knee_e)
    discon  = [rois['hip_r'] is None, rois['hip_l'] is None]
    leg_max = [None, None]
    hip = [None if discon[0] else rois['hip_r'][:2],
           None if discon[1] else rois['hip_l'][:2]]
    leg_pkp = [None, None]
    finished = [False, False]
    d_to_h = [None, None]
    vrbMsg(5, 'findProfileKnee() neck_pelvis_dst = ' + str(neck_pelvis_dst) +
        ', pelvis_ln = ' + str(pelvis_ln) +
        ', knee_e = ' + str(knee_e) + 
        ', knee_s = ' + str(knee_s))
    srt_ln = pelvis_ln
    end_ln = pelvis_ln + knee_e
    for ln in range(srt_ln, end_ln, 1): #{
      vrbMsg(5, 'findProfileKnee() ln = ' + str(ln) +
          ', leg_pkp = ' + str(leg_pkp))
      x = ary[ln]
      pkp,pkw = findPeaks(x, prm_peak_sigma, prm_peak_minval, show = False)
      if(len(pkp) < 2): #{
        discon = [not finished[0], not finished[1]]
        break
      #}
      # Get highest pair of peaks
      hi_pkp = findHighestPkPair(x, pkp)
      # Check positions of the highest peaks are contiguous with the previous
      for i in range(0, 2): #{
        if(not (discon[i] or finished[i])): #{
          j = hi_pkp[0] if(len(hi_pkp) == 1) else hi_pkp[i]
          leg_pkp[i] = pkp[j]
          d = abs(pkp[j] - leg_pkp[i])
          if((not discon[i]) and (ln >= knee_s) and (d > prm_peak_sigma)): #{
            # Discontinuity
            discon[i] = True
            vrbMsg(1, 'findProfileKnee() Leg ' + str(i) +
                  ' discontinuity at [' + str(pkp[j]) + ',' + str(ln) + '].')
          else: #}{
            d_to_h[i] = dist(hip[i], [leg_pkp[i], ln])
            if(d_to_h[i] > knee_e): #{
              finished[i] = True
            elif(d_to_h[i] > knee_s): #}{
              p = x[leg_pkp[i]]
              if((leg_max[i] is None) or (p > leg_max[i][2])): #{
                leg_max[i] = [leg_pkp[i], ln, p]
              #}
            #}
          #}
        #}
        vrbMsg(5, 'findProfileKnee() ln = ' + str(ln) + 
            ' leg_max = ' + str(leg_max))
      #}
      if((discon[0] and discon[1]) or (finished[0] and finished[1])): #{
        break
      #}
    #}
    if(not (discon[0] or (leg_max[0] is None))): #{
      rois['knee_r'] = [leg_max[0][0], leg_max[0][1]]
    #}
    if(not (discon[1] or (leg_max[1] is None))): #{
      rois['knee_l'] = [leg_max[1][0], leg_max[1][1]]
    #}
  #}
  vrbMsg(3, 'findProfileKnee( rois[\'knee_l\'] = ' + str(rois['knee_l']))
  vrbMsg(3, 'findProfileKnee( rois[\'knee_r\'] = ' + str(rois['knee_r']))
#}

def findProfileAnkle(rois, ary): #{
  # Have found at least one hip and the corresponding knee, the centre of
  # the neck and pelvis. # Find angles by tracking down the legs from the
  # knees looking for the minimum height peak that's between 1/2 and 1/1
  # of the upper leg length down from the knee.
  vrbMsg(3, 'findProfileAnkle()')
  if(bool(rois['hip_l']) and bool(rois['hip_r']) and
     (bool(rois['knee_l']) or bool(rois['knee_r']))): #{
    knee_min = 0
    up_leg_len = [0, 0]
    knee_s = [None, None]
    ankle_e = [None, None]
    leg_pkp = [None, None]
    leg_min = [None, None]
    discon  = [False, False]
    finished = [False, False]
    if(bool(rois['knee_r']) and bool(rois['hip_r'])): #{
      knee_s[0] = rois['knee_r'][:2]
      knee_min = knee_s[0][1]
      up_leg_len[0] = dist(rois['hip_r'][:2], knee_s[0])
      ankle_e[0] = rois['knee_r'][1] + up_leg_len[0]
      ankle_max = ankle_e[0]
    #}
    if(bool(rois['knee_l']) and bool(rois['hip_l'])): #{
      knee_s[1] = rois['knee_l'][:2]
      up_leg_len[1] = dist(rois['hip_l'][:2], knee_s[1])
      ankle_e[1] = rois['knee_l'][1] + up_leg_len[1]
      if(bool(rois['knee_r'])): #{
        knee_min = min(knee_min, knee_s[1][1])
        ankle_max = max(ankle_max, ankle_e[1])
      else: #}{
        knee_min = knee_s[1][1]
        ankle_max = ankle_e[1]
      #}
    #}
    knee_min = int(knee_min)
    ankle_max = int(ankle_max)
    last_leg_pkp = knee_s[:]
    for ln in range(knee_min, ankle_max, 1): #{
      vrbMsg(5, 'findProfileAnkle() ln = ' + str(ln))
      x = ary[ln]
      pkp,pkw = findPeaks(x, prm_peak_sigma, prm_peak_minval, show = False)
      if(len(pkp) < 2): #{
        discon = [True, True]
        break
      #}
      # Check position of first and last peaks are contiguous with
      # the previous
      for i in range(0, 2): #{
        if(not ((knee_s[i] is None) or discon[i])): #{
          j = i * (len(pkp) - 1)
          leg_pkp[i] = [pkp[j], ln]
          d0 = dist(last_leg_pkp[i], leg_pkp[i])
          vrbMsg(5, 'findProfileAnkle() leg ' + str(i) + 
              ' leg_pkp ' + str(leg_pkp[i]) + ' d0 ' + str(d0))
          if(d0 > prm_peak_sigma): #{
            # Discontinuity
            vrbMsg(1, 'findProfileAnkle() Leg ' + str(i) +
                ' discontinuity at [' + str(pkp[j]) + ',' + str(ln) + '].')
            discon[i] = True
            break
          else: #}{
            d1 = dist(knee_s[i], leg_pkp[i])
            vrbMsg(5, 'findProfileAnkle() leg ' + str(i) +
                ' d1 ' + str(d1) + ' leg_len ' + str(up_leg_len[i]))
            if(d1 > up_leg_len[i]): #{
              finished[i] = True
            elif(d1 > up_leg_len[i] * 0.5): #}{
              if(leg_min[i] is None): #{
                leg_min[i] = [leg_pkp[i][0], leg_pkp[i][1],  x[leg_pkp[i][0]]]
              else: #}{
                m = x[leg_pkp[i][0]]
                if(m < leg_min[i][2]): #{
                  leg_min[i] = [leg_pkp[i][0], leg_pkp[i][1], m]
                #}
              #}
            #}
            last_leg_pkp[i] = leg_pkp[i]
          #}
        #}
        if(finished[0] and finished[1]): #{
          break
        #}
      #}
      vrbMsg(5, 'findProfileAnkle() ln = ' + str(ln) +
          ' leg_min = ' + str(leg_min))
    #}
    if((not discon[0]) and (not (leg_min[0] is None))): #{
      rois['ankle_r'] = leg_min[0][:2]
    #}
    if((not discon[1]) and (not (leg_min[1] is None))): #{
      rois['ankle_l'] = leg_min[1][:2]
    #}
  #}  
  vrbMsg(3, 'findProfileAnkle( rois[\'ankle_l\'] = ' + str(rois['ankle_l']))
  vrbMsg(3, 'findProfileAnkle( rois[\'ankle_r\'] = ' + str(rois['ankle_r']))
#}

def findROICentres(smt_file): #{
  """
  Find the centres of regions of interest in the scan image.

  Args:
  smt_file:    File with smoothed (Gaussian convolved) scan image

  Returns:      Regions of interest as an ordered array of integer pairs.
  """
  rois = []
  ary = None
  state = 0
  roi_centres = { 
      'head_ct': None,
      'neck_cb': None,
      'shoulder_l': None,
      'shoulder_r': None,
      'elbow_l': None,
      'elbow_r': None,
      'pelvis_c': None,
      'wrist_l': None,
      'wrist_r': None,
      'hip_l': None,
      'hip_r': None,
      'knee_l': None,
      'knee_r': None,
      'ankle_l': None,
      'ankle_r': None,
      'toes_e': None}

  vrbMsg(1, 'findROICentres() smt_file = ' + smt_file)
  # Read the smoothed Woolz object and create a NumPy array from it
  err_num, smt_obj = readWoolzObj(smt_file)
  if(not bool(err_num)): #{
	  err_num, org, ary = wlzObjToNP(smt_obj)
  #}
  sz = np.shape(ary)
  vrbMsg(5, 'findROICentres() object size = ' + str(sz))
  # Work down the scan finding coordinates, ordering has dependency
  # but this is checked in the individual functions
  if((sz[0] >= img_size_min) and (sz[1] >= img_size_min)): #{
    findProfileHeadCT(roi_centres, ary)
    findProfileToes(roi_centres, ary)
  #}
  if(bool(roi_centres['toes_e']) and bool(roi_centres['head_ct'])): #{
    findProfileShoulders(roi_centres, ary)
    findProfilePelvis(roi_centres, ary)
    findProfileHip(roi_centres, ary)
    findProfileElbow(roi_centres, ary)
    findProfileWrist(roi_centres, ary)
    findProfileKnee(roi_centres, ary)
    findProfileAnkle(roi_centres, ary)
    rois = roi_centres
    for cen in rois: #{
      # numpy gives int64 which is not always handled (eg by json) so convert
      pos = rois[cen]
      if(not (pos is None)): #{
        rois[cen] = [int(pos[0]), int(pos[1])]
      #}
    #}
  #}
  vrbMsg(1, 'findROICentres() rois = ' + str(rois))
  return rois 
#}

def saveROICentres(roc_file, roi_centres): #{
  vrbMsg(1, 'saveROICentres() file = ' + roc_file + ', roi_centres = {...}')
  status = 0
  try: #{
    with open(roc_file, 'w') as f: #{
      json.dump(roi_centres, f)
  #}
  except Exception as e: # }{
    vrbMsg(1, 'saveROICentres() exception: ' + str(e))
    status = 1
  #}
  vrbMsg(1, 'saveROICentres() status = ' + str(status))
  return status 
#}

def extractROIFromCentres(roi_file_base, obj_file, roi_centres): #{
  vrbMsg(1, 'extractROIFromCentres() roi_file_base = ' + roi_file_base +
         ', obj_file = ' + obj_file +
         ', roi_centres = {...}')
  status = 0
  # Read Woolz object
  err_num, base_obj = readWoolzObj(obj_file)
  if(bool(err_num)): #{
    vrbMsg(1,
        'extractROIFromCentres() Failed to read Woolz object from file ' +
        obj_file)
    status = 1
  #}
  # For each ROI centre extract region about this centre
  if(status == 0): #{
    for roc in roi_centres: #{
      cen = roi_centres[roc]
      vrbMsg(3, 'extractROIFromCentres() ' + str(roc) + ' ' + str(cen))
      if(cen is None): #{
        vrbMsg(3, 'extractROIFromCentres() ' + str(roc) + ' ROI not found')
      else: #}{
        status, roi_obj = extractROI(base_obj, roc, cen, roi_obj_size,
           prm_filter_sobel, 1 if prm_filter_window else 0)
        # Save region to file
        if((status == 0) and bool(roi_obj)): #{
          roi_file = roi_file_base + '-' + roc + '.wlz'
          err_num = writeWoolzObj(roi_file, roi_obj)
          if(bool(err_num)): #{
            status = 1
          #}
          if(bool(roi_obj)): #{
            dummy = w.WlzFreeObj(roi_obj)
          #}
        #}
        if(status != 0): #{
          vrbMsg(1,
              'extractROIFromCentres() Failed to write Woolz object to file ' +
              roi_file)
        #}
      #}
    #}
  #}
  if(bool(base_obj)): #{
    dummy = w.WlzFreeObj(base_obj)
  #}
  vrbMsg(1, 'extractROIFromCentres() status = ' + str(status))
  return status
#}

def extractROI(base_obj, roc, cen, roi_sz,
    filter_sobel=False, filter_smooth=True, filter_window=0): #{
  """
  Extract a filtered region of interest (ROI) from the given Woolz object.

  Args:
  base_obj:     Given Woolz object.
  roc:          Centre of ROI within the given object.
  roi_sz:       ROI width ( = ROI height)
  filter_sobel: Use Sobel filter to enhance edges if True.
  filter_smooth: Apply Gaussian smoothing after Sobel filter.
  filter_window: Use (Parzen) window function to enhance centre.
                This is a bit mask: with value: 0, no window function;
                bit 0 (1) set, window at extraction and bit 1 (2) set
                window after optional Sobel filter.

  Returns:      Status, ROI object
  """
  vrbMsg(1, 'extractROI() base_obj = ' + str(base_obj) +
         ', roc = ' + str(roc) +
         ', cen = ' + str(cen) +
         ', roi_sz = ' + str(roi_sz) +
         ', filter_sobel = ' + str(filter_sobel) +
         ', filter_smooth = ' + str(filter_smooth) +
         ', filter_window = ' + str(filter_window))
  status = 0
  regions_to_sep = ['ankle_l', 'ankle_r', 'elbow_l', 'elbow_r',
                    'hip_l', 'hip_r', 'knee_l', 'knee_r',
                    'wrist_l', 'wrist_r']
  cutBox = w.WlzIBox2()
  cutBox.xMin = int(cen[0] - (roi_sz / 2))
  cutBox.yMin = int(cen[1] - (roi_sz / 2))
  cutBox.xMax = roi_sz + cutBox.xMin - 1
  cutBox.yMax = roi_sz + cutBox.yMin - 1
  err_num = c.c_int(w.WLZ_ERR_NONE)
  roi_obj = w.WlzAssignObject(
      w.WlzCutObjToBox2D(base_obj, cutBox,
          c.c_int(w.WLZ_GREY_UBYTE), c.c_int(0),
          c.c_double(0.0), c.c_double(0.0), c.byref(err_num)), None)
  vrbMsg(5, 'extractROI() roi_obj = ' + str(roi_obj) +
      ', err_num = ' + str(err_num))
  # Window the region to reduce edge effects
  if(not bool(err_num)): #{
    if(roc in regions_to_sep): #{
      # If this ROI is appropriate select just that part of the object
      # which contains the ROI centre
      # Smooth ROI object and normalise it
      vrbMsg(5, 'extractROI() roc = ' + roc + ' - separating regions')
      obj0 = None
      obj1 = None
      obj2 = None
      lbl_ary = None
      max_rgn = 1000
      n_rgn = c.c_int(0)
      sigma = w.WlzDVertex3()
      sigma.vtX = c.c_double(prm_peak_sigma / 2)
      sigma.vtY = c.c_double(prm_peak_sigma / 2)
      sigma.vtZ = c.c_double(1.0)
      order = w.WlzIVertex3()
      order.vtX = c.c_int(0)
      order.vtY = c.c_int(0)
      order.vtZ = c.c_int(0)
      direc = w.WlzIVertex3()
      direc.vtX = c.c_int(1)
      direc.vtY = c.c_int(1)
      direc.vtX = c.c_int(0)
      obj0 = w.WlzAssignObject(
          w.WlzGaussFilter(roi_obj, sigma, order, direc,
             c.c_int(w.WLZ_GREY_INT),
             c.c_int(w.ALG_PAD_ZERO), c.c_double(0.0), c.c_int(0),
             c.byref(err_num)), None)
      vrbMsg(5, 'extractROI() obj0 = ' + str(obj1) +
          ', err_num = ' + str(err_num))
      if(not bool(err_num)): #{
        err_num = w.WlzGreyNormalise(obj0, 0)
        vrbMsg(5, 'extractROI() err_num = ' + str(err_num))
      #}
      if((not bool(err_num)) and bool(filter_window & 1)): #{
        org = w.WlzIVertex2()
        org.vtX = int(cen[0])
        org.vtY = int(cen[1])
        rad = w.WlzIVertex2()
        rad.vtX = int((roi_sz / 2))
        rad.vtY = int((roi_sz / 2))
        err_num = c.c_int(w.WLZ_ERR_NONE)
        obj1 = w.WlzAssignObject(
            w.WlzWindow(obj0, c.c_int(w.WLZ_WINDOWFN_PARZEN),
                org, rad, c.byref(err_num)), None)
        vrbMsg(5, 'extractROI() WlzWindow() err_num = ' + str(err_num))
        if(not bool(err_num)): #{
          dummy = w.WlzFreeObj(obj0)
          obj0 = obj1
        #}
      #}
      # Threshold the object to get foreground
      if(not bool(err_num)): #{
        hilo = w.WLZ_THRESH_HIGH
        tv = w.WlzPixelV()
        tv = c.c_int(prm_peak_minval * 2)
        err_num = c.c_int(w.WLZ_ERR_NONE)
        obj1 = w.WlzAssignObject(
            w.WlzThresholdI(obj0, tv, hilo, c.byref(err_num)), None)
        vrbMsg(5, 'extractROI() obj2 = ' + str(obj1) +
            ', err_num = ' + str(err_num))
      #}
      dummy = w.WlzFreeObj(obj0)
      obj0 = None
      # Perform connected component labeling
      if(not bool(err_num)): #{
        lbl_ary = c.cast(w.AlcCalloc(max_rgn, c.sizeof(w.WlzObject)),
        c.POINTER(c.POINTER(w.WlzObject)))
        err_num = w.WlzLabel(obj1, c.byref(n_rgn), c.byref(lbl_ary),
            max_rgn, 1, c.c_int(w.WLZ_8_CONNECTED))
      #}
      dummy = w.WlzFreeObj(obj1)
      obj1 = None
      # Find which component contains the ROI centre
      n_rgn = n_rgn.value
      for i in range(0, n_rgn): #{
        vrbMsg(5, 'extractROI() lbl_ary[' + str(i) + '] = ' +
            str(lbl_ary[i]))
        err_num = c.c_int(w.WLZ_ERR_NONE)
        hit = w.WlzInsideDomain(lbl_ary[i], 0, cen[1], cen[0],
            c.byref(err_num))
        if(bool(hit)): #{
          vrbMsg(5, 'extractROI() hit')
          obj2 = w.WlzAssignObject(
              w.WlzGreyTransfer(lbl_ary[i], roi_obj, 0,
                  c.byref(err_num)), None)
          dummy = w.WlzFreeObj(roi_obj)
          roi_obj = w.WlzAssignObject(
              w.WlzCutObjToBox2D(obj2, cutBox,
                  c.c_int(w.WLZ_GREY_UBYTE), c.c_int(0),
                  c.c_double(0.0), c.c_double(0.0),
                  c.byref(err_num)), None)
          dummy = w.WlzFreeObj(obj2)
          obj2 = None
          break
        #}
      #}
      if(not(lbl_ary is None)): #{
        for i in range(0, n_rgn): #{
          if(bool(lbl_ary[i])): #{
            dummy = w.WlzFreeObj(lbl_ary[i])
          #}
        #}
        w.AlcFree(lbl_ary)
      #}
    #}
  #}
  # Apply Sobel filter so that we match to edges
  if((not bool(err_num)) and filter_sobel): #{
    err_num = c.c_int(w.WLZ_ERR_NONE)
    obj1 = w.WlzAssignObject(
        w.WlzSobel(roi_obj, 1, 1, c.byref(err_num)), None)
    vrbMsg(5, 'extractROI() WlzSobel() err_num = ' + str(err_num))
    if(not bool(err_num)): #{
      dummy = w.WlzFreeObj(roi_obj)
      roi_obj = obj1
    #}
  #}
  if((not bool(err_num)) and filter_smooth): #{
    sigma = w.WlzDVertex3()
    sigma.vtX = 1.0
    sigma.vtY = 1.0
    sigma.vtZ = 0.0
    order = w.WlzIVertex3()
    order.vtX = 0
    order.vtY = 0
    order.vtZ = 0
    direc = w.WlzIVertex3()
    direc.vtX = 1
    direc.vtY = 1
    direc.vtZ = 0
    err_num = c.c_int(w.WLZ_ERR_NONE)
    obj1 = w.WlzAssignObject(
        w.WlzGaussFilter(roi_obj, sigma, order, direc,
            c.c_int(w.WLZ_GREY_ERROR), c.c_int(w.ALG_PAD_NONE), 0.0, 0,
            c.byref(err_num)), None)
    vrbMsg(5, 'extractROI() WlzGaussFilter() err_num = ' + str(err_num))
    if(not bool(err_num)): #{
      dummy = w.WlzFreeObj(roi_obj)
      roi_obj = obj1
    #}
  #}
  # Apply window function again if required
  if((not bool(err_num)) and bool(filter_window & 2)): #{
    org = w.WlzIVertex2()
    org.vtX = int(cen[0])
    org.vtY = int(cen[1])
    rad = w.WlzIVertex2()
    rad.vtX = int((roi_sz / 2))
    rad.vtY = int((roi_sz / 2))
    err_num = c.c_int(w.WLZ_ERR_NONE)
    obj0 = w.WlzAssignObject(
        w.WlzWindow(roi_obj, c.c_int(w.WLZ_WINDOWFN_PARZEN),
            org, rad, c.byref(err_num)), None)
    vrbMsg(5, 'extractROI() WlzWindow() err_num = ' + str(err_num))
    if(not bool(err_num)): #{
      dummy = w.WlzFreeObj(roi_obj)
      roi_obj = obj0
    #}
  #}
  # Shift the region to the origin
  if(not bool(err_num)): #{
    x_shift = -(roi_obj.contents.domain.i.contents.kol1)
    y_shift = -(roi_obj.contents.domain.i.contents.line1)
    err_num = c.c_int(w.WLZ_ERR_NONE)
    obj1 = w.WlzAssignObject(
        w.WlzShiftObject(roi_obj,
            x_shift, y_shift, 0, c.byref(err_num)), None)
    vrbMsg(5, 'extractROI() WlzShiftObject() err_num = ' + str(err_num))
    if(not bool(err_num)): #{
      dummy = w.WlzFreeObj(roi_obj)
      roi_obj = obj1
    #}
  #}
  vrbMsg(1, 'extractROI() status = ' + str(status) +
         ', roi_obj = ' + str(roi_obj))
  return status, roi_obj
#}

def processImage(model, image, fmt): #{
  vrbMsg(1, 'processImage() model = ' + model +
         ', image = ' + image + 
         ', fmt = ' + fmt)
  status = 0
  # Convert image to a temporary Woolz image
  tmp_file = tmpFileBase()
  bfi_file = bfiFileBase(image)
  hst_file = hstFileBase(model)
  smt_file = smtFileBase(image)
  roc_file = roiCentresFile(image)
  roi_file_base = roiFileBase(image)
  vrbMsg(5, 'processImage() tmp_file = ' + tmp_file)
  status,dummy = runCmd('convert ' + args.scan_dir + '/' + 
      image + '.' + fmt + ' ' + tmp_file + '0.tif')
  if(status == 0): #{
    status,dummy = runCmd('WlzExtFFConvert -f tif -F wlz ' +
        '-o ' + tmp_file + '0.wlz ' +
        tmp_file + '0.tif')
  #}
  # Check Woolz image grey type, if required convert from colour
  status,stats = runCmd('WlzGreyStats ' + tmp_file + '0.wlz ')
  if(status == 0): #{
    gType = stats.split(' ')[1]
    if(gType == 'WLZ_GREY_RGBA'): #{
      status,dummy = runCmd('WlzRGBAConvert ' +
          wlz_rgba_channel[args.colour_channel] + ' ' +
          tmp_file + '0.wlz > ' + tmp_file + '1.wlz')
      if(status == 0): #{
        status,dummy = runCmd('mv ' + tmp_file + '1.wlz ' + tmp_file + '0.wlz')
      #}
    #}
  #}
  # Make sure the image is a normalised unsigned byte image
  if(status == 0): #{
    status,dummy = runCmd('WlzGreyNormalise -u ' + 
        tmp_file + '0.wlz > ' + tmp_file + '1.wlz')
    if(status == 0): #{
      status,dummy = runCmd('mv ' + tmp_file + '1.wlz ' + tmp_file + '0.wlz')
    #}
  #}
  # Fill background of image with black rather than white
  if(args.background_valid): #{
    status,dummy = runCmd('mv ' + tmp_file + '0.wlz ' + tmp_file + '1.wlz')
  else: #}{
    status = fillBackground(tmp_file, 0, 1)
  #}
  # Move the filled image to the working directory.
  status,dummy = runCmd('mv ' + tmp_file + '1.wlz ' + bfi_file + '.wlz')
  if(status == 0): #{
    dummy = runCmd('rm -f ' + tmp_file + '*.*')
  #}
  # Histogram match to the model
  if(status == 0): #{
    if(image == model): #{
      status,dummy = runCmd('WlzHistogramObj -o ' + hst_file + '.wlz ' +
          bfi_file + '.wlz')
    else: #}{
      status,dummy = runCmd('WlzHistogramMatchObj -m ' + hst_file + '.wlz ' +
          ' -o ' + tmp_file + '1.wlz ' + bfi_file + '.wlz')
      if(status == 0): #{
        status,dummy = runCmd('mv ' + tmp_file + '1.wlz ' + bfi_file + '.wlz')
      #}
    #}
  #}
  # Create NIfTI file
  if(status == 0): #{
    status, dummy = runCmd('WlzExtFFConvert -f wlz -F nii -o ' +
        bfi_file + '.nii ' + bfi_file + '.wlz')
  #}
  # Smooth the image
  if(status == 0): #{
    status, dummy = runCmd('WlzSepFilterObj -g u -m ' +
        str(prm_peak_sigma) + ',' + str(prm_peak_sigma) + ' -Pe -t x,y ' +
        bfi_file + '.wlz | WlzGreyNormalise -u > ' + smt_file + '.wlz')
  #}
  # Find ROI centres
  if(status == 0): #{
    roi_centres = findROICentres(smt_file + '.wlz')
  #}
  # Save ROI centres
  if(status == 0): #{
    status = saveROICentres(roc_file, roi_centres)
  #}
  # Extract ROI images as Woolz and NIfTI
  if(status == 0): #{
    status = extractROIFromCentres(roi_file_base,
        bfi_file + '.wlz', roi_centres)
  #}
  vrbMsg(1, 'processImage() status = ' + str(status))
  return status 
#}

def registerAssay(model, assay): #{
  vrbMsg(1, 'registerAssay() model = ' + model + ', assay = ' + assay)
  status = 0
  roc_lists = [None, None]
  # Read roi centres for the model and assay
  assay_pair = [model, assay]
  for idx in range(0, len(assay_pair)): #{
    roc_file = roiCentresFile(assay_pair[idx])
    vrbMsg(5, 'registerAssay() roc_file = ' + roc_file)
    status,roc_lists[idx] = readJsnFile(roc_file)
    if(bool(status)): #{
      errMsg('Failed to read region of interest centres file: ' + roc_file)
      break
    #}
  #}
  # Register each region of interest pair using ANTs
  if(status == 0): #{
    for roc in roc_lists[0]: #{
      roi_file_pair = [None, None]
      vrbMsg(1, 'registerAssay() model = ' + model +
          ', assay = ' + assay + ', roc = ' + roc)
      if((roc_lists[0][roc] is None) or (roc_lists[1][roc] is None)): #{
        vrbMsg(1, 'registerAssay() no matching pair for model = ' + model +
            ', assay = ' + assay + ', roc = ' + roc)
      else: #}{
        # Convert Woolz file to NIfTI
        for idx in range(0, len(assay_pair)): #{
          roi_file_pair[idx] = roiFileBase(assay_pair[idx]) + '-' + roc
          status,dummy = runCmd('WlzExtFFConvert -f wlz -F nii -o ' + 
              roi_file_pair[idx] + '.nii ' + roi_file_pair[idx] + '.wlz')
          if(status > 0): #{
            errMsg('Failed to convert ' + roi_file_pair[idx] + '.wlz to NIfTI')
            break
          #}
        #}
        if(status > 0): #{
          break
        #}
        trx_base = trxFileBase(model, assay, roc)
        reg_cmd = 'ANTS 2'
        if(isVerbose(3)): #{
          reg_cmd = reg_cmd + ' -v'
        #}
        reg_cmd = (reg_cmd + 
            ' -o ' + trx_base +
            ' -m MI\\[' + roi_file_pair[0] + '.nii,' +
                roi_file_pair[1] + '.nii, 1, 5\\]' +
            ' -t SyN\\[0.5\\] -r Gauss\\[10,5\\]' +
            ' -i 500x100x100 --subsampling-factors 16x4x1' +
            ' --gaussian-smoothing-sigmas 8x4x1' +
            ' --rigid-affine')
        status,dummy = runCmd(reg_cmd)
        if(status > 0): #{
          errMsg('Failed to register regions of interest ' +
              roi_file_pair[0] + '.nii,' + roi_file_pair[1] + '.nii')
          break
        #}
        # Create registered image for checking registration
        reg_cmd = 'antsApplyTransforms'
        if(isVerbose(3)): #{
          reg_cmd = reg_cmd + ' -v'
        #}
        reg_file_base = regFileBase(model, assay, roc)
        reg_cmd = (reg_cmd +
            ' -d 2' +
            ' -i ' + roi_file_pair[1] + '.nii' +
            ' -r ' + roi_file_pair[0] + '.nii' +
            ' -o ' + reg_file_base + '.nii' +
            ' -t ' + trx_base + 'Affine.txt' +
            ' -t ' + trx_base + 'Warp.nii.gz') # 
        status,dummy = runCmd(reg_cmd)
        if(status > 0): #{
          errMsg('Failed to apply registration for ' + reg_file_base + '.nii')
          break
        #}
        reg_cmd = ('WlzExtFFConvert -o - -F wlz ' + reg_file_base + '.nii ' +
            ' | WlzGreyNormalise -u > ' + reg_file_base + '.wlz')
        status,dummy = runCmd(reg_cmd)
        if(status == 0): #{
          reg_cmd = ('WlzExtFFConvert -o ' + reg_file_base + '.jpg ' +
              reg_file_base + '.wlz')
          status,dummy = runCmd(reg_cmd)
        #}
        if(status > 0): #{
          errMsg('Failed to convert ' + reg_file_base +
              '.nii to Woolz and Jpeg')
          break
        #}
      #}
    #}
  #}
  vrbMsg(1, 'registerAssay() status = ' + str(status))
  return status 
#}

def readModelPoints(model, model_pts_file): #{
  vrbMsg(1, 'readModelPoints() model = ' + model + 
      ', model_pts_file = ' + model_pts_file)
  status = 0
  status,model_pts = readJsnFile(model_pts_file)
  if(status == 0): #{
    if(not (('points' in model_pts) and
            (len(model_pts['points']) > 1))): #{
      vrbMsg(1, 'readModelPoints() model_pts_file invalid')
      status = 1
    #}
  #}
  # Ensure points have 4 components x, y, confidence value (1.0 for model)
  # and model number
  if(status == 0): #{
    z = [0, 0, 1.0, model]
    points = model_pts['points']
    for pnt in points: #{
      l = len(points[pnt])
      if(l > 4): #{
        points[pnt] = points[pnt][:4]
      elif(l < 4): #}{
        points[pnt].extend(z[l:4])
      #}
    #}
  #}
  vrbMsg(1, 'readModelPoints() model_pts = ' + str(model_pts) + 
      ', status = ' + str(status))
  return model_pts,status 
#}

  """
  Return measurements read from the file with the given name. If the file
  extension is  '.cvs' then a comma seperated file in which the first record's
  fields are the measurement keys and subsequent records are the measurement
  values (corresponding to the keys) for each assay otherwise the file is
  assumed to be a valid json file.
  Example for csv file:
    assay,pixelsz_x,pixelsz_y,ulna_rd_x,ulna_rd_y,ulna_rd_c,ulna_rd_m,...
    1020654,2.0,2.0,62.0,420.0,0.707068367486,4775892,...
    775892,2.0,2.0,211.0,656.0,0.703338260521,4775892,...
    ...
  
  Args:
  filename:     Full path of file to read measurements from.
  force_mod:    If bool(force_mod) is true the force the model fields
                to this value.

  Returns:      Measurements and status.
  """
def readMeasurements(filename, mod): #{
  status = 0
  meas = []
  vrbMsg(1, 'readMeasurements() filename = ' + filename)
  try: #{
    ln = 0
    z = [0, 0, 1.0, 0]
    s0 = filename.split('.')
    if((len(s0) > 1) and (s0[-1] == 'csv')): #{
      keys = []
      components = {'x': 0, 'y': 1, 'c': 2, 'm': 3}
      with open(filename) as f: #{
        reader = csv.reader(f, delimiter=',')
        for row in reader: #{
          if(ln == 0): #{
            keys =  row
          else: #}{
            meas_asy = {'assay': 0, 'pixelsz': [-1.0,-1.0], 'points': {}}
            for i in range(0, len(keys)): #{
              if(keys[i] == 'assay'): #{
                meas_asy['assay'] = row[i]
              elif(keys[i] == 'pixelsz_x'): #}{
                meas_asy['pixelsz'][0] = row[i]
              elif(keys[i] == 'pixelsz_y'): #}{
                meas_asy['pixelsz'][1] = row[i]
              else: #}{
                ks = keys[i].split('_')
                kc = ks[-1]
                ky = '_'.join(ks[:-1])
                if(not (ky in meas_asy['points'])): #{
                  m = int(mod) if bool(mod) else 0
                  meas_asy['points'][ky] = [0.0, 0.0, 0.0, m]
                #} 
                if(kc == 'm'): #{
                  m = int(mod) if bool(mod) else int(row[i])
                  meas_asy['points'][ky][components[kc]] = m
                else: #}{
                  meas_asy['points'][ky][components[kc]] = float(row[i])
                #}
              #}
            #}
            meas.append(meas_asy)
          #}
          ln += 1
        #}
      #}
    else: #}{
      errMsg('Reading measurements from json files not implemented.') # TODO
      status = 1
    #}
  except (IOError, cvs.Error): #}{
    status = 1
  #}
  vrbMsg(1, 'readMeasurements() meas = ' + str(meas) +
      ', status = ' + str(status))
  return meas,status
#}

def modelMeasurements(model, model_pts_file, scan_dir, in_fmt): #{
  vrbMsg(1, 'modelMeasurements() model = ' + model +
      ', model_pts_file = ' + model_pts_file +
      ', scan_dir = ' + scan_dir +
      ',in_fmt = ' + in_fmt)
  model_meas,status = readModelPoints(model, model_pts_file)
  if(status == 0): #{
    pixelsz = getPixelSz(model, scan_dir, in_fmt)
    model_meas.update({'assay': model, 'pixelsz': pixelsz})
  #}
  vrbMsg(1, 'modelMeasurements() model_meas = ' + str(model_meas) +
      ', status = ' + str(status))
  return model_meas,status
#}

def combineMeasurements(gvn_meas, min_cnf): #{
  vrbMsg(1, 'combineMeasurements() gvn_meas, min_cnf = ' +str(min_cnf))
  status = 0
  cmb_meas = []
  n_cmb_meas = 0
  n_gvn_meas = len(gvn_meas)
  if(n_gvn_meas > 0): #{
    for ig in range(0, n_gvn_meas): #{ for each assay of the given measurements
      vrbMsg(3, 'combineMeasurements() measurements index = ' + str(ig))
      gvn_m = gvn_meas[ig]    # current measurements
      gvn_a = gvn_m['assay']  # current given assay current measurements
      gvn_p = gvn_m['points'] # current given points current measurements
      cmb_m = None            # combined measurements
      # Find combined assay corresponding to the current given one
      cmb_a = None
      n_cmb_meas = len(cmb_meas)
      for ic in range(0, n_cmb_meas): #{
        cmb_m = cmb_meas[ic]
        cmb_a = cmb_m['assay']
        if(cmb_a == gvn_a): #{
          break
        #}
      #}
      vrbMsg(3, 'combineMeasurements() assay = ' + str(gvn_a))
      if(cmb_a == gvn_a): #{
        # Found assay in combined measurements so update it's values
        vrbMsg(3, 'combineMeasurements() assay ' + str(gvn_a) + ' in combined')
        vrbMsg(7, 'combineMeasurements() given points = ' + str(gvn_p))
        cmb_p = cmb_m['points']
        for p in cmb_p: #{
          if(p in gvn_p): #{
            # Point of combined assay is in given assay
            gvn_x = gvn_p[p] # Given point
            cmb_x = cmb_p[p] # Combined point
            # If confidence value greater than threshold and greater than
            # that of the combined point update the combined point
            if((not (gvn_x[2] < min_cnf)) and (cmb_x[2] < gvn_x[2])): #{
              cmb_p[p] = cp.copy(gvn_x)
            #}
          else: #}{
            # Point not in combined measurements for assay so add it
            gvn_x = gvn_p[p] # Given point
            cmb_p[p] = thresholdMeasuredPoints(gvn_x, min_cnf)
          #}
        #}
        vrbMsg(7, 'combineMeasurements() combined points = ' + str(cmb_p))
      else: #}{
        # Given measurements assay not in combined measurements so add it
        vrbMsg(3, 'combineMeasurements() assay ' + str(gvn_a) +
            ' not in combined')
        vrbMsg(7, 'combineMeasurements() given points = ' + str(gvn_p))
        cmb_z = cp.copy(gvn_m['pixelsz'])
        cmb_p = cp.copy(gvn_p)
        for p in cmb_p: #{
          cmb_p[p] = thresholdMeasuredPoints(cmb_p[p], min_cnf)
        #}
        vrbMsg(7, 'combineMeasurements() combined points = ' + str(cmb_p))
        cmb_meas.append({'assay': gvn_a, 'pixelsz': cmb_z, 'points': cmb_p})
      #}
    #}
  #}
  vrbMsg(1, 'combineMeasurements()')
  return cmb_meas
#}

def thresholdMeasuredPoints(gvn_x, min_cnf): #{
  thr_x = [0.0, 0.0, 0.0, 0]
  if(not (gvn_x is None)): #{
    if(not (gvn_x[2] < min_cnf)): #{
      thr_x = cp.copy(gvn_x)
    #}
  #}
  return thr_x
#}

def mapMeasurements(model, model_pts_file, assaylist, scan_dir, in_fmt,
    comp_cnf = True): #{
  vrbMsg(1, 'mapMeasurements() model = ' + model +
      ', model_pts_file = ' + model_pts_file +
      ', assaylist = ' + str(assaylist) +
      ', scan_dir = ' + scan_dir +
      ', in_fmt = ' + in_fmt +
      ', comp_cnf = ' + str(comp_cnf))
  meas = []
  roc_list = []
  idx_obj_base = tmpFileBase() + '-roi-idx'
  tmp_reg_base = tmpFileBase() + '-reg-idx'
  # Read model points
  model_pts,status = readModelPoints(model, model_pts_file)
  if(status != 0): #{
    errMsg('Failed to read model points from file: ' + model_pts_file)
  #}
  # Read model centres file
  if(status == 0): #{
    model_roc_file = roiCentresFile(model)
    status,model_roc = readJsnFile(model_roc_file)
    if(bool(status)): #{
      errMsg('Failed to read model region centres from file: ' + model_pts_file)
    #}
  #}
  # Create index object for all assay point locations
  if(status == 0): #{
    cmd = ('WlzMakeRect -x 0,' + str(roi_obj_size) +
           ' -y 0,' + str(roi_obj_size) +
           ' | WlzGreySetIncValues -o ' + idx_obj_base + '.wlz')
    status,dummy = runCmd(cmd)
    if(status == 0): #{
      cmd = ('WlzExtFFConvert -f wlz -F nii ' +
             '-o ' + idx_obj_base + '.nii ' + idx_obj_base + '.wlz')
      status,dummy = runCmd(cmd)
    #}
    if(bool(status)): #{
      errMsg('Failed to create index object ' + idx_obj_base + '.(nii)|(wlz)')
    #}
  #}
  # Find roc for each point (closest to centre of the region of interest and
  # within the region of interest image.
  if(status == 0): #{
    for pnt in model_pts['points']: #{
      vrbMsg(7, 'mapMeasurements() model_pts[\'points\'][' + str(pnt) + '] = ' +
          str(model_pts['points'][pnt]))
      cen_match = [None, None]
      for cen in model_roc: #{
        d0 = dist(model_pts['points'][pnt], model_roc[cen])
        if((d0 >= 0.0) and (d0 < (roi_obj_size / 4.0))): #{
          if((cen_match[1] is None) or (d0 < cen_match[1])): #{
            vrbMsg(9, 'mapMeasurements() cen = ' + str(cen) + 
                ', cen_match = ' + str(cen_match))
            cen_match = [cen, d0]
          #}
        #}
      #}
      vrbMsg(7, 'mapMeasurements() min dist pnt = ' + str(pnt) + 
          ', cen_match = ' + str(cen_match))
      roc_list.append([pnt, cen_match[0]])
    #}
    vrbMsg(7, 'mapMeasurements() roc_list = ' + str(roc_list))
  #}
  # For each assay map the points
  if(status == 0): #{
    rof = roi_obj_size / 2
    for assay in assaylist: #{
      assay_pts = {'assay': {}, 'pixelsz': [], 'points': {}}
      assay_pts['pixelsz'] = getPixelSz(assay, scan_dir, in_fmt)
      assay_pts['assay'] = assay
      # Read assay centres file
      assay_roc_file = roiCentresFile(assay)
      status,assay_roc = readJsnFile(assay_roc_file)
      if(status != 0): #{
        errMsg('Failed to read assay region centres from file: ' +
            assay_roc_file)
        break
      #}
      for idx in range(0, len(roc_list)): #{
        pnt = roc_list[idx][0]
        roi = roc_list[idx][1]
        vrbMsg(7, 'mapMeasurements() pnt = ' + str(pnt) +
            ', roi = ' + str(roi))
        if(not ((pnt is None) or (roi is None) or
                (model_roc[roi] is None) or (assay_roc[roi] is None))): #{
          vrbMsg(7, 'mapMeasurements() model_roc[roi] = ' +
              str(model_roc[roi]) + ', assay_roc[roi] = ' + str(assay_roc[roi]))
          trx_base = trxFileBase(model, assay, roi)
          roi_off_mod = [model_roc[roi][0] - rof, model_roc[roi][1] - rof]
          roi_off_asy = [assay_roc[roi][0] - rof, assay_roc[roi][1] - rof]
          idx_val = None
          asy_pos = None
          mod_pos = model_pts['points'][pnt]
          asy_pos_rel = None
          mod_pos_rel = [mod_pos[0] - roi_off_mod[0],
                         mod_pos[1] - roi_off_mod[1]]
          vrbMsg(7, 'mapMeasurements() mod_pos = ' + str(mod_pos) + 
              ', mod_pos_rel = ' + str(mod_pos_rel))
          # Transform the index object
          cmd = ('antsApplyTransforms -d 2 -n MultiLabel'
                 ' -o ' + tmp_reg_base + '.nii' +
                 ' -i ' + idx_obj_base + '.nii' +
                 ' -r ' + idx_obj_base + '.nii' +
                 ' -t ' + trx_base + 'Affine.txt' +
                 ' -t ' + trx_base + 'Warp.nii.gz')
          status,dummy = runCmd(cmd)
          if(status == 0): #{
            cmd = ('WlzExtFFConvert -f nii -F wlz ' + tmp_reg_base + '.nii |' +
                   'WlzConvertPix -t 1 > ' + tmp_reg_base + '.wlz')
            status,dummy = runCmd(cmd)
          #}
          if(status != 0): #{
            errMsg('Failed to register index object for ' + trx_base)
            break
          #}
          # Find location of the model point in the assay image. Do this by
          #  1. Get model point location
          #  2. Find index value in the warped index object at model location
          #  3. Find location in un-warped index object with this value.
          # We do this rather than just use antsApplyTransformsToPoints
          # because it just doesn't work when you have multiple transforms!
          #
          # 1. Find index value
          idx_val = 0
          cmd = ('WlzGreyValue -x ' + str(mod_pos_rel[0]) +
                 ' -y ' + str(mod_pos_rel[1]) + ' ' + tmp_reg_base + '.wlz')
          status,rtnstr = runCmd(cmd)
          if((status != 0) or (len(rtnstr) < 1)): #{
            status = 1
          else: #}{
            try: #{
              idx_val = int(rtnstr)
            except: #}{
              status = 1
            #}
          #}
          # 2. Find location
          if(status == 0): #{
            cmd = ('WlzThreshold -E -v ' + 
                str(idx_val) + ' ' + idx_obj_base + '.wlz | WlzCentreOfMass')
            status,rtnstr = runCmd(cmd)
          #}
          if(status == 0): #{
            cen = rtnstr.split(' ')
            if(len(cen) != 4): #{
              status = 1
            else: #}{
              asy_pos_rel = [int(float(cen[1])), int(float(cen[2]))]
            #}
          #}
          if(status != 0): #{
            errMsg('Failed to find assay point position for assay = ' +
                assay + ', pnt = ' + pnt + ', roi = ' + roi +
                ', mod_pos = ' + str(mod_pos))
            break
          #}
          asy_pos = [asy_pos_rel[0] + roi_off_asy[0],
                     asy_pos_rel[1] + roi_off_asy[1],
                     0.0, model]
          vrbMsg(5, 'mapMeasurements() assay = ' + assay + ', pnt = ' + pnt +
              ', roi = ' + roi + ' mod_pos = ' + str(mod_pos) +
              ', asy_pos = ' + str(asy_pos))
          if(comp_cnf and (status == 0)): #{
            mean_val = [0.0, 0.0]
            sigma_val = [0.0, 0.0]
            cnf_roi_obj = [None, None]
            for i in range(0, 2): #{
              base = model if(i == 0) else assay
              f = bfiFileBase(base) + '.wlz'
              err_num, obj0 = readWoolzObj(f)
              if(bool(err_num)): #{
                status = 1
              else: #}{
                pos = mod_pos if(i == 0) else asy_pos
                status, cnf_roi_obj[i] = extractROI(obj0, roi, pos, 64,
                    filter_sobel=True, filter_window=2)
                if(status == 0): #{
                  mn = c.c_double(0)
                  sg = c.c_double(0)
                  dummy = w.WlzGreyStats(cnf_roi_obj[i], None, None, None,
                    None, None, c.byref(mn), c.byref(sg), c.byref(err_num))
                  if(bool(err_num)): #{
                    status = 1
                  else: #}{
                    mean_val[i] = mn.value
                    sigma_val[i] = sg.value
                  #}
                #}
              #}
              if(bool(obj0)): #{
                dummy = w.WlzFreeObj(obj0)
                obj0 = None
              #}
              if(status != 0): #{
                errMsg('mapMeasurements() ' +
                    'Failed to read Woolz object from file ' + f)
                status = 1
                break
              #}
            #}
            # Compute confidence value using model (m) and registered
            # assay (a) ROI images extracted using Sobel, Gaussian
            # and Parzen filters.
            # Compute sum image (m + a)
            # Find area (n) of union of domains with values >= 1 for (m + a).
            # Compute mean sum value (sum)
            # Compute difference image (m - a).
            # Compute rms of differences image (dif)
            # Create confidence value cnf = (sum - dif) / sum
            if(status == 0): #{
              area = 0
              nrm_sum = 0.0
              nrm_dif = 0.0
              sumObj = None
              difObj = None
              err_num = c.c_int(w.WLZ_ERR_NONE)
              sumObj = w.WlzAssignObject(
                  w.WlzImageArithmetic(cnf_roi_obj[0], cnf_roi_obj[1],
                      w.WLZ_BO_ADD, 0, c.byref(err_num)), None)
              if(not bool(err_num)): #{
                tmpObj = w.WlzAssignObject(
                    w.WlzThresholdI(sumObj, 1, w.WLZ_THRESH_HIGH,
                    c.byref(err_num)), None)
                w.WlzFreeObj(sumObj);
                sumObj = tmpObj
              #}
              if(not bool(err_num)): #{
                area = w.WlzArea(sumObj, c.byref(err_num))
              #}
              if((not bool(err_num)) and (area > 0)): #{
                tcd = c.c_double(0.0)
                err_num = c.c_int(w.WLZ_ERR_NONE)
                dummy = w.WlzGreyStats(sumObj, None, None, None,
                    c.byref(tcd), None, None, None, c.byref(err_num))
                if(not bool(err_num)): #{
                  nrm_sum = tcd.value / float(1 + area)
                #}
              #}
              w.WlzFreeObj(sumObj);
              if(not bool(err_num)): #{
                difObj = w.WlzAssignObject(
                    w.WlzImageArithmetic(cnf_roi_obj[0], cnf_roi_obj[1],
                        w.WLZ_BO_SUBTRACT, 0, c.byref(err_num)), None)
              #}
              if(debugValue('confidence') > 0): #{
                writeWoolzObj(debugFile('cnf_0') + '.wlz', cnf_roi_obj[0])
                writeWoolzObj(debugFile('cnf_1') + '.wlz', cnf_roi_obj[1])
                writeWoolzObj(debugFile('cnf_c') + '.wlz', difObj)
                if(debugValue('exit') > 0): #{
                  exit(0)
                #}
              #}
              if(not bool(err_num)): #{
                tcd = c.c_double(0.0)
                err_num = c.c_int(w.WLZ_ERR_NONE)
                # Get sum not sum of squares because image values squared
                dummy = w.WlzGreyStats(difObj, None, None, None, None,
                    c.byref(tcd), None, None, c.byref(err_num))
                if(not bool(err_num)): #{
                  nrm_dif = m.sqrt(tcd.value / float(1 + area))
                #}
                vrbMsg(3, 'mapMeasurements() area = ' + str(area) + 
                    ', nrm_sum = ' + str(nrm_sum) +
                    ', nrm_dif = ' + str(nrm_dif))
              #}
              w.WlzFreeObj(difObj);
              if(bool(err_num)): #{
                cnf_val = 0.0
                status = 1
              else: #}{
                eps = 0.0001
                if((mean_val[0] < eps) or (mean_val[1] < eps) or
                   (sigma_val[0] < eps) or (sigma_val[1] < eps) or
                   (area == 0)): #{
                  cnf_val = 0.0
                else: #}{
                  cnf_val = (nrm_sum - nrm_dif) / nrm_sum
                  vrbMsg(3, 'mapMeasurements() cnf_val = ' + str(cnf_val))
                #}
              #}
              asy_pos[2] = cnf_val
            #}
            for i in range(0, 2): #{
              if(bool(cnf_roi_obj[i])): #{
                dummy = w.WlzFreeObj(cnf_roi_obj[i])
                cnf_roi_obj[i] = None
              #}
            #}
          #}
          assay_pts['points'].update({pnt: asy_pos})
        #}
      #}
      if(status == 0): #{
        meas.append(assay_pts)
      #}
    #}
  #}
  if(status != 0): #{
    meas = []
  #}
  vrbMsg(1, 'mapMeasurements() meas = ' + str(meas) + 
      ', status = ' + str(status))
  return meas,status 
#}

def saveMeasurements(filename, measurements): #{
  vrbMsg(1, 'saveMeasurements() filename = ' + filename +
      ', measurements = ' + str(measurements))
  status = 0
  # Open file
  f = None
  try: #{
    if(filename == '-'): #{
      f = sys.stdout
    else: #}{
      f = open(filename, 'w')
    #}
  except: #}{
    status = 1
  #}
  # Save measurements either as jsn or csv (default is csv)
  if(status == 0): #{
    s0 = filename.split('.')
    if((len(s0) > 1) and (s0[-1] == 'jsn')): #{
      # Save as jsn
      vrbMsg(5, 'saveMeasurements() saving to jsn file.')
      json.dump(measurements, f)
    else: #}{
      vrbMsg(5, 'saveMeasurements() saving to csv file.')
      # Collect points
      points = []
      for meas in measurements: #{             
        for p in meas['points']: #{                 
          if(not p in points): #{
            points.append(p)
          #}
        #}
      #}
      # Save as csv 
      s = 'assay,pixelsz_x,pixelsz_y'
      for p in points: #{
        sp = str(p)
        s = s + ',' + sp + '_x,' + sp + '_y,' + sp + '_c,' + sp + '_m'
      #}
      print(s, file=f)
      for meas in measurements: #{
        z = meas['pixelsz']
        s = str(meas['assay']) + ',' + str(z[0]) + ',' + str(z[1])
        for idx in range(0, len(points)): #{
          p = points[idx]
          if(p in meas['points']): #{
            q = meas['points'][p]
            s = (s + ',' + str(q[0]) + ',' + str(q[1]) + ',' +
                 str(q[2]) + ',' + str(q[3]))
          else: #}{
            s = s + ',0,0,0.0,0'
          #}
        #}
        print(s, file=f)
      #}
    #}
  #}
  # Close file
  if((not f is None) and (filename != '-')): #{
    f.close()
  #}
  vrbMsg(1, 'saveMeasurements() status = ' + str(status))
  return status
#}

def plotMeasurements(meas, file_base): #{
  import matplotlib.pyplot as plt
  vrbMsg(1, 'plotMeasurements() meas = ' + str(meas))
  status = 0
  UPV = c.POINTER(c.c_void_p)
  UPP = c.POINTER(c.POINTER(c.c_ubyte))
  img = None
  sz = w.WlzIVertex2()
  org = w.WlzIVertex2()
  bfi_file = bfiFileBase(meas['assay']) + '.wlz'
  err_num, bfi_obj = readWoolzObj(bfi_file)
  if(not bool(err_num)): #{
    aryc = c.cast(0,UPV)
    err_num = c.c_int(w.WLZ_ERR_NONE)
    bbox = w.WlzBoundingBox2I(bfi_obj, c.byref(err_num))
    if(not bool(err_num)): #{
      org.vtX = bbox.xMin
      org.vtY = bbox.yMin
      sz.vtX = bbox.xMax - bbox.xMin + 1
      sz.vtY = bbox.yMax - bbox.yMin + 1
      err_num = w.WlzToArray2D(c.byref(aryc), bfi_obj, sz, org, 0,
                               c.c_int(w.WLZ_GREY_UBYTE))
    #}
    if(not bool(err_num)): #{
      aryc = c.cast(aryc, UPP)
      img = np.ctypeslib.as_array(aryc.contents, (sz.vtY, sz.vtX))
    #}
  #}
  if(bool(err_num)): #{
    status = 1
  #}
  if(status == 0): #{
    fig, ax = plt.subplots(1, 1, figsize=(3, 9))
    ax.axis('off')
    ax.set_title('assay = ' + str(meas['assay']))
    ax.imshow(img, cmap='gray', interpolation='nearest')
    points = meas['points']
    for pnt in points: #{
      p = points[pnt]
      # Ignore points with confience value 0.0
      if(p[2] > 0.0): #{
        p = [p[0] - org.vtX, p[1] - org.vtY]
        ax.plot(p[0], p[1], 'x', mec='g')
      #}
    #}
    if((file_base is None) or (len(file_base) == 0)): #{
      plt.show()
    else: #}{
      f = file_base + '_' + str(meas['assay']) + '.png'
      plt.savefig(f)
    #}
  #}
  vrbMsg(1, 'plotMeasurements() status = ' + str(status))
  return status
#}

#======================================================================#

def ParseArgs(): #{
  status = 0
  parser = argparse.ArgumentParser(add_help=False, description= 
  '''Computes, combines and plots long bone measurements from DXA scan images.
  Measurements made on a model image may be transformed to an assay image
  through region of interest extraction and non-linear registration. These
  measurements may then also be plotted and combined.''')
  parser.add_argument('-f', '--file-format', 
      type=str, default=def_format,
      help='Input model/assay file format, eg: dcm, jpg, tif, ....')
  parser.add_argument('-d', '--debug', 
      type=int, default=0,
      help='For debuging only. Setting this may give incorrect measurements' +
           ', cause the program to exit, write large files, ....')
  parser.add_argument('-v', '--verbose', 
      type=int, default=0,
      help='Verbose output (0 = none (default), 9 = max' +
           ' (probably only useful for debugging)).')
  parser.add_argument('-E', '--everything',
      action='store_true', default=False,
      help='Do (almost) everything: Process model and assay images, ' +
           'register assays, map measurements, save measurements.')
  parser.add_argument('-M', '--process-model',
      action='store_true', default=False,
      help='Process model image.')
  parser.add_argument('-A', '--process-assay',
      action='store_true', default=False,
      help='Process assay images.')
  parser.add_argument('-C', '--combine-measurements',
      action='store_true', default=False,
      help='Combine the input measurements by selecting highest confidence ' +
           'measurement values.')
  parser.add_argument('-P', '--map-measurements',
      action='store_true', default=False,
      help='Map measurements from model to assay.')
  parser.add_argument('-L', '--plot-measurements',
      action='store_true', default=False,
      help='Plot point possitions over image.')
  parser.add_argument('-R', '--register-assays',
      action='store_true', default=False,
      help='Register assays.')
  parser.add_argument('-c', '--force-to-model',
      action='store_true', default=False,
      help = 'When reading input measurements force them to be for the ' +
             'given model.')
  parser.add_argument('-i', '--infile', 
      type=str, default='', 
      help='Input file (or list of files) previously computed measurements.')
  parser.add_argument('-l', '--plot-file',
      type=str, default='', 
      help='Filebase for plots, default is plots not saved but shown ' +
          ' using matplotlib.pyplot.show().')
  parser.add_argument('-a', '--colour-channel',
      type=str, default='modulus',
      help='Colour channel to use if input model/assay is colour. ' +
          'Options are: ' + str(wlz_rgba_channel.keys()) + '.')
  parser.add_argument('-b', '--background-valid',
      action='store_true', default=False,
      help='Assays have valid background.')
  parser.add_argument('-m', '--model', 
      type=str, default=def_model, 
      help='Target model for the assays.')
  parser.add_argument('-o', '--outfile', 
      type=str, default='', 
      help='Output file for measurements. ' +
           'The default is that no measurements are output. ' +
           'Use \'-\' for the standard output.')
  parser.add_argument('-p', '--model-points',
      type=str, default='',
      help='File with known model points (default = pts-<model>.jsn).')
  parser.add_argument('-r', '--threshold-cnf', 
      type=float, default=min_cnf_value, 
      help='Minimum confidence value accepted when combining measurements ' +
           '(range [0.0 - 1.0], default ' + str(min_cnf_value) + ').')
  parser.add_argument('-s', '--scan-dir',
      type=str, default=def_scan_dir,
      help='Input scan directory.')
  parser.add_argument('-t', '--max-threads', 
      type=int, default=0, 
      help='Maximum number of threads to use (0 implies maximum available).')
  parser.add_argument('-w', '--work-dir',
      type=str, default=def_work_dir,
      help='Working files directory.')
  parser.add_argument('-h', '--help',
      action='store_true', default=False,
      help='show this usage message and then exit.')
  parser.add_argument('assaylist',
      default='', nargs='?',
      help='Input assay list.')
  args = parser.parse_args()
  if(args.help): #{
    parser.print_help()
    print('\nVersion: ' + version)
    print(
   '''\nExamples:
  1 - compute measurements using model 1848 for assays 2617, 2674 
  and 2787 putting the output into m-1848.csv with verbosity set to 1
  so as to get some progress output:
    bonemeasurement.py -E -v1 -o m-1848.csv -m 1848 2617,2674,2787
  2 - combine measurements for m-1848.csv, m-2674.csv and m-4775.csv and
  put the resulting measurements into m-cmb.csv
    bonemeasurement.py -C -i m-1848.csv,m-2674.csv,m-4775.csv -o m-cmb.csv
  3 - plot measurements in file m-cmb.csv
    bonemeasurement.py -L -i m-cmb.csv''')
  #}
  if(args.combine_measurements and args.map_measurements): #{
    errMsg('Measurements can either be mapped of combined but not both.')
    status = 1
  #}
  if(not (args.colour_channel in wlz_rgba_channel)): #{
    errMsg('Colour channel invalid, valid channels are: ' +
        str(wlz_rgba_channel.keys()) + '.')
    status = 1
  #}
  assaylist = []
  global debug_mask
  debug_mask = args.debug
  if(len(args.assaylist) > 0): #{
    assaylist = args.assaylist.split(',')
  #}
  if(len(args.model_points) == 0): #{
    args.model_points = 'pts-' + args.model + '.jsn'
  #}
  if(len(args.outfile) == 0): #{
    args.outfile = None
  #}
  if(args.threshold_cnf < 0.0): #{
    args.threshold_cnf = 0.0
  elif(args.threshold_cnf > 1.0): #}{
    args.threshold_cnf = 1.0
  #}
  return status,args,assaylist
#}

def Main(): #{
  f = None
  status = 0
  infile_list = []
  measurements = []
  vrbMsg(3, 'model = ' + str(args.model))
  if(args.max_threads > 0): #{
    os.environ['OMP_NUM_THREADS'] = str(args.max_threads)
  #}
  if(args.everything or args.process_model): #{
    status = processImage(args.model, args.model, args.file_format)
  #}
  if(status == 0): #{
    for assay in assaylist: #{
      if(args.everything or args.process_assay): #{
        status = processImage(args.model, assay, args.file_format)
        if(status != 0): #{
          break
        #}
      #}
      if(args.everything or args.register_assays): #{
        status = registerAssay(args.model, assay)
        if(status != 0): #{
          break
        #}
      #}
    #}
  #}
  if((status == 0) and (len(args.infile) != 0)): #{
    infile_list = args.infile.split(',')
  #}
  if((status == 0) and (len(infile_list) > 0)): #{
      mod = None
      if(args.force_to_model): #{
        mod = args.model
      #}
      for infile in infile_list: #{
        meas,status = readMeasurements(infile, mod)
        if(status == 0): #{
          measurements.extend(meas)
        else: #}{
          break
        #}
      #}
    #}
  if((status == 0) and (args.everything or args.map_measurements)): #{
    meas, status = mapMeasurements(args.model, args.model_points, assaylist,
        args.scan_dir, args.file_format)
    if(status == 0): #{
      measurements.extend(meas)
      model_meas, status = modelMeasurements(args.model, args.model_points,
          args.scan_dir, args.file_format)
      if(status == 0): #{
        measurements.append(model_meas)
      #}
    #}
  #}
  if((status == 0) and args.combine_measurements): #{
    measurements = combineMeasurements(measurements, args.threshold_cnf)
  #}
  if((status == 0) and bool(args.outfile)): #{
    status = saveMeasurements(args.outfile, measurements)
  #}
  if(status == 0): #{
    if(args.plot_measurements): #{
      for meas in measurements: #{
        status = plotMeasurements(meas, args.plot_file)
      #}
    #}
  #}
  return status 
#}

if __name__ == '__main__': #{
  status = 0
  a = 'a-zA-Z0-9_'
  prog = re.compile('(^[^' + a + ']+)|(\.[' + a + ']+$)').sub('', sys.argv[0])
  status,args,assaylist = ParseArgs()
  vrb_level = args.verbose
  vrbMsg(3, 'args = ' + str(args))
  vrbMsg(3, 'assaylist = ' + str(assaylist))
  if(status == 0): #{
    status = Main()
  #}
  exit(status)
#}
