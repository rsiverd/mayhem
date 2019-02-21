import numpy as np
import os, sys, time

## OLD DELIVERY DOC GIVES:
## 190.0 x 184.2 x 130.0

## NEWER DELIVERY DOC GIVES:
## 190.0 x 205.7 x 13 x 55deg

long_side_full = 193.5      # including beveled corners
long_side_flat = 192.67     # width of usable flat face
apex_angle_deg =  55.0

apex_length_mm = 174.0      # length of bisector from short side to apex

short_side_flat = 177.32    # very this!

## Side length ratio for isosceles traingle with 55deg apex:
ls_ratio = 1.0 / (2.0 * np.sin(np.radians(0.5 * apex_angle_deg)))
sys.stderr.write("ls_ratio = %10.6f\n" % ls_ratio)

## Larger angles are equal (isosceles triangle):
large_angle_deg = 0.5 * (180.0 - apex_angle_deg)

## Law of sines provides short edge:
apex_angle_rad  = np.radians(apex_angle_deg)
large_angle_rad = np.radians(large_angle_deg)
short_side_calc = np.sin(apex_angle_rad) \
                * long_side_flat / np.sin(large_angle_rad)

calc_apex_length = np.sin(np.radians(large_angle_deg)) * long_side_full

## Report:
sys.stderr.write("short_side_flat (from drawings): %10.5f\n" % short_side_flat)
sys.stderr.write("short_side_calc (from geometry): %10.5f\n" % short_side_calc)

