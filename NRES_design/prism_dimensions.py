import numpy as np
import os, sys, time

long_side_full = 193.5      # including beveled corners
long_side_flat = 192.67     # width of usable flat face
apex_angle_deg =  55.0

apex_length_mm = 174.0      # length of bisector from short side to apex

short_side_flat = 177.32    # very this!

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

