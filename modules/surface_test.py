import os, sys, time
import numpy as np

def vec_length(vector):
    return np.sqrt(np.sum(vector**2))

def make_unitvec(x, y, z):
    temp = np.array([x, y, z]).astype('float')
    return temp / vec_length(temp)

def angle_sep_rad(vec1, vec2):
    return np.arccos(np.dot(vec1, vec2))

def angle_sep_deg(vec1, vec2):
    return np.degrees(angle_sep_rad(vec1, vec2))

def calc_incident_angle(v_incident, surf_norm):
    return angle_sep_rad(-v_incident, surf_norm)

def calc_tangent_vector(v_incident, surf_norm):
    temp = v_incident - surf_norm * np.dot(v_incident, surf_norm)
    return temp / vec_length(temp)  # as unit vector

def calc_refract_angle(incid_rad, n1n2_ratio):
    return np.arcsin(n1n2_ratio * np.sin(incid_rad))

def calc_refract_vector(v_incident, surf_norm, n1n2_ratio):
    theta_1 = calc_incident_angle(v_incident, surf_norm)
    theta_2 = calc_refract_angle(theta_1, n1n2_ratio)
    tgt_vec = calc_tangent_vector(v_incident, surf_norm)
    refrvec = tgt_vec * np.sin(theta_2) - surf_norm * np.cos(theta_2)
    return refrvec

def test_tgt_vec(v_incident, surf_norm):
    temp = v_incident + surf_norm * np.dot(v_incident, -surf_norm)
    return temp

## Calculate vectors of reflection and reflaction (Wikipedia):
#def wiki_calc_refract_vector(v_incident, surf_norm, n1_n2_ratio):
def calc_surface_vectors(v_incident, surf_norm, n1_n2_ratio):
    cti = -1.0 * np.dot(v_incident, surf_norm)     # cos(theta_i)
    nnr = n1_n2_ratio
    #sys.stderr.write("cti: %10.5f\n" % cti)
    v_reflect = v_incident + 2. * cti * surf_norm
    #sys.stderr.write("v_incident: %s\n" % str(v_incident))
    #sys.stderr.write("v_reflect:  %s\n" % str(v_reflect))
    smult = nnr*cti - np.sqrt(1.0 - nnr**2 * (1.0 - cti**2))
    v_refract = nnr * v_incident + smult * surf_norm
    return v_reflect, v_refract
    #return v_refract



vecsize = 1.0

surface_outer = np.array([0.0, 0.0, 1.0])
surface_inner = -1.0 * surface_outer

almost_normal = make_unitvec(0.1, 0, 1)
almost_faceon = -almost_normal

#fortyfive_xz = np.array([1.0, 0.0, 1.0])
fortyfive_xz = make_unitvec(1., 0., 1.)
fortyfive_yz = make_unitvec(0., 1., 1.)

#inbound_45 = -fortyfive_xz
incident_vec = almost_faceon
test_result = calc_refract_vector(incident_vec, surface_outer, 0.9)

v_reflect, v_refract = \
    calc_surface_vectors(incident_vec, surface_outer, 0.9)

