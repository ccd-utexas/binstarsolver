#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utilities for calculating physical quantities from observed quantities.

"""


# Import standard packages.
from __future__ import absolute_import, division, print_function
import sys
import warnings
# Import installed packages.
import numpy as np
import scipy.optimize as sci_opt
import scipy.constants as sci_con
import matplotlib.pyplot as plt


def calc_flux_intg_ratio_from_mags(mag_1, mag_2):
    """Calculate the ratio of integrated fluxes from two magnitudes.
    
    Parameters
    ----------
    mag_1 : float
        Magnitude of source 1. Unit is magnitudes.
    mag_2 : float
        Magnitude of source 2. Unit is magnitudes.
        
    Returns
    -------
    flux_intg_ratio : float
        Ratio of fluxes for sources 1 and 2 as flux_1 / flux_2. Unitless.
    
    Notes
    -----
    flux_1 / flux_2 = 100**((mag_2 - mag_1)/5)
    From equation 3.3 in section 3.2 of [1]_.
    
    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics
    
    """
    flux_intg_ratio = 100.0**((mag_2 - mag_1)/5.0)
    return flux_intg_ratio


def calc_fluxes_intg_rel_from_light(light_oc, light_ref=1.0):
    """Calculate integrated fluxes of the of both binary stars relative to the total integrated flux.
    
    Parameters
    ----------
    light_oc : float
        Light level during occultation (primary eclipse). Light level is the average of the flat region at the bottom
        of one of the eclipse minima. The event is a total eclipse of the star with the smaller radius.
    light_ref : {1.0}, float, optional
        Reference for light levels outside of eclipse. Typically total system light is normalized to 1.0.
        
    Returns
    -------
    flux_intg_rel_s : float
        Integrated flux of greater-sized star as a fraction of total integrated flux.
    flux_intg_rel_g : float
        Integrated flux of greater-sized star as a fraction of total integrated flux.
    
    Notes
    -----
    Note: Fluxes are returned as a tuple: (fs, fg)
    flux_intg_rel_s = (light_ref - light_oc) / light_ref
    flux_intg_rel_g = light_oc / light_ref
    From equation 7.2 in section 7.3 of [1]_.
    
    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry
    
    """
    flux_intg_rel_s = (light_ref - light_oc) / light_ref
    flux_intg_rel_g = light_oc / light_ref
    return (flux_intg_rel_s, flux_intg_rel_g)


def calc_phase_orb_from_time_period(time_event, period, time_mideclipse=0.0):
    """Calculate orbital phase angle in radians for an event relative
    to mid-eclipse, given period.
    
    Parameters
    ----------
    time_event : float
        Time of event. Units are seconds.
    period : float
        Period of eclipse. Units are seconds.
    time_mideclipse : {0.0}, float, optional
        Mid-eclipse time. Units are seconds. Typically `time_event`
        is relative to `time_mideclipse` such that `time_mideclipse` == 0.0.
        
    Returns
    -------
    phase_orb : float
        Orbital phase angle in radians for `time_event`.
    
    Notes
    ------
    phase_orb = 2*pi*(t - t0) / period
    where t = `time_event`, t0 = `time_mideclipse`
    From equation 7.3 in section 7.3 of [1]_.

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry
    
    """
    phase_orb = 2.0 * np.pi * (time_event - time_mideclipse) / period
    return phase_orb


def calc_sep_proj_from_incl_phase(incl, phase_orb):
    """Calculate the projected separation of the two star centers in units of the
    star-star separation distance. Assumes orbit is circular about the system center of mass.
    
    Parameters
    ----------
    incl : float
        Orbital inclination. Angle between line of sight and the axis of the orbit. Unit is radians.
    phase_orb : float
        Orbital phase angle of event in radians.
    
    Returns
    -------
    sep_proj : float
        Separation of the two star centers projected onto the tangent plane of the sky.
        Unit is orbit radius.
    
    Notes
    -----
    sep_proj^2 = cos^2(incl) + sin^2(phase_orb)*sin^2(incl)
    From equation 7.3 in section 7.3 of [1]_.

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry
    
    """
    sep_proj = np.sqrt((np.cos(incl)**2.0) + (np.sin(phase_orb)**2.0)*(np.sin(incl)**2.0))
    return sep_proj


def calc_radii_ratio_from_light(light_oc, light_tr, light_ref=1.0):
    """Calculate ratio of radii of smaller-sized star to greater-sized star from light levels 
    during occultation and transit. Assumes L \propto r^2, no limb darkening.
    
    Method may not be valid for radii ratios > 10. See "Notes".
    
    Parameters
    ----------
    light_oc : float
        Light level during occultation (primary eclipse). Light level is the average
        of the flat region at the bottom of one of the eclipse minima. The event is
        a total eclipse of the star with the smaller radius.
    light_tr : float
        Light level during transit (secondary eclipse). Light level is the average
        of the rounded region at the bottom of one of the eclipse minima. The event is
        an annular eclipse of the star with the larger radius.
    light_ref : {1.0}, float, optional
        Reference for light levels outside of eclipse. Typically total system light is normalized to 1.0.
        
    Returns
    -------
    radii_ratio_lt : float
        Ratio of radii of smaller-sized star to greater-sized star calculated from light levels
        during occultation and transit.
        radii_ratio = radius_s / radius_g

    Notes
    -----
    Note: Method may not be valid for stars in different stages of evolution or for radii ratios > 10
    (e.g. a binary system with main sequence star and a red giant)
    radii_ratio = radius_s / radius_g
                = sqrt((light_ref - light_tr) / light_oc)
    From equation 7.8 in section 7.3 of [1]_.

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry
    
    """
    radii_ratio_lt = np.sqrt((light_ref - light_tr) / light_oc)
    return radii_ratio_lt


def calc_radii_sep_from_seps(sep_proj_ext, sep_proj_int):
    """Calculate the radii of both binary stars from projections of star separations
    at external and internal tangencies. Method is independent of any inclination assumptions.
    
    Parameters
    ----------
    sep_proj_ext : float
        Projected separation of star centers at external tangencies.
        (e.g. begin ingress, end egress). Unit is star-star separation distance.
    sep_proj_int : float
        Projected separation of star centers at internal tangencies.
        (e.g. end ingress, begin egress). Unit is star-star separation distance.
    
    Returns
    -------
    radius_sep_s : float
        Radius of smaller-sized star. Unit is star-star separation distance.
    radius_sep_g : float
        Radius of greater-sized star. Unit is star-star separation distance.
    
    See Also
    --------
    calc_radii_ratio_from_light, calc_radius_from_velrs_times 
    
    Notes
    -----
    Note: Radii are returned as a tuple: (rs, rg)
    Note: Method does not assume an inclination.
    radii_ratio = radius_sep_s / radius_sep_g
    sep_proj_ext = radius_sep_g * (1 + radii_ratio)
    sep_proj_int = radius_sep_g * (1 - radii_ratio)
    => sep_proj_ext - sep_proj_int = 2 * radius_sep_g * radii_ratio = 2 * radius_sep_s
       radius_sep_s = (sep_proj_ext - sep_proj_int) / 2
    => sep_proj_ext + sep_proj_int = 2 * radius_sep_g
       radius_sep_g = (sep_proj_ext + sep_proj_int) / 2
    From equations 7.8, 7.9, 7.10 in section 7.3 of [1]_.
    
    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry
    
    """
    radius_sep_s = (sep_proj_ext - sep_proj_int) / 2.0    
    radius_sep_g = (sep_proj_ext + sep_proj_int) / 2.0
    return (radius_sep_s, radius_sep_g) 


def calc_radii_ratio_from_rads(radius_sep_s, radius_sep_g):
    """Calculate ratio of radii of smaller-sized star to greater-sized star.
    
    Parameters
    ----------
    radius_sep_s : float
        Radius of smaller-sized star. Unit is star-star separation distance.
    radius_sep_g : float
        Radius of greater-sized star. Unit is star-star separation distance.
    
    Returns
    -------
    radii_ratio_rad : float
        Ratio of radii of smaller-sized star to greater-sized star calculated from star radii.
        radii_ratio = radius_s / radius_g
    
    Notes
    -----
    radii_ratio = radius_sep_s / radius_sep_g
    From equations 7.8 in section 7.3 of [1]_.
    
    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry
    
    """
    radii_ratio_rad = radius_sep_s / radius_sep_g
    return radii_ratio_rad


def calc_incl_from_radii_ratios_phase_incl(
    radii_ratio_lt, phase_orb_ext, phase_orb_int, tol=1e-4, maxiter=10, show_plots=False):
    """Calculate inclination angle by minimizing difference of ratios of stellar radii as calulated
    from light levels and light curve events (tangencies) for various values of inclination.
    
    Parameters
    ----------
    radii_ratio_lt : float
        Ratio of radii of smaller-sized star to greater-sized star
        calculated from light levels during occultation and transit.
        radii_ratio = radius_s / radius_g
    phase_orb_ext : float
        Orbital phase angle at external tangencies (e.g. begin ingress, end egress). Unit is radians.
    phase_orb_int : float
        Orbital phase angle at internal tangencies (e.g. end ingress, begin egress). Unit is radians.
    tol : {1e-4}, float, optional
        Maximum tolerance for difference in radii ratios at a self-consistent solution for inclination,
        i.e. at solved inclination:
        abs('radii ratio from light levels' - 'radii ratio from eclipse events') < tol.
    maxiter : {10}, int, optional
        Maximum number of iterations to perform when solving for inclination.
        For `tol = 1e-4`, the solution typically converges within 5 iterations.
    show_plots : {True}, bool, optional
        Create and show diagnostic plots of difference in independent radii ratio values vs inclination angle.
        Use to check solution in case initial guess for inclination angle caused convergence to wrong solution
        for inclination angle.
    
    Returns
    -------
    incl : float
        Orbital inclination. Angle between line of sight and the axis of the orbit. Unit is radians.
        If solution does not exist or does not converge, `incl = numpy.nan`
    
    See Also
    --------
    calc_sep_proj_from_incl_phase, calc_radius_sep_s_from_sep, calc_radius_sep_g_from_sep, calc_radii_ratio_from_light
    
    Notes
    -----
    - Method uses calc_radii_ratio_from_light, which may not be valid for stars in different stages of evolution
        or for radii ratios > 10 (e.g. a binary system with main sequence star and a red giant)
    - Compares two independent ways of calculating radii_ratio in order to infer inclination angle. 
    - Equations from section 7.3 of [1]_.
    
    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    # Make radii_ratio_rad and all dependencies functions of inclination.
    # Note: stdout and stderr are delayed if called within a loop in an IPython Notebook.
    sep_proj_ext = lambda incl: calc_sep_proj_from_incl_phase(incl=incl, phase_orb=phase_orb_ext)
    sep_proj_int = lambda incl: calc_sep_proj_from_incl_phase(incl=incl, phase_orb=phase_orb_int)
    radii_sep = lambda incl: calc_radii_sep_from_seps(
        sep_proj_ext=sep_proj_ext(incl=incl), sep_proj_int=sep_proj_int(incl=incl))
    radii_ratio_rad = lambda incl: calc_radii_ratio_from_rads(*radii_sep(incl=incl))
    diff_radii_ratios = lambda incl: radii_ratio_lt - radii_ratio_rad(incl=incl)
    fmt_parameters =  \
      ("    radii_ratio_lt = {rrl}\n" +
       "    phase_orb_ext  = {poe}\n" +
       "    phase_orb_int  = {poi}\n" +
       "    tol            = {tol}").format(
           rrl=radii_ratio_lt, poe=phase_orb_ext, poi=phase_orb_int, tol=tol)
    # Minimize difference between independent radii_ratio values to within a tolerance.
    # Note: A naive implentation of scipy.optimize.minimize will not find the solution
    # for some parameters due to the non-differentiability of `abs(diff_radii_ratios)`
    # Note: diff_radii_ratios is monotonically decreasing with a zero
    # at the solution for self-consistent inclination:
    # - For incl < incl_soln, diff_radii_ratios > 0.
    # - For incl > incl_soln, diff_radii_ratios < 0.
    # Find the solution to within `tol` by recursively zooming in where
    # `diff_radii_ratios` changes sign.
    incls = np.deg2rad(np.linspace(start=0.0, stop=90.0, num=10, endpoint=True))
    diffs = np.asarray(map(diff_radii_ratios, incls))
    if not (np.all(np.diff(diffs) <= 0.0)):
        raise AssertionError(
            "Program error. `diff_radii_ratios` must be monotonically decreasing.")
    if (diffs[0] > 0.0) and (diffs[-1] < 0.0):
        itol = diffs[0] - diffs[-1]
        inum = 0
        while (itol > tol) and (inum < maxiter):
            if inum > 0:
                incls = \
                    np.linspace(
                        start=incls[idx_diff_least_pos], stop=incls[idx_diff_least_neg],
                        num=10, endpoint=True)
                diffs = np.asarray(map(diff_radii_ratios, incls))
            idx_diff_least_pos = len(diffs[diffs > 0.0]) - 1
            idx_diff_least_neg = -1 * len(diffs[diffs < 0.0])
            incl = incls[idx_diff_least_pos]
            itol = abs(diff_radii_ratios(incl=incl))
            inum += 1
        # Check exit condition and solution.
        if (itol > tol) and (inum >= maxiter):
            incl = np.nan
            warnings.warn(
                ("\n" +
                "    Difference in radii ratios did not converge to within tolerance.\n" +
                "    Input parameters:\n" +
                fmt_parameters))
        if radii_ratio_rad(incl=incl) < 0.1:
            warnings.warn(
                ("\n" +
                "    From eclipse timing events, ratio of smaller star's radius\n" +
                "    to greater star's radius is < 0.1. The radii ratio as calculated\n" +
                "    from light levels may not be valid (e.g. for a binary system with\n" +
                "    a main sequence star and a red giant).\n" +
                "    VALID:\n" +
                "    radii_ratio_rad = radius_s / radius_g from eclipse timings = {rtime}\n" +
                "    MAYBE INVALID:\n" +
                "    radii_ratio_lt  = radius_s / radius_g from light levels = {rlt}").format(
                    rtime=radii_ratio_rad(incl=incl),
                    rlt=radii_ratio_lt))
    else:
        incl = np.nan
        warnings.warn(
            ("\n" +
             "    Inclination does not yield self-consistent solution for model.\n" +
             "    Input parameters cannot be fit by model:\n" +
             fmt_parameters))
    # Create and show diagnostic plots.
    if show_plots:
        incls_out = np.deg2rad(np.linspace(start=0, stop=90, num=100))
        diff_radii_ratios_out = map(diff_radii_ratios, incls_out)
        plt.plot(np.rad2deg(incls_out), diff_radii_ratios_out)
        plt.axhline(0.0, color='black', linestyle='--')
        plt.title("Difference between independent radii ratio values\n" +
                  "vs inclination angle (zoomed out view)")
        xlabel = "inclination angle (degrees)"
        ylabel = "radii ratio from light levels - radii ratio from eclipse events"
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()
        if incl is not np.nan:
            incls_in = \
                np.deg2rad(
                    np.linspace(
                        start=np.rad2deg(incl) - 1.0,
                        stop=min(np.rad2deg(incl) + 1.0, 90.0),
                        num=100))
            diff_radii_ratios_in = map(diff_radii_ratios, incls_in)
            plt.plot(np.rad2deg(incls_in), diff_radii_ratios_in)
            plt.axhline(0.0, color='black', linestyle='--')
            plt.title("Difference between independent radii ratio values\n" +
                      "vs inclination angle (zoomed in view)")
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.show()
        else:
            warnings.warn("\nNo inclination solution found.")
    return incl


def calc_semimaj_axis_from_period_velr_incl(period, velr, incl):
    """Calculate semi-major axis of a star's orbit from observed period, radial velocity,
    and orbital inclination. Assumes no orbital eccentricity.
    
    Parameters
    ----------
    period : float
        Period of eclipse. Unit is seconds.
    velr : float
        Observed radial velocity of star. Unit is m/s.
    incl : float
        Orbital inclination. Angle between line of sight and the axis of the orbit. Unit is radians.

    Returns
    -------
    axis : float
        Semi-major axis of star's orbit. Unit is meters.
    
    Notes
    -----
    v = 2*pi*a / P, where v is orbital velocity, a is semi-major axis, P is orbital period.
      = vr / sin(i), where vr is observed radial orbital velocity, i is orbital inclination.
    => a = (P / (2*pi)) * (vr / sin(i))
    From section 7.3 of [1]_.
    
    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics
    
    """
    axis = (period / (2.0*np.pi)) * (velr / np.sin(incl))
    return axis


def calc_sep_from_semimaj_axes(axis_1, axis_2):
    """Calculate separation distance between binary stars from semi-major axes.
    Assumes circular orbits.
    
    Parameters
    ----------
    axis_1 : float
        Semi-major axis of star 1's orbit. Unit is meters.
    axis_2 : float
        Semi-major axis of star 2's orbit. Unit is meters.
    
    Returns
    -------
    sep : float
        Separation distance between stars. Unit is meters.
        
    Notes
    -----
    From section 7.3 of [1]_.
    
    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics
    
    """
    sep = axis_1 + axis_2
    return sep


def calc_radius_from_radius_sep(radius_sep, sep):
    """Calculate the radius of a star converting from units of star-star separation distance
    to meters.
    
    Parameters
    ----------
    radius_sep : float
        Radius of star. Units are star-star separation distance.
    sep : float
        Separation distance between stars. Units are meters.
        
    Returns
    -------
    radius : float
        Radius of star. Units are meters.
    
    Notes
    -----
    radius = radius_sep * sep
    
    """
    radius = radius_sep * sep
    return radius


def calc_radius_from_velrs_times(velr_1, velr_2, time_1, time_2):
    """Calculate the radius of a star from the radial velocities of the stars and\
    relative times of eclipse events. Assumes inclination = 90 deg, no eccentricity. 
    
    Parameters
    ----------
    velr_1 : float
        Observed radial velocity of star 1. Unit is m/s.
    velr_2 : float
        Observed radial velocity of star 2. Unit is m/s.
    time_1 : float
        Relative time of beginning of ingress (first contact). Unit is seconds.
    time_2 : float
        Relative time of end of ingress or beginning of egress.
        during primary minimum (occultation event). Unit is seconds.
        If eclipse is not total, relative time of end of ingress == relative time of beginning of egress.
        For radius of smaller-sized star: time_2 = relative time of end of ingress.
        For radius of greater-sized star: time_2 = relative time of beginning of egress.
    
    Returns
    -------
    radius : float
        Radius of star. Unit is meters.
        For radius of smaller-sized star: time_2 = relative time of end of ingress.
        For radius of greater-sized star: time_2 = relative time of beginning of egress.
        
    See Also
    --------
    calc_radii_ratio_from_light, calc_radius_sep_g_from_sep, calc_radius_sep_s_from_sep
        
    Notes
    -----
    Note: Calculated radius may not agree with radii from methods that do not  
    rs = (vr1 + vr2)/2 * (t2 - t1), t2 is end of ingress
    rg = (vr1 + vr2)/2 * (t2 - t1), t2 is beginning of egress
    From equations 7.8, 7.9 in section 7.3 of [1]_.
    
    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics
    
    """
    radius = ((velr_1 + velr_2) / 2.0) * (time_2 - time_1)
    return radius


def calc_mass_ratio_from_velrs(velr_1, velr_2):
    """Calculate ratio of stellar masses from observed radial velocities.
    Assumes circular orbit.
    
    Parameters
    ----------
    velr_1 : float
        Observed radial velocity of star 1. Unit is m/s.
    velr_2 : float
        Observed radial velocity of star 2. Unit is m/s.
    
    Returns
    -------
    mass_ratio : float
        Ratio of stellar masses for stars 1 and 2 as mass1 / mass2. Unitless.
        
    Notes
    -----
    m1 / m2 = a2 / a1, where a is semi-major axis of low-eccentricty orbit.
    v = 2*pi*a / P, where v is orbital velocity, P is orbital period.
      = vr / sin(i), where vr is observed radial orbital velocity, i is orbital inclination.
    => m1 / m2 = v2r / v1r
    From equation 7.5 in section 7.3 of [1]_.
    
    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    mass_ratio = velr_2 / velr_1
    return mass_ratio


def calc_mass_sum_from_period_velrs_incl(period, velr_1, velr_2, incl):
    """Calculate the sum of stellar masses from observed period, radial velocities,
    and orbital inclination. Assumes orbital excentricity << 1.
    
    Parameters
    ----------
    period : float
        Period of eclipse. Unit is seconds.
    velr_1 : float
        Observed radial velocity of star 1. Unit is m/s.
    v2r : float
        Observed radial velocity of star 2. Unit is m/s.
    incl : float
        Orbital inclination. Angle between line of sight and the axis of the orbit. Unit is radians.

    Returns
    -------
    mass_sum : float
        Sum of stellar masses for stars 1 and 2. Unit is kg.
    
    Notes
    -----
    a = (P/(2*pi)) * (v1 + v2), where a is semi-major axis of low-eccentricity orbit,
        P is orbital period, v is orbital velocity.
    P**2 = ((4*pi**2) / (G*(m1 + m2))) * a**3, where G is gravitational constant,
    m is stellar mass. Kepler's Third Law.
    v = vr / sin(i), where vr is observed radial orbital velocity, i is orbital inclination.
    => m1 + m2 = P/(2*pi*G) * ((v1r + v2r) / sin(i))**3    
    From equation 7.6 in section 7.3 of [1]_.
    
    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics
    
    """
    mass_sum = (period / (2.0*np.pi*sci_con.G)) * ((velr_1 + velr_2) / np.sin(incl))**3.0
    return mass_sum


def calc_masses_from_ratio_sum(mass_ratio, mass_sum):
    """Calculate the individual stellar masses from their sum and ratio.
    
    Parameters
    ----------
    mass_ratio : float
        Ratio of stellar masses for stars 1 and 2 as mass1 / mass2. Unitless.
    mass_sum : float
        Sum of stellar masses for stars 1 and 2. Unit is kg.
    
    Returns
    -------
    mass_1 : float
        Mass of star 1. Unit is kg.
    mass_2 : float
        Mass of star 2. Unit is kg.
    
    Notes
    -----
    Note: Masses are returned as a tuple: (m1, m2)
    m1 + m2 = mass_sum
    m1 / m2 = mass_ratio
    => m1 = mass_ratio*m2
       m1 + m2 = m2*(mass_ratio + 1) = mass_sum
       m2 = mass_sum / (mass_ratio + 1)


    """
    mass_2 = mass_sum / (mass_ratio + 1.0)
    mass_1 = mass_ratio * mass_2
    return (mass_1, mass_2)


def calc_flux_rad_ratio_from_light(light_oc, light_tr, light_ref=1.0):
    """Calculate ratio of the radiative fluxes of the smaller star to greater star
    from their light levels  during occultation and transit.
    Assumes L \propto r^2, no limb darkening.
    
    Parameters
    ----------
    light_oc : float
        Light level during occultation (primary eclipse). Light level is the average
        of the flat region at the bottom of one of the eclipse minima. The event is
        a total eclipse of the star with the smaller radius.
    light_tr : float
        Light level during transit (secondary eclipse). Light level is the average
        of the rounded region at the bottom of one of the eclipse minima. The event is
        an annular eclipse of the star with the larger radius.
    light_ref : {1.0}, float, optional
        Reference for light levels outside of eclipse. Typically total system light is normalized to 1.0.
        
    Returns
    -------
    flux_rad_ratio : float
        Ratio of radiative fluxes of smaller star to greater star.
        flux_rad_ratio = flux_rad_s / flux_rad_g

    Notes
    -----
    flux_rad_ratio = flux_rad_s / flux_rad_g
                   = (light_ref - light_oc) / (light_ref - light_tr)
    From equation 7.10 in section 7.3 of [1]_
    
    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics
    
    """
    flux_rad_ratio = (light_ref - light_oc) / (light_ref - light_tr)
    return flux_rad_ratio


def calc_teff_ratio_from_flux_rad_ratio(flux_rad_ratio):
    """Calculate ratio of the effective temperatures of the smaller star to greater star
    from the ratios of the radiative fluxes.
    
    Parameters
    ----------
    flux_rad_ratio : float
        Ratio of radiative fluxes of smaller star to greater star.
        flux_rad_ratio = flux_rad_s / flux_rad_g
        
    Returns
    -------
    teff_ratio : float
        Ratio of effective temperatures of smaller-sized star to greater-sized star.
        teff_ratio = teff_s / teff_g

    Notes
    -----
    teff_ratio = teff_s / teff_g
               = flux_rad_ratio**0.25 = (flux_rad_s / flux_rad_g)**0.25
    From equation 7.10, 7.11 in section 7.3 of [1]_
    
    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics
    
    """
    teff_ratio = flux_rad_ratio**0.25
    return teff_ratio


def calc_lum_ratio_from_radii_teff_ratios(radii_ratio, teff_ratio):
    """Calculate ratio of the luminosities of the smaller star to greater star
    from the ratios of their radii and effective temperature.
    
    Parameters
    ----------
    radii_ratio : float
        Ratio of radii of smaller-sized star to greater-sized star.
        radii_ratio = radius_s / radius_g
    teff_ratio : float
        Ratio of effective temperatures of smaller-sized star to greater-sized star.
        teff_ratio = teff_s / teff_g
        
    Returns
    -------
    lum_ratio : float
        Ratio of luminosities of smaller-sized star to greater-sized star.
        lum_ratio = lum_s / lum_g

    Notes
    -----
    lum_ratio = lum_s / lum_g
              = (4 * pi * radius_s**2 * sigma * teff_s**4) / (4 * pi * radius_g**2 * sigma * teff_g**4)
    radii_ratio = radius_s / radius_g
    teff_ratio = teff_s / teff_g
    => lum_ratio = radii_ratio**2 * teff_ratio**4
    From equation 3.17 in section 3.4 of [1]_
    
    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics
    
    """
    lum_ratio = radii_ratio**2.0 * teff_ratio**4.0
    return lum_ratio
