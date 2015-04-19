#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pytest style tests for binstarsolver/utils.py.

Notes
-----
Test that calculations are correct using:
- example from section 7.3, page 261 of [1]_
- examples 7.3.1, 7.3.2 of [2]_

References
----------
.. [1] Budding, 2007, Introduction to Astronomical Photometry
.. [2] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

"""


# Import standard packages.
from __future__ import absolute_import, division, print_function
import sys
sys.path.insert(0, '.') # Test the code in this repository.
# Import installed packages.
import numpy as np
# Import local packages.
import binstarsolver as bss


def test_calc_flux_intg_ratio_from_mags(
    mag_1=9.6, mag_2=6.3,
    flux_intg_ratio=0.0478630092323):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_flux_intg_ratio_from_mags(
            mag_1=mag_1, mag_2=mag_2),
        flux_intg_ratio)
    return None


def test_calc_fluxes_intg_rel_from_light(
    light_oc=0.898, light_ref=1.0,
    flux_intg_rel_s=0.102, flux_intg_rel_g=0.898):
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(
        bss.utils.calc_fluxes_intg_rel_from_light(
            light_oc=light_oc, light_ref=light_ref),
        (flux_intg_rel_s, flux_intg_rel_g)).all()
    return None


def test_calc_phase_orb_from_time_period(
    time_event=12.3, period=360.0, time_mideclipse=0.0,
    phase_orb=np.deg2rad(12.3)):
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(
        bss.utils.calc_phase_orb_from_time_period(
            time_event=time_event, period=period, time_mideclipse=time_mideclipse),
        phase_orb)
    return None


def test_calc_sep_proj_from_incl_phase(
    incl=1.5514042883817927, phase_orb=0.21467549799530256,
    sep_proj=0.21387118950583997):
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(
        bss.utils.calc_sep_proj_from_incl_phase(
            incl=incl, phase_orb=phase_orb),
        sep_proj)
    return None


def test_calc_radii_ratio_from_light(
    light_oc=0.898, light_tr=0.739, light_ref=1.0,
    radii_ratio_lt=0.53911583146179209):
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(
        bss.utils.calc_radii_ratio_from_light(
            light_oc=light_oc, light_tr=light_tr, light_ref=light_ref),
        radii_ratio_lt)
    return None


def test_calc_radii_sep_from_seps(
    sep_proj_ext=0.213871189506, sep_proj_int=0.0640431640294,
    radius_sep_s=0.0749140127382, radius_sep_g=0.138957176768):
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(
        bss.utils.calc_radii_sep_from_seps(
            sep_proj_ext=sep_proj_ext, sep_proj_int=sep_proj_int),
        (radius_sep_s, radius_sep_g)).all()
    return None


def test_calc_radii_ratio_from_rads(
    radius_sep_s=0.0749140127382, radius_sep_g=0.138957176768,
    radii_ratio=0.53911582316839601):
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(
        bss.utils.calc_radii_ratio_from_rads(
            radius_sep_s=radius_sep_s, radius_sep_g=radius_sep_g),
        radii_ratio)
    return None


def test_calc_incl_from_radii_ratios_phase_incl(
    radii_ratio_lt=0.53911583146179209, phase_orb_ext=0.21467549799530256,
    phase_orb_int=0.061086523819801536, tol=1e-4, maxiter=10, show_plots=False,
    incl=1.5514042883817927):
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(
        bss.utils.calc_incl_from_radii_ratios_phase_incl(
            radii_ratio_lt=radii_ratio_lt, phase_orb_ext=phase_orb_ext,
            phase_orb_int=phase_orb_int, tol=tol, maxiter=maxiter, show_plots=show_plots),
        incl)
    return None


# Test special cases for test_calc_incl_from_radii_ratios_phase_incl
test_calc_incl_from_radii_ratios_phase_incl(
    radii_ratio_lt=0.4381670461247158, phase_orb_ext=0.0469912,
    phase_orb_int=0.01681132, tol=1e-4, maxiter=10, show_plots=False,
    incl=1.5628010760257987)


def test_calc_semimaj_axis_from_period_velr_incl(
    period=271209600.0, velr=33000.0, incl=1.5708021113113511,
    semimaj_axis=1424423498981.198):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_semimaj_axis_from_period_velr_incl(
            period=period, velr=velr, incl=incl),
        semimaj_axis)
    return None


def test_calc_sep_from_semimaj_axes(
    axis_1=1424423498981.198, axis_2=133809480209.56334,
    sep=1558232979214.5923):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_sep_from_semimaj_axes(
            axis_1=axis_1, axis_2=axis_2),
        sep)
    return None


def test_calc_radius_from_radius_sep(
    radius_sep=0.16388077362590259, sep=1558232979214.5923,
    radius=255364426123.08237):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_radius_from_radius_sep(
            radius_sep=radius_sep, sep=sep),
        radius)
    return None


def test_calc_radius_from_velrs_times(
    velr_1=33000.0, velr_2=3100.0, time_1=-7126920.0, time_2=7084800.0,
    radius=256521546000.0):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_radius_from_velrs_times(
            velr_1=velr_1, velr_2=velr_2, time_1=time_1, time_2=time_2),
        radius)
    return None


def test_calc_mass_ratio_from_velrs(
    velr_1=33000.0, velr_2=3100.0,
    mass_ratio=0.09393939393939393):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_mass_ratio_from_velrs(
            velr_1=velr_1, velr_2=velr_2),
        mass_ratio)
    return None


def test_calc_mass_sum_from_period_velrs_incl(
    period=271209600.0, velr_1=33000.0, velr_2=3100.0, incl=1.5708021113113511,
    mass_sum=3.0427831666779509e+31):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_mass_sum_from_period_velrs_incl(
            period=period, velr_1=velr_1, velr_2=velr_2, incl=incl),
        mass_sum)
    return None


def test_calc_masses_from_ratio_sum(
    mass_ratio=0.09393939393939393, mass_sum=3.0427831666779509e+31,
    mass_1=2.6129162927151373e+30, mass_2=2.7814915374064368e+31):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_masses_from_ratio_sum(
            mass_ratio=mass_ratio, mass_sum=mass_sum),
        (mass_1, mass_2)).all()
    return None


def test_calc_flux_rad_ratio_from_light(
    light_oc=0.04786300923226385, light_tr=0.7585775750291839, light_ref=1.0,
    flux_rad_ratio=3.94386308928358):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_flux_rad_ratio_from_light(
            light_oc=light_oc, light_tr=light_tr, light_ref=light_ref),
        flux_rad_ratio)
    return None


def test_calc_teff_ratio_from_flux_rad_ratio(
    flux_rad_ratio=3.94386308928358,
    teff_ratio=1.409225384334092):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_teff_ratio_from_flux_rad_ratio(
            flux_rad_ratio=flux_rad_ratio),
        teff_ratio)
    return None


def test_calc_lum_ratio_from_radii_teff_ratios(
    radii_ratio=1.06/2.2, teff_ratio=5940.0/9800.0,
    lum_ratio=0.03133342331313779):
    """Test that calculations are correct using quanties for two main-sequence stars,
    type A0 and G0, from Appendix G of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(
        bss.utils.calc_lum_ratio_from_radii_teff_ratios(
            radii_ratio=radii_ratio, teff_ratio=teff_ratio),
        lum_ratio)
    return None
