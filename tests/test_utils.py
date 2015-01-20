"""Tests for binstarsolver/utils.py.

Test that calculations are correct using example from section 7.3 of [1]_

References
----------
.. [1] Budding, 2007, Introduction to Astronomical Photometry

"""


from __future__ import absolute_import, division, print_function
import sys
import numpy as np
import binstarsolver as bss


def test_calc_flux_intg_rel_g_from_light(light_oc=0.898, light_ref=1.0, flux_intg_rel_g=0.898):
    """Test that calculations are correct using example from section 7.3 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_flux_intg_rel_g_from_light(light_oc=light_oc,
                                                                light_ref=light_ref),
                      flux_intg_rel_g)
    return None


def test_calc_flux_intg_rel_s_from_light(light_oc=0.898, light_ref=1.0, flux_intg_rel_s=0.102): 
    """Test that calculations are correct using example from section 7.3 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_flux_intg_rel_s_from_light(light_oc=light_oc,
                                                                light_ref=light_ref),
                      flux_intg_rel_s)
    return None


def test_calc_phase_orb_from_time_period(time_event=12.3, period=360.0, time_mideclipse=0.0,
                                         phase_orb=np.deg2rad(12.3)):
    """Test that calculations are correct using example from section 7.3 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_phase_orb_from_time_period(time_event=time_event, period=period,
                                                                time_mideclipse=time_mideclipse),
                      phase_orb)
    return None


def test_calc_sep_proj_from_incl_phase(incl=1.5514042883817927,
                                       phase_orb=0.21467549799530256,
                                       sep_proj=0.21387118950583997):
    """Test that calculations are correct using example from section 7.3 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_sep_proj_from_incl_phase(incl=incl, phase_orb=phase_orb),
                      sep_proj)
    return None


def test_calc_radii_ratio_from_light(light_oc=0.898, light_tr=0.739,
                                     light_ref=1.0, radii_ratio_lt=0.53911583146179209):
    """Test that calculations are correct using example from section 7.3 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_radii_ratio_from_light(light_oc=light_oc,
                                                            light_tr=light_tr,
                                                            light_ref=light_ref),
                      radii_ratio_lt)
    return None


def test_calc_radius_sep_g_from_sep(sep_proj_ext=0.213871189506,
                                    sep_proj_int=0.0640431640294,
                                    radius_sep_g=0.138957176768):
    """Test that calculations are correct using example from section 7.3 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_radius_sep_g_from_sep(sep_proj_ext=sep_proj_ext,
                                                           sep_proj_int=sep_proj_int),
                      radius_sep_g)
    return None


def test_calc_radius_sep_s_from_sep(sep_proj_ext=0.213871189506,
                                    sep_proj_int=0.0640431640294,
                                    radius_sep_s=0.0749140127382):
    """Test that calculations are correct using example from section 7.3 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_radius_sep_s_from_sep(sep_proj_ext=sep_proj_ext,
                                                           sep_proj_int=sep_proj_int),
                      radius_sep_s)
    return None


def test_calc_radii_ratio_from_rads(radius_sep_s=0.0749140127382,
                                    radius_sep_g=0.138957176768,
                                    radii_ratio=0.53911582316839601):
    """Test that calculations are correct using example from section 7.3 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_radii_ratio_from_rads(radius_sep_s=radius_sep_s,
                                                           radius_sep_g=radius_sep_g),
                      radii_ratio)
    return None


def test_calc_incl_from_radii_ratios_phase_incl(radii_ratio_lt=0.53911583146179209,
                                                phase_orb_ext=0.21467549799530256,
                                                phase_orb_int=0.061086523819801536,
                                                incl_init=np.deg2rad(85.0),
                                                show_plot=False,
                                                incl=1.5514042883817927):
    """Test that calculations are correct using example from section 7.3 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_incl_from_radii_ratios_phase_incl(radii_ratio_lt=radii_ratio_lt,
                                                                       phase_orb_ext=phase_orb_ext,
                                                                       phase_orb_int=phase_orb_int,
                                                                       incl_init=incl_init,
                                                                       show_plot=show_plot),
                      incl)
    return None
