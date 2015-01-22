"""Tests for binstarsolver/utils.py.

Test that calculations are correct using:
- example from section 7.3, page 261 of [1]_
- examples 7.3.1, 7.3.2 of [2]_

Notes
-----
Tests are executed using pytest.

References
----------
.. [1] Budding, 2007, Introduction to Astronomical Photometry
.. [2] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

"""


from __future__ import absolute_import, division, print_function
import numpy as np
import binstarsolver as bss


def test_calc_flux_intg_ratio_from_mags(mag_1=9.6, mag_2=6.3, flux_intg_ratio=0.0478630092323):
    """Test that calculations are correct using examples 7.3.1, 7.3.2 of [1]_

    References
    ----------
    .. [1] Carroll and Ostlie, 2007, An Introduction to Modern Astrophysics

    """
    assert np.isclose(bss.utils.calc_flux_intg_ratio_from_mags(mag_1=mag_1, mag_2=mag_2),
                      flux_intg_ratio)
    return None


def test_calc_flux_intg_rel_g_from_light(light_oc=0.898, light_ref=1.0, flux_intg_rel_g=0.898):
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_flux_intg_rel_g_from_light(light_oc=light_oc,
                                                                light_ref=light_ref),
                      flux_intg_rel_g)
    return None


def test_calc_flux_intg_rel_s_from_light(light_oc=0.898, light_ref=1.0, flux_intg_rel_s=0.102): 
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

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
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

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
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

    References
    ----------
    .. [1] Budding, 2007, Introduction to Astronomical Photometry

    """
    assert np.isclose(bss.utils.calc_sep_proj_from_incl_phase(incl=incl, phase_orb=phase_orb),
                      sep_proj)
    return None


def test_calc_radii_ratio_from_light(light_oc=0.898, light_tr=0.739,
                                     light_ref=1.0, radii_ratio_lt=0.53911583146179209):
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

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
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

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
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

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
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

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
    """Test that calculations are correct using example from section 7.3, page 261 of [1]_

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


def test_calc_semimaj_axis_from_period_velr_incl(period, velr, incl):
    pass


def test_calc_sep_from_semimaj_axes(axis_1, axis_2):
    pass


def test_calc_radius_from_radius_sep(radius_sep, sep):
    pass


def test_calc_radius_from_velrs_times(velr_1, velr_2, time_1, time_2):
    pass


def test_calc_mass_ratio_from_velrs(velr_1, velr_2):
    pass


def test_calc_mass_sum_from_period_velrs_incl(period, velr_1, velr_2, incl):
    pass


def test_calc_masses_from_ratio_sum(mass_ratio, mass_sum):
    pass


def test_calc_flux_rad_ratio_from_light(light_oc, light_tr, light_ref=1.0):
    pass


def test_calc_teff_ratio_from_flux_rad_ratio(flux_rad_ratio):
    pass
