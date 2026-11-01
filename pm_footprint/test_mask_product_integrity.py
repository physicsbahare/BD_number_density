"""Offline provenance and raster-integrity tests; no FITS or network needed."""
import copy
import hashlib
import tempfile
import unittest
from pathlib import Path

from .mask_product_integrity import (
    OFFICIAL_ARCHIVE_SHA256, WCS_CONVENTION, validate_mask_product
)


def fixture(sample_bytes=b'x_lambert_rad,y_lambert_rad\n0.0,0.0\n'):
    tile = {
        'native_area_deg2': .0007,
        'wcs_center_area_ratio_to_nominal_30mas': 1.,
        'wcs_pc_off_diagonal': [0., 0.],
        'grid_0.3_area_ratio_to_native': .993,
        'grid_0.6_area_ratio_to_native': .999,
    }
    return {
        'wcs_convention': WCS_CONVENTION,
        'source_archive_sha256': OFFICIAL_ARCHIVE_SHA256,
        'tiles': [copy.deepcopy(tile) for _ in range(20)],
        'sum_tile_native_masked_area_deg2_NOT_UNION': .014,
        'grid_results': {
            '0.3': {'union_lower_bound_deg2': .0126, 'union_upper_bound_deg2': .0128},
            '0.6': {'union_lower_bound_deg2': .0127, 'union_upper_bound_deg2': .0129},
        },
        'weighted_mc_sample': {
            'source_grid_arcsec': .6,
            'n': 100,
            'sha256': hashlib.sha256(sample_bytes).hexdigest(),
        }
    }


class TestMaskProductIntegrity(unittest.TestCase):
    def test_fresh_product_and_matching_sample(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'points.csv'
            path.write_bytes(b'x_lambert_rad,y_lambert_rad\n0.0,0.0\n')
            self.assertTrue(validate_mask_product(fixture(), path, 100))

    def test_legacy_product_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Stale/unverified'):
            validate_mask_product({'grids': {'0.6arcsec': {}}})

    def test_wrong_source_archive_rejected(self):
        s = fixture(); s['source_archive_sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'archive SHA256'):
            validate_mask_product(s)

    def test_missing_tile_wcs_metadata_rejected(self):
        s = fixture(); del s['tiles'][0]['wcs_pc_off_diagonal']
        with self.assertRaisesRegex(ValueError, 'PC off-diagonal'):
            validate_mask_product(s)

    def test_area_sum_mismatch_rejected(self):
        s = fixture(); s['sum_tile_native_masked_area_deg2_NOT_UNION'] = .01
        with self.assertRaisesRegex(ValueError, 'does not conserve'):
            validate_mask_product(s)

    def test_union_larger_than_tile_sum_rejected(self):
        s = fixture(); s['grid_results']['0.6']['union_upper_bound_deg2'] = .02
        with self.assertRaisesRegex(ValueError, 'union area'):
            validate_mask_product(s)

    def test_corrupted_sample_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'points.csv'
            path.write_text('wrong contents')
            with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
                validate_mask_product(fixture(), path, 100)

    def test_wrong_sample_count_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'points.csv'
            path.write_bytes(b'x_lambert_rad,y_lambert_rad\n0.0,0.0\n')
            with self.assertRaisesRegex(ValueError, 'row count mismatch'):
                validate_mask_product(fixture(), path, 99)

    def test_missing_sample_hash_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'points.csv'
            path.write_bytes(b'x_lambert_rad,y_lambert_rad\n0.0,0.0\n')
            s = fixture(); del s['weighted_mc_sample']['sha256']
            with self.assertRaisesRegex(ValueError, 'missing SHA256'):
                validate_mask_product(s, path, 100)


if __name__ == '__main__':
    unittest.main()
