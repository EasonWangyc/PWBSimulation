import unittest

import numpy as np

from pwb_core import (
    PWBParameters,
    centerline_to_segments,
    generate_pwb_path_3,
    generate_radius_profile_3,
    normalized_path_position,
)


class ComplexGeometryTests(unittest.TestCase):
    def test_complex_path_is_planar_and_reaches_expected_endpoints(self):
        params = PWBParameters()
        path = generate_pwb_path_3(params)

        self.assertEqual(path.shape[1], 3)
        self.assertGreater(len(path), params.curve_points)
        np.testing.assert_allclose(path[0], [0.0, 0.0, 0.0], atol=1e-15)
        np.testing.assert_allclose(path[-1], [params.L, 0.0, -params.h], atol=1e-12)
        np.testing.assert_allclose(path[:, 1], 0.0, atol=1e-15)
        self.assertTrue(np.all(np.diff(path[:, 0]) >= -1e-12))
        self.assertGreater(np.max(path[:, 2]), 0.0)

    def test_bend_lift_controls_actual_arch_height(self):
        params = PWBParameters()
        path = generate_pwb_path_3(params)
        expected_arch_height = params.bend_lift * params.h

        self.assertGreater(np.max(path[:, 2]), 0.8 * expected_arch_height)

    def test_normalized_path_position_runs_from_zero_to_one(self):
        params = PWBParameters()
        path = generate_pwb_path_3(params)
        s = normalized_path_position(path)

        self.assertEqual(len(s), len(path))
        self.assertAlmostEqual(s[0], 0.0)
        self.assertAlmostEqual(s[-1], 1.0)
        self.assertTrue(np.all(np.diff(s) >= -1e-15))

    def test_complex_path_has_smooth_tangent_after_input_straight(self):
        params = PWBParameters()
        path = generate_pwb_path_3(params)
        joint_idx = params.curve_points - 1

        before = path[joint_idx] - path[joint_idx - 1]
        after = path[joint_idx + 1] - path[joint_idx]
        before = before / np.linalg.norm(before)
        after = after / np.linalg.norm(after)
        angle_deg = np.degrees(np.arccos(np.clip(np.dot(before, after), -1.0, 1.0)))

        self.assertLess(angle_deg, 2.0)

    def test_radius_profile_transitions_from_input_to_core_to_output(self):
        params = PWBParameters()
        path = generate_pwb_path_3(params)
        radii = generate_radius_profile_3(params, path)
        s = normalized_path_position(path)

        self.assertEqual(len(radii), len(path))
        self.assertAlmostEqual(radii[0], params.r1)
        self.assertAlmostEqual(radii[-1], params.r2)
        core_idx = int(np.argmin(np.abs(s - 0.5)))
        self.assertAlmostEqual(radii[core_idx], params.r, delta=0.05e-6)
        self.assertTrue(np.all(radii > 0.0))

    def test_centerline_segments_follow_adjacent_points(self):
        params = PWBParameters()
        path = generate_pwb_path_3(params)
        radii = generate_radius_profile_3(params, path)
        segments = centerline_to_segments(path, radii)

        self.assertEqual(len(segments), len(path) - 1)
        first = segments[0]
        np.testing.assert_allclose(first["center"], (path[0] + path[1]) / 2, atol=1e-15)
        self.assertGreater(first["length"], 0.0)
        self.assertGreater(first["radius"], 0.0)
        self.assertAlmostEqual(np.linalg.norm(first["direction"]), 1.0)


if __name__ == "__main__":
    unittest.main()
