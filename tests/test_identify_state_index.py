import unittest
from types import SimpleNamespace
from unittest.mock import Mock

import matplotlib.pyplot as plt
import numpy as np

from fluxonium_classes import Double_Junc_Fluxonium


class IdentifyStateIndexTests(unittest.TestCase):
    def test_returns_best_matching_eigenstate_and_overlap(self):
        device = Double_Junc_Fluxonium.__new__(Double_Junc_Fluxonium)
        eigenstates = np.eye(3, dtype=complex)
        states_to_identify = [
            1j * eigenstates[:, 2],
            0.6 * eigenstates[:, 0] + 0.8j * eigenstates[:, 1],
        ]

        indices, overlaps = device.identify_state_index(
            states_to_identify, eigenstates
        )

        self.assertEqual(indices, [2, 1])
        np.testing.assert_allclose(overlaps[0], [0, 0, 1])
        np.testing.assert_allclose(overlaps[1], [0.6, 0.8, 0])

    def test_identification_uses_best_non_excluded_eigenstate(self):
        device = Double_Junc_Fluxonium.__new__(Double_Junc_Fluxonium)
        eigenstates = np.eye(4, dtype=complex)
        reference_state = (
            0.6 * eigenstates[:, 1] + 0.8 * eigenstates[:, 2]
        )

        indices, overlaps = device.identify_state_index(
            [reference_state], eigenstates, excluded_indices=[2]
        )

        self.assertEqual(indices, [1])
        np.testing.assert_allclose(overlaps[0], [0, 0.6, 0.8, 0])

    def test_identification_rejects_weak_matches(self):
        device = Double_Junc_Fluxonium.__new__(Double_Junc_Fluxonium)
        eigenstates = np.eye(5, dtype=complex)
        reference_state = np.ones(5, dtype=complex) / np.sqrt(5)

        with self.assertRaisesRegex(RuntimeError, "below the required minimum"):
            device.identify_state_index(
                [reference_state], eigenstates, min_overlap=0.5
            )

    def test_full_spectrum_plot_highlights_tracked_states(self):
        device = Double_Junc_Fluxonium.__new__(Double_Junc_Fluxonium)
        device.resonators = [SimpleNamespace(get_freq_eff=lambda: 5.0)]
        device.find_chi_full_over_flux = Mock(
            return_value=(
                [0.0, 0.0],
                [100.0, 200.0],
                [300.0, 400.0],
                [4.0, 3.0],
                np.array([[0.0, 1.0], [0.0, 2.0]]),
            )
        )

        fig, ax = device.plot_spectrum_full_over_flux(0, [2.0, 1.0])
        self.addCleanup(plt.close, fig)

        self.assertEqual(
            [line.get_label() for line in ax.lines[-3:]],
            ["Resonator |0,1>", "Qubit |1,0>", "Qubit-resonator |1,1>"],
        )
        np.testing.assert_allclose(ax.lines[-3].get_ydata(), [5.2, 5.1])
        np.testing.assert_allclose(ax.lines[-2].get_ydata(), [3.0, 4.0])
        np.testing.assert_allclose(ax.lines[-1].get_ydata(), [8.4, 9.3])


if __name__ == "__main__":
    unittest.main()
