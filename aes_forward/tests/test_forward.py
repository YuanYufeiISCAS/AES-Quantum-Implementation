"""Fixed-circuit regressions and rejection of incomplete or changed witnesses."""
from copy import deepcopy
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from forward_replay.common import ROOT,load_components,read,write
from forward_replay.circuit import assemble,CONTRACT
from forward_replay.verify import verify,verify_components
from forward_replay import mixcolumns,transport,schedule,key_geometry as key
from forward_replay.geometry import verify_global


class ForwardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.components = load_components()
        cls.program = assemble(cls.components)

    def test_complete_forward_reproduction(self):
        # None of the retained search entry points may be used in reproduction.
        with patch('sbox_compile.parallel.search',side_effect=AssertionError('search invoked')):
            report = verify(self.program,random_trials=8)
        self.assertTrue(report['passed'])
        self.assertEqual(report['timing']['latency'],5075)
        self.assertEqual(report['timing']['round_durations'],[4,504]+[511]*8+[479])
        self.assertEqual(report['geometry']['reserved_patches'],56918)
        self.assertEqual(report['resources']['linear_cnot'],92292)
        self.assertEqual(report['resources']['explicit_cnot'],158052)
        self.assertEqual(report['resources']['ccz_states'],10080)
        self.assertEqual(report['resources']['h'],49200)
        self.assertEqual(report['functional']['known_ciphertexts'][0],
                         '69c4e0d86a7b0430d8cdb78070b4c55a')
        self.assertTrue(report['functional']['scratch_zero'])
        self.assertTrue(report['functional']['retained_lookahead_checked'])

    def test_saved_circuit_replay(self):
        with tempfile.TemporaryDirectory(prefix='aes_forward_test_') as directory:
            path = Path(directory)/'circuit.json.gz'
            write(path,self.program)
            saved = read(path)
            self.assertEqual(saved,self.program)
            self.assertTrue(verify(saved,random_trials=0)['passed'])
            with self.assertRaises(FileExistsError):
                write(path,self.program)

    def test_no_external_workspace_references(self):
        serialized = str(self.program)
        self.assertNotIn('/data_600G/',serialized)
        self.assertNotIn('/home/',serialized)

    def test_mixcolumns_fixed_interfaces(self):
        for block in range(4):
            stats = mixcolumns.verify(self.components['mixcolumns'][block%2],block)
            self.assertEqual((stats['cnots'],stats['layers'],stats['surface_depth']),(169,14,28))
            self.assertEqual(stats['clean_outputs'],20)

    def test_route_direction_corruption(self):
        circuit = deepcopy(self.components['mixcolumns'][0])
        circuit['layers'][0]['operations'][0]['path'].reverse()
        with self.assertRaises(ValueError):
            mixcolumns.verify(circuit,0)

    def test_free_output_permutation_rejected(self):
        circuit = deepcopy(self.components['mixcolumns'][0])
        circuit['output_permutation'][0:2] = [1,0]
        with self.assertRaises(ValueError):
            mixcolumns.verify(circuit,0)

    def test_deferred_frame_rejected(self):
        circuit = deepcopy(self.components['mixcolumns'][0])
        circuit['h_frame'] = 1
        with self.assertRaises(ValueError):
            mixcolumns.verify(circuit,0)

    def test_missing_shiftrows_transfer(self):
        transport_data = deepcopy(self.components['shiftrows'])
        transport_data['routes'].pop()
        with self.assertRaises(AssertionError):
            transport.verify(transport_data)

    def test_wrong_adapter_slot(self):
        adapter = deepcopy(self.components['adapter'])
        adapter['schedule']['paths'][0]['target_slot'] = 7
        with self.assertRaises(ValueError):
            transport.verify_adapter(adapter)

    def test_wrong_key_word_routes(self):
        routes = deepcopy(self.components['key_words'])
        routes['layers'][0],routes['layers'][1] = routes['layers'][1],routes['layers'][0]
        with self.assertRaises(ValueError):
            key.verify_round_local_routes(self.components['key_interface'],routes)

    def test_missing_lookahead_copy(self):
        components = deepcopy(self.components)
        components['lookahead']['load_groups'].pop()
        components['lookahead']['load_paths'].pop()
        with self.assertRaises(ValueError):
            verify_components(components)

    def test_missing_dependency(self):
        nodes = deepcopy(self.program['macros'])
        next(n for n in nodes if n['id'] == 'ARK1')['predecessors'].pop()
        with self.assertRaises(ValueError):
            schedule.check(nodes,self.components)

    def test_wrong_timestamp(self):
        nodes = deepcopy(self.program['macros'])
        next(n for n in nodes if n['id'] == 'SB10')['start'] -= 1
        with self.assertRaises(ValueError):
            schedule.check(nodes,self.components)

    def test_missing_global_reservation(self):
        records = deepcopy(self.program['reservations'])
        records.pop()
        with self.assertRaises(ValueError):
            verify_global(self.components,self.program['macros'],records)

    def test_contract_is_not_shared_mutable_state(self):
        program = assemble(self.components)
        program['contract']['CNOT_cycles'] = 1
        self.assertEqual(CONTRACT['CNOT_cycles'],2)
        with self.assertRaises(ValueError):
            verify(program,random_trials=0)

    def test_false_gate_count_rejected(self):
        program = deepcopy(self.program)
        program['resources']['conditional_cnot'] = 0
        with self.assertRaises(ValueError):
            verify(program,random_trials=0)

    def test_wrong_sbox_timing_rejected(self):
        program = deepcopy(self.program)
        program['components']['state']['cost']['latency'] -= 1
        with self.assertRaises(ValueError):
            verify(program,random_trials=0)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory(prefix='aes_forward_cli_') as directory:
            marker = Path(directory)/'preserve.txt'
            marker.write_text('unchanged')
            result = subprocess.run([sys.executable,'-B',str(ROOT/'reproduce.py'),
                                     '--output',directory],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertEqual(marker.read_text(),'unchanged')


if __name__ == '__main__':
    unittest.main()
