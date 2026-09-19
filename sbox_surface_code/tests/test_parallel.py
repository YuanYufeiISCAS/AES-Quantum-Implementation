"""Regression and corruption checks for concurrent primitive schedules."""
from copy import deepcopy
import unittest

from sbox_compile.model import CASES, ROOT, read
from sbox_compile.parallel import paper_cases, reproduce, schedule
from sbox_compile.parallel_model import make_model
from sbox_compile.parallel_verify import audit
from sbox_compile.verify import verify


class PaperSchedules(unittest.TestCase):
    def test_all_paper_results_and_replay(self):
        expected = dict(zip(CASES, (495, 478, 472, 512, 520, 459, 440, 447, 428, 439)))
        count = 0
        for entry, circuit, report in paper_cases():
            with self.subTest(case=entry['file']):
                latency = 412 if entry['role'] == 'refinement' else expected[entry['case']]
                self.assertEqual(report['cost']['latency'], latency)
                replayed = reproduce(circuit)
                self.assertEqual(replayed['cost'], circuit['cost'])
                self.assertEqual(replayed['source'], circuit['source'])
                self.assertGreater(report['timing']['H_CNOT_overlap_cycles'], 0)
                self.assertEqual(circuit['space_time']['reserved_patch_cycles'],
                    circuit['space_time']['reserved_module_patches'] * latency)
                self.assertNotIn('/data_600G/', str(circuit))
                count += 1
        self.assertEqual(count, 11)

    def test_small_fresh_portfolio(self):
        stored = read(ROOT / 'paper_results/cstar_4row_refinement.json.gz')
        source = deepcopy(stored['source'])
        result, report = schedule(source, starts=4)
        self.assertEqual(source, stored['source'])
        self.assertLessEqual(result['cost']['latency'], source['cost']['latency'])
        self.assertEqual(report['starts_per_variant'], 4)
        self.assertTrue(verify(result)['passed'])
        self.assertEqual({k:v for k,v in result['cost'].items() if k.startswith('N_')},
                         {k:v for k,v in source['cost'].items() if k.startswith('N_')})

    def test_split_model(self):
        circuit = read(ROOT / 'paper_results/cstar_4row.json.gz')
        model = make_model(circuit['source'], split_linear=True)
        self.assertFalse(any(n['kind'] == 'linear_cnot_layer' for n in model['nodes']))
        self.assertEqual(sum(n['kind'] == 'linear_cnot' for n in model['nodes']),
                         circuit['cost']['N_CNOT_linear'])


class CorruptSchedules(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.circuit = read(ROOT / 'paper_results/cstar_4row_refinement.json.gz')
        cls.model = make_model(cls.circuit['source'], split_linear=cls.circuit['split_linear'])
        cls.nodes = {n['id']: n for n in cls.model['nodes']}

    def reject_timing(self, change):
        timeline = deepcopy(self.circuit['schedule'])
        change(timeline)
        with self.assertRaises((AssertionError, ValueError)):
            audit(self.circuit['source'], timeline, split_linear=self.circuit['split_linear'])

    def test_missing_correction(self):
        self.reject_timing(lambda t: t['events'].remove(next(
            e for e in t['events'] if self.nodes[e['id']]['kind'] == 'conditional_cnot')))

    def test_missing_initial_port_reset(self):
        self.reject_timing(lambda t: t['events'].remove(next(
            e for e in t['events'] if e['id'] == 'boundary:port_reset')))

    def test_missing_reused_port_reset(self):
        self.reject_timing(lambda t: t['events'].remove(next(
            e for e in t['events'] if self.nodes[e['id']]['kind'] == 'port_reset')))

    def test_duplicate_event(self):
        self.reject_timing(lambda t: t['events'].append(deepcopy(t['events'][0])))

    def test_wrong_duration(self):
        def change(t):
            e = next(e for e in t['events'] if self.nodes[e['id']]['kind'] == 'h')
            e['end'] -= 1
        self.reject_timing(change)

    def test_wrong_makespan(self):
        self.reject_timing(lambda t: t.update(latency=t['latency'] - 1))

    def test_early_feedback(self):
        def change(t):
            consumer = next(iter(self.model['outcomes_by_consumer']))
            e = next(e for e in t['events'] if e['id'] == consumer)
            e.update(start=0, end=self.nodes[consumer]['duration'])
        self.reject_timing(change)

    def test_h_and_cnot_footprint_conflict(self):
        def change(t):
            by_id = {e['id']: e for e in t['events']}
            h = next(n for n in self.model['nodes'] if n['kind'] == 'h')
            footprint = set(map(tuple, h['footprint']))
            cnot = next(n for n in self.model['nodes'] if n.get('cnot_count', 0)
                        and footprint.intersection(map(tuple, n['footprint'])))
            e = by_id[h['id']]
            e.update(start=by_id[cnot['id']]['start'], end=by_id[cnot['id']]['start']+h['duration'])
        self.reject_timing(change)

    def test_zero_duration_order(self):
        def change(t):
            by_id = {e['id']: e for e in t['events']}
            a, b = next((a,b) for a,b in self.model['dependencies']
                        if by_id[a]['start'] == by_id[a]['end'] == by_id[b]['start'])
            i, j = t['events'].index(by_id[a]), t['events'].index(by_id[b])
            t['events'][i], t['events'][j] = t['events'][j], t['events'][i]
        self.reject_timing(change)

    def test_wrong_reported_cost(self):
        bad = deepcopy(self.circuit)
        bad['cost']['latency'] -= 1
        self.assertFalse(verify(bad)['passed'])

    def test_incomplete_source_correction(self):
        bad = deepcopy(self.circuit)
        op = next(o for o in bad['source']['operations'] if o['kind'] == 'conditional_cnot')
        bad['source']['operations'].remove(op)
        self.assertFalse(verify(bad)['passed'])

    def test_source_h_region_cannot_be_shrunk(self):
        bad = deepcopy(self.circuit)
        op = next(o for o in bad['source']['operations'] if o['kind'] == 'h')
        op['footprint']['reserved_patch_coords'].pop()
        self.assertFalse(verify(bad)['passed'])


if __name__ == '__main__':
    unittest.main()
