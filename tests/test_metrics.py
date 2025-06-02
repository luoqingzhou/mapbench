import unittest
import scanpy as sc
from mapbench.metrics import metric_template

class TestMetrics(unittest.TestCase):
    def test_template_metric(self):
        ref = sc.AnnData()
        query = sc.AnnData()
        score = metric_template(ref, query, 'batch', 'label', 'predicted_label')
        self.assertIsInstance(score, float)

if __name__ == '__main__':
    unittest.main()
