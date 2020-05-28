import unittest
from pyspecdata import *
from captured_output import captured_output
from logging import shutdown
class TestLogging(unittest.TestCase):
    def setUp(self):
        pass
    def test_basic_log(self):
        # {{{ give a case that should match
        with captured_output() as (out, err):
            self.mylogger = init_logging("info")
            self.mylogger.info("this is a test")
            shutdown() # need this only so I can run multiple tests
        def this_check(x):
            self.assertRegex(x[0],
                    'test_logging\.py.*root.*test_basic_log')
            self.assertRegex(x[1],
                    'INFO: this is a test')
        out = out.getvalue().split('\n')
        this_check(out)
        log_filename = os.path.join(os.path.expanduser('~'),'pyspecdata.log')
        file_out = []
        with open(log_filename, encoding='utf-8') as fp:
            for line in fp:
                file_out.append(line.rstrip())
        this_check(file_out)
        # }}}
    def test_basic_log_mismatch(self):
        # {{{ and one that should fail
        with captured_output() as (out, err):
            self.mylogger = init_logging("info")
            self.mylogger.info("this is NOT a test")
            shutdown() # need this only so I can run multiple tests
        def this_check(x):
            self.assertRegex(x[0],
                    'test_logging\.py.*root.*test_basic_log')
            self.assertNotRegex(x[1],
                    'INFO: this is a test')
        out = out.getvalue().split('\n')
        this_check(out)
        log_filename = os.path.join(os.path.expanduser('~'),'pyspecdata.log')
        file_out = []
        with open(log_filename, encoding='utf-8') as fp:
            for line in fp:
                file_out.append(line.rstrip())
        this_check(file_out)
        # }}}
if __name__ == '__main__':
    unittest.main()
# This can go inside or outside the `with` block
