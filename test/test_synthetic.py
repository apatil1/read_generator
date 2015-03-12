import unittest
from cStringIO import StringIO

from .. import synthetic

class Test_Tools(unittest.TestCase):
    
    def test_empty_file(self):
        self.assertRaises(IOError, parser.parse_tool_names, StringIO(""))
    
 
    
    def test_read_size(self):
        
