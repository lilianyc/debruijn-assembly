import os
import sys
# Path of debruijn package
#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
#                                                '..')))
# Path for core
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '../debruijn')))
import debruijn
#import debruijn_comp