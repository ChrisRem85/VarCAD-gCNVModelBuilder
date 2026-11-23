import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

# Behave like pandas 3 (copy-on-write by default: .to_numpy() views are read-only), so code that
# modifies such arrays in place fails here too, not only on servers with a newer pandas.
try:
    pd.set_option("mode.copy_on_write", True)
except (KeyError, ValueError, pd.errors.OptionError):   # option removed once CoW is the only mode
    pass
