

import re, mmap, logging, sys
from VariantCallFixer.Functions import openVCF

LINEWIDTH = 80
LOGGER = logging.Logger(__name__, level=logging.FATAL)