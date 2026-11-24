# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""What a downstream ``mypy --strict`` run sees of Bio.PDB.kdtrees.KDTree.

Type-checked with the mypy.ini beside it, never run. KDTree takes its
arguments by position or by keyword.
"""

import numpy as np

from Bio.PDB.kdtrees import KDTree

coordinates = np.zeros((4, 3))
KDTree(coordinates, 2)
KDTree(coordinates, bucket_size=2)
KDTree(coordinates=coordinates, bucket_size=2)
