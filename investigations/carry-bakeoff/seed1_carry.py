import sys
import ttm.carry as cr
_orig = cr.CarryParams
cr.CarryParams = lambda: _orig(rebuild_seed=1)
from ttm.cli.sort import main
main(sys.argv[1:])
