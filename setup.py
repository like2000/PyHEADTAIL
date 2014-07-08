'''
@author: Kevin Li, Danilo Quartullo

'''

import numpy as np

import os
import sys
import subprocess
import cython_gsl

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

# Remember that you have to launch this setup.py script from console with the exact following syntax "python setup.py cleanall build_ext --inplace"

args = sys.argv[1:]

# Make a `cleanall` rule to get rid of previously created cython files

if "cleanall" in args:
    print "Deleting cython files..."
    if "lin" in sys.platform:
        subprocess.Popen("rm -rf build", shell = True, executable = "/bin/bash")
        subprocess.Popen("rm -rf *.c", shell = True, executable = "/bin/bash")
        subprocess.Popen("rm -rf *.so", shell = True, executable = "/bin/bash")
        sys.argv[1] = "clean"
    elif "win" in sys.platform:
        os.system('rd /s/q '+ os.getcwd() +'\\build')
        os.system('del /s/q '+ os.getcwd() +'\\*.c')
        os.system('del /s/q '+ os.getcwd() +'\\*.html')
        os.system('del /s/q '+ os.getcwd() +'\\*.pyd')
        os.system('del /s/q '+ os.getcwd() +'\\*.h5')
        os.system('del /s/q '+ os.getcwd() +'\\*.h5part')
        sys.argv[1] = "clean"
    else:
        print "You have not a Windows or Linux operating system. Aborting..."
        sys.exit()

# We want to always use build_ext --inplace
if args.count("build_ext") > 0 and args.count("--inplace") == 0:
    sys.argv.insert(sys.argv.index("build_ext") + 1, "--inplace")

# Set up extension and build
cy_ext = [
        
        ]

cy_ext_options = {"compiler_directives": {"profile": True}, "annotate": True}

setup(cmdclass={'build_ext': build_ext},
      ext_modules=cythonize(cy_ext, **cy_ext_options),)


