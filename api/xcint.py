"""
Python handle to XCint library.
"""

import os
import subprocess
from cffi import FFI

BUILD_DIR = os.path.abspath(os.path.dirname(__file__))

ffi = FFI()

ffi.cdef(
    subprocess.Popen(
        [
            'cc',
            '-E',
            '-DXCINT_API=',
            '-DXCINT_NOINCLUDE',
            os.path.join(BUILD_DIR, 'include', 'xcint.h')
        ],
        stdout=subprocess.PIPE).communicate()[0].decode('utf-8'))

lib = ffi.dlopen(os.path.join(BUILD_DIR, 'lib', 'libxcint_shared.so'))
