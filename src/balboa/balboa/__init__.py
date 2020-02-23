import os
from .cffi_helpers import get_lib_handle


_this_path = os.path.dirname(os.path.realpath(__file__))

_library_dir = os.getenv('BALBOA_LIBRARY_DIR')
if _library_dir is None:
    _library_dir = os.path.join(_this_path, 'lib')

_include_dir = os.getenv('BALBOA_INCLUDE_DIR')
if _include_dir is None:
    _include_dir = os.path.join(_this_path, 'include')

_lib = get_lib_handle(
    ['-DBALBOA_API=', '-DBALBOA_NOINCLUDE'],
    'balboa.h',
    'balboa',
    _library_dir,
    _include_dir
)

new_context = _lib.balboa_new_context
free_context = _lib.balboa_free_context
set_basis = _lib.balboa_set_basis
get_buffer_len = _lib.balboa_get_buffer_len
get_ao = _lib.balboa_get_ao
get_num_aos = _lib.balboa_get_num_aos
get_ao_center = _lib.balboa_get_ao_center
get_geo_offset = _lib.balboa_get_geo_offset
