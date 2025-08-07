import sys
import types

import numpy as np
import h5py
import tables

from conftest import load_module

# ensure optional compiled dependency is stubbed
sys.modules.setdefault('_nnls', types.ModuleType('_nnls'))
# numpy.core.rec was removed in newer numpy, but nddata.hdf5_write expects it
rec = types.ModuleType('rec')
try:
    rec.fromarrays = np.core.records.fromarrays
except Exception:
    rec.fromarrays = np.core.records.fromarrays
sys.modules['numpy.core.rec'] = rec
np.core.rec = rec

# load nddata using helper to avoid requiring full dependencies
core = load_module('core', use_real_pint=True)
nddata = core.nddata


def test_hdf5_layout(tmp_path):
    a = nddata(np.arange(6).reshape(2, 3), ['t', 'f'])
    a.labels({'t': np.linspace(0.0, 1.0, 2), 'f': np.array([10, 20, 30])})
    a.set_units('t', 's')
    a.set_units('f', 'Hz')
    a.set_units('V')
    a.name('test_nd')
    a.meta = {'level1': {'level2': 5, 'level2list': [1, 2, 3]}}

    a.hdf5_write('sample.h5', directory=str(tmp_path))

    with h5py.File(tmp_path / 'sample.h5', 'r') as f:
        assert 'test_nd' in f
        g = f['test_nd']
        assert 'data' in g
        dimlabels = [x[0].decode() for x in g.attrs['dimlabels']]
        assert dimlabels == ['t', 'f']

        axes_group = g['axes']
        for name, coords, unit in [
            ('t', a.getaxis('t'), 's'),
            ('f', a.getaxis('f'), 'Hz'),
        ]:
            assert name in axes_group
            ds = axes_group[name]
            np.testing.assert_allclose(ds['data'], coords)
            assert ds.attrs['axis_coords_units'].decode() == unit

        meta = g['meta']
        assert 'level1' in meta
        assert meta['level1'].attrs['level2'] == 5
        np.testing.assert_array_equal(
            meta['level1'].attrs['level2list']['LISTELEMENTS'], [1, 2, 3]
        )
