from distutils.core import setup
from distutils.extension import Extension


# silo_root = '/usr/gapps/silo/4.9/chaos_5_x86_64_ib'
# hdf5_root = '/usr/gapps/silo/hdf5/1.8.16/chaos_5_x86_64_ib_icc'
#silo_root = '/usr/gapps/silo/4.10.3/toss_3_x86_64-gcc-4.9'
#hdf5_root = '/usr/gapps/silo/hdf5/1.10.2/toss_3_x86_64'
silo_root = '/Users/zli2/Codes/silo/silo-4.10'
hdf5_root = '/usr/local/Cellar/hdf5/1.10.5_1/'


setup(name='geos_post',
      version='2.2.0',
      description='GEOS post processing package',
      author='Chris Sherman',
      author_email='sherman27@llnl.gov',
      ext_modules=[Extension('silo_parser',
                             sources=['silo_parser.cpp'],
                             extra_compile_args=['-I%s/include' % (silo_root)],
                             library_dirs=['%s/lib' % (silo_root),
                                           '%s/lib' % (hdf5_root)],
                             runtime_library_dirs=['%s/lib' % (silo_root),
                                                   '%s/lib' % (hdf5_root)],
                             libraries=['hdf5','siloh5'])])

