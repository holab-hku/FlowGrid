from distutils.core import setup
setup(
  name = 'FlowGrid',         
  packages = ['FlowGrid'],   
  version = 'v1.0.1',      
  license='MIT',        
  description = 'FlowGrid implemented in Scanpy',   
  author = 'Xiunan Fang',                   
  author_email = 'xiunan@hku.hk',      
  url = 'https://github.com/holab_hku/FlowGrid',   
  download_url = 'https://github.com/holab_hku/FlowGrid/archive/v1.0.1.tar.gz',    
  keywords = ['Single Cell RNA-seq', 'Scanpy', 'Clustering'],   
  install_requires=[           
          'numpy',
          'pandas',
          'anndata',
          'natsort',
          'typing',
          'sklearn',
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      
    'Intended Audience :: Developers',     
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3',      
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
)
