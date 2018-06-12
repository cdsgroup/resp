from setuptools import setup

setup(
    name='resp',
    url='https://github.com/cdsgroup/resp',
    license='BSD license',
    packages=['resp'],
    install_requires=['numpy'],
    package_data={'resp': ['tests/*.py'],
                  },
    zip_safe=False,
    classifiers=[
          "Intended Audience :: Developers",
          "Operating System :: OS Independent",
          "Topic :: Software Development",
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
          "Development Status :: 5 - Production/Stable",
          "License :: OSI Approved :: BSD License"
    ],
)
