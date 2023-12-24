import setuptools

with open("README.md", "r",encoding="utf-8") as fh:
  long_description = fh.read()

setuptools.setup(
  name="ervdetective",
  version="1.0.1",
  author="Zhi-Jian Zhou",
  author_email="zjzhou@hnu.edu.cn",
  description="An efficient pipeline for identification and annotation of endogenous retroviruses (ERVs)",
  keywords="endogenous retroviruses, virus, evolution",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/ZhijianZhou01/ervdetective",
  packages=setuptools.find_packages(),
  install_requires=["biopython>=1.78", "psutil>=5.9.1", "colorama>=0.4.6"],
  include_package_data=True,

  classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)",
        "Operating System :: OS Independent",
    ],
  entry_points={
             'console_scripts': [
                 'ervdetective = ervdetective.main:starts',
             ],
    }
)
