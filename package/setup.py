from setuptools import setup, find_packages

# Get readme content
with open("README.rst") as fd:
    LONG_DESCRIPTION = fd.read()

setup(
    name = "gdee",
    version = "1.0.0",
    description = "Python package necessary to run the Gene Disovery and Enzyme Engineering Platform",
    long_description = LONG_DESCRIPTION,
    long_description_content_type = "text/x-rst",
    url = "https://gdee",
    project_urls={
        "Documentation": "https://gdee/",
        "Source": "https://gdee/",
        "Tracker": "https://gdee/issues/",
    },
    author = "Caio S. Souza",
    author_email = "cssouza@itqb.unl.pt",
    license = "Copyrighted",
    classifiers = [
        "Development Status :: 5 - Production/Stable",

        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",

        "Programming Language :: Python :: 3",

        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
    ],
    keywords = "enzyme mutation gene docking modeling",
    packages = find_packages(),
    package_dir = {"gdee": "gdee"},
    ext_package = "gdee",
    # ext_modules = extensions,
    install_requires = ["numpy>=1.14", "mdanalysis>=0.20", "mpi4py>=3.0", "path>=13.0", "biopython>=1.75", "oddt>=0.7", "openbabel-wheel>=3.1.1.1"],
    package_data={"gdee": ["data/blosum_all.npz"]},
)
